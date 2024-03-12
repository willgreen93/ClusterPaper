## MCMC function
MCMC_function <- function(x){
  n_iterations <- x$n_iterations
  parameters <- x$parameters
  vl_kinetics_type <- x$vl_kinetics_type
  knots <- x$knots
  population <- x$population
  data_array_raw <- x$data_array
  test_pop <- x$test_pop
  MCMC_sds <- x$MCMC_sds
  prior_mean <- x$prior_mean
  prior_sds <- x$prior_sds
  ncores <- x$ncores
  proposal_type <- x$proposal_type
  cov_matrix <- x$cov_matrix
  form <- x$form
  tag <- x$tag
  strain <- x$strain
  max_inf <- x$max_inf
  symp_delay_lim <- x$symp_delay_lim
  max_day <- x$max_day
  ignored_spline_par <- x$ignored_spline_par
  cov_start <- x$cov_start
  gene <- x$gene
  days <- x$days
  dates <- x$dates
  n_days <- x$n_days
  vls <- x$vls
  n_vls <- x$n_vls
  n_test_days <- x$n_test_days
  test_days <- x$test_days
  n_parameters <- x$n_parameters
  ignored_pars <- x$ignored_pars 
  vg <- x$vg
  prior_cov_final <- x$prior_cov_final
  lower_bound <- x$lower_bound
  
  #Rcpp::sourceCpp("will3p_inc_vl_reduced.cpp")
  Rcpp::sourceCpp("will3p_peak_vl_reduced.cpp")
  Rcpp::sourceCpp("will3p_thresh_vl_reduced2.cpp")
  Rcpp::sourceCpp("will3p_thresh_peak_vl_reduced2.cpp")
  Rcpp::sourceCpp("will3p_thresh_vl_tdist_reduced2.cpp")
  
  source("key_functions_simplified.R")
  
  if(length(dim(data_array_raw))==3) data_array <- weekly_aggregator(daily_array=data_array_raw)
  if(length(dim(data_array_raw))==2) data_array <- data_array_raw
  
  ## posterior function
  posterior <- function(parameters, prior_mean, prior_cov_final, lower_bound, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars){
    a <- prior(parameters, prior_mean, prior_cov_final, lower_bound, form, ignored_pars) 
    
    if(a==-Inf) return(-Inf)
    
    if(length(dim(data_array))==3) b <- likelihood_function3(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim)
    if(length(dim(data_array))==2) b <- likelihood_REACT(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim)
    
    if(is.nan(b)) b <- -Inf
    
    return(c(a+b,b))
  }
  
  ## proposal function
  proposal <- function(parameters, iteration, MCMC_sds, proposal_type, cov_matrix, ignored_pars, ignored_spline_par){
    if(proposal_type == "MH"){
      param_interest <- (iteration-1) %% length(parameters) + 1
      
      new_param <- rnorm(n=1, mean=parameters[param_interest], sd=MCMC_sds[param_interest])
      parameters[param_interest] <- new_param
    }
    
    if(proposal_type == "Cov"){
      output <- c(mvtnorm::rmvnorm(n=1, mean=parameters, sigma=cov_matrix))
      names(output) <- names(parameters)
      output[ignored_spline_par] <- parameters[ignored_spline_par]
      output[ignored_pars] <- parameters[ignored_pars]
      #if(test_form=='empirical') output[c('test1', 'test2')] <- parameters[c('test1', 'test2')]
      
      parameters <- output
    }
    
    return(parameters)
  }
  
  MCMC_output <- matrix(nrow = n_iterations + 1, ncol = n_parameters)
  MCMC_output[1,] <- parameters
  colnames(MCMC_output) <- names(parameters)
  
  vk <- 1
  
  p <- proc.time()
  current_posterior <- posterior(parameters, prior_mean, prior_cov_final, lower_bound, knots, vk=vk, vg=vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars)
  current_posterior
  print((proc.time()-p)["elapsed"])
  
  MCMC_posteriors <- matrix(nrow=n_iterations + 1, ncol=1)
  MCMC_likelihoods <- matrix(nrow=n_iterations + 1, ncol=1)
  
  MCMC_posteriors[1] <- current_posterior[1]
  MCMC_likelihoods[1] <- current_posterior[2]
  
  Acceptances <- vector(length=n_iterations)
  time_vec <- vector(length=n_iterations)
  
  for(i in 1:n_iterations){
    if(i > 10000) test_pop <- 1e7
    p <- proc.time()
    cat("\n")
    param_interest <- (i-1) %% n_parameters + 1
    print(i) 
    print(param_interest)
    
    if(proposal_type=="MH" && (i %% n_parameters %in% c(ignored_spline_par,ignored_pars))){
      print("Rejected", quote=F)
      MCMC_output[i+1,] <- MCMC_output[i,]
      Acceptances[i] <- NA
      MCMC_posteriors[i+1] <- MCMC_posteriors[i]
      MCMC_likelihoods[i+1] <- MCMC_likelihoods[i]
      next
    }
    
    if(i >= cov_start) proposal_type <- "Cov"
    
    if(i == cov_start) cov_matrix <- cov(MCMC_output[(i-10000):i,])*(2.38^2)/n_parameters
      
    if((proposal_type == "Cov") & (i %% 10000 == 0) & (i > 20000)) cov_matrix <- cov(MCMC_output[(i-10000):i,])*(2.38^2)/n_parameters
    
    current_parameters <- MCMC_output[i,]
    current_posterior <- MCMC_posteriors[i]
    
    proposed_parameters <- proposal(parameters=current_parameters,iteration=i,MCMC_sds=MCMC_sds,proposal_type=proposal_type, cov_matrix=cov_matrix, ignored_pars=ignored_pars, ignored_spline_par=ignored_spline_par)
    print(c("proposed parameters = ", proposed_parameters), quote=F)
    
    #saveRDS(proposed_parameters, paste0("parameters_break_",form,".rds"))
    
    #proposed_parameters <- readRDS("parameters_break.rds")
    posterior_calc <- posterior(parameters=proposed_parameters, prior_mean, prior_cov_final, lower_bound, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars)
    
    proposed_posterior <- posterior_calc[1]
    proposed_likelihood <- posterior_calc[2]
    
    #print(c("current posterior = ", current_posterior), quote=FALSE)
    print(c("proposed posterior = ", proposed_posterior), quote=FALSE)
    print(c("proposed likelihood = ", proposed_likelihood), quote=FALSE)
    
    likelihood_ratio <- ifelse(proposed_posterior == -Inf || is.nan(proposed_posterior), 0, exp(proposed_posterior-current_posterior))
    
    if(runif(1) < likelihood_ratio){
      print("Accepted", quote=F)
      MCMC_output[i+1,] <- proposed_parameters
      MCMC_posteriors[i+1] <- proposed_posterior
      MCMC_likelihoods[i+1] <- proposed_likelihood
      Acceptances[i] <- 1
    }
    
    else{
      print("Rejected", quote=F)
      MCMC_output[i+1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      if(runif(1) < 0.01){
        calc_posterior <- posterior(current_parameters, prior_mean, prior_cov_final, lower_bound, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, ignored_pars)
        MCMC_posteriors[i+1] <- calc_posterior[1]
        MCMC_likelihoods[i+1] <- calc_posterior[2]
      } 
      else{
        MCMC_posteriors[i+1] <- MCMC_posteriors[i]
        MCMC_likelihoods[i+1] <- MCMC_likelihoods[i]
      } 
    }
    
    print(c("Last accepted = ", max(which(Acceptances==1))), quote=F)
    print(c("Acceptance rate =", round(sum(Acceptances[1:i], na.rm=T)/i,2)), quote=F)
    
    if(i %% 5000 == 0){
      Acceptances_vec <- vector(length=n_parameters)
      for(j in 1:length(Acceptances_vec)){
        Acceptances_param <- Acceptances[seq(j,i,n_parameters)]
        n_accept_param <- sum(Acceptances_param, na.rm=T)
        n_it_param <- length(Acceptances_param[is.na(Acceptances_param)==FALSE])
        Acceptances_vec[j] <- n_accept_param/n_it_param
      }
      names(Acceptances_vec) <- names(parameters)
      print(Acceptances_vec, quote=F)
      
      folder <- ifelse(length(dim(data_array))==2, "react", "")
      folder2 <- ifelse(length(dim(data_array))==2, "react_", "")
      
      saveRDS(list(MCMC_output=MCMC_output, MCMC_posteriors=MCMC_posteriors, MCMC_likelihoods=MCMC_likelihoods, Acceptances=Acceptances, cov_matrix=cov_matrix, time_vec=time_vec, x=x), 
        file=paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i,scientific=F),"_",Sys.Date(),"_mult2.rds"))
      
      unlink(paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i-5000,scientific=F),"_",Sys.Date(),"_mult2.rds"))
      unlink(paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i-10000,scientific=F),"_",Sys.Date(),"_mult2.rds"))
      unlink(paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i-15000,scientific=F),"_",Sys.Date(),"_mult2.rds"))
      unlink(paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_MH",            "_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i- 5000,scientific=F),"_",Sys.Date(),"_mult2.rds"))
      unlink(paste0(strain,"/",folder,"/",folder2,strain,"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_kinetics_type,"_",tag,"vg=",vg,"_",format(i- 5000,scientific=F),"_",Sys.Date()-1,"_mult2.rds"))
    }
    print((proc.time()-p)["elapsed"])
    time_vec[i] <- ((proc.time()-p)["elapsed"])
  }
}

