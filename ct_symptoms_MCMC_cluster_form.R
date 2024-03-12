## MCMC function
MCMC_function <- function(x){
  n_iterations <- x$n_iterations
  parameters <- x$parameters
  vl_growth_type <- x$vl_growth_type
  vl_kinetics_type <- x$vl_kinetics_type
  vl_truncation <- x$vl_truncation
  delay_dist <- x$delay_dist
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
  test_form <- x$test_form
  ignored_pars <- x$ignored_pars
  spline_burnin <- x$spline_burnin
  strain <- x$strain
  max_inf <- x$max_inf
  symp_delay_lim <- x$symp_delay_lim
  day_agg <- x$day_agg
  max_day <- x$max_day
  ignored_spline_par <- x$ignored_spline_par
  cov_start <- x$cov_start
  gene <- x$gene
  
  sourceCpp("will3p_vlmax_all_tstoch_vl.cpp")
  sourceCpp("will3p_inc.cpp")
  source("ct_symptoms_key_functions.R")
  
  weekly_aggregator <- function(daily_array){
    dimnames2 <- seq(3.5, dim(daily_array)[2]-3.5, 7)
    
    if(dim(daily_array)[2] %% 7 != 0) print("non integer number of weeks", quote=F)
    
    weekly_array <- array(dim=c(dim(daily_array)[1], ceiling(dim(daily_array)[2]/7), dim(daily_array)[3]))
    for(j in 1:ncol(weekly_array)) weekly_array[,j,] <- apply(daily_array[,(7*(j-1)+1):(7*j),],c(1,3),sum)
    
    dimnames(weekly_array) <- list(dimnames(daily_array)[[1]], dimnames2, dimnames(daily_array)[[3]])
    
    return(weekly_array)
  }
  
  if(day_agg == "weekly") data_array <- weekly_aggregator(daily_array=data_array_raw)
  else data_array <- data_array_raw
  
  ## prior function
  prior <- function(parameters, prior_mean, prior_sds, form, delay_dist, vk){
    theta_start <- which(names(parameters)=="theta1")
    parameters[theta_start:(length(parameters)-1)] <- parameters[theta_start:(length(parameters)-1)]*parameters['multiplier']
    if(vk==1){
      if(form=="incidence") prior_val <- sum(log(dtruncnorm(x=parameters, a=c(-Inf, 0,    0, 0.05, -Inf, 0,    0, 0,   0,  0, rep(-Inf, length(parameters)-11)), mean=prior_mean, sd=prior_sds)))
      if(form=="peak")      prior_val <- sum(log(dtruncnorm(x=parameters, a=c(-Inf, 0, -Inf, 0, -Inf, 0,    0, 0, -Inf, 0, rep(-Inf, length(parameters)-11)), mean=prior_mean, sd=prior_sds)))
    }
    if(vk==0){
      if(form=="incidence") prior_val <- sum(log(dtruncnorm(x=parameters, a=c(-Inf, 0, -Inf, 0, -Inf, 0, -Inf, 0,    0, 0, -Inf, rep(-Inf, length(parameters)-11)), mean=prior_mean, sd=prior_sds)))
      if(form=="peak")      prior_val <- sum(log(dtruncnorm(x=parameters, a=c(-Inf, 0, -Inf, 0, -Inf, 0, -Inf, 0, -Inf, 0, -Inf, rep(-Inf, length(parameters)-11)), mean=prior_mean, sd=prior_sds)))
    }  
    return(prior_val)
  }
  
  ## likelihood function
  likelihood_function3 <- function(parameters, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg){
    rdisp <- parameters['rdisp']
    
    parameters[1:8] = c(log(3), 0, log(2), 0, 9, 0, 0, 0)
    p_array <- p_array_func(parameters, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, stoch=0.5)
    print(sum(p_array, na.rm=T))
    if(day_agg == "weekly") p_array <- weekly_aggregator(p_array)
    index_start <- which(is.na(p_array[1,,1])==F)[1]
    likelihood <- sum(dnbinom(x=data_array[,index_start:ncol(p_array),], size=1/rdisp, mu=p_array[,index_start:ncol(p_array),], log=T))
    
    return(likelihood)
  }
  
  ## posterior function
  posterior <- function(parameters, prior_mean, prior_sds, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg){
    a <- prior(parameters, prior_mean, prior_sds, form, delay_dist, vk) 
    b <- ifelse(a != -Inf, likelihood_function3(parameters, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg), 0)
    if(is.nan(b)) b <- -Inf
    
    return(a+b)
  }
  
  ## proposal function
  proposal <- function(parameters, iteration, MCMC_sds, proposal_type, cov_matrix, ignored_pars, ignored_spline_par, test_form){
    if(proposal_type == "MH"){
      param_interest <- (iteration-1) %% length(parameters) + 1
      
      new_param <- rnorm(n=1, mean=parameters[param_interest], sd=MCMC_sds[param_interest])
      parameters[param_interest] <- new_param
    }
    
    if(proposal_type == "Cov"){
      output <- c(rmvnorm(n=1, mean=parameters, sigma=cov_matrix))
      names(output) <- names(parameters)
      output[ignored_spline_par] <- parameters[ignored_spline_par]
      if(is.na(ignored_pars)==F) output[ignored_pars] <- parameters[ignored_pars]
      #if(test_form=='empirical') output[c('test1', 'test2')] <- parameters[c('test1', 'test2')]
      
      parameters <- output
    }
    
    return(parameters)
  }
  
  MCMC_output <- matrix(nrow = n_iterations + 1, ncol = length(parameters))
  MCMC_output[1,] <- parameters
  colnames(MCMC_output) <- names(parameters)
  
  vg <- ifelse(vl_growth_type == "lognorm", 1, 0)
  vk <- ifelse(vl_kinetics_type == "non-differentiable", 1, 0)
  vt <- ifelse(vl_truncation == "truncated", 1, 0)
  
  #for(i in 1: 100){
  p <- proc.time()
  current_posterior <- posterior(parameters, prior_mean, prior_sds, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg)
  current_posterior
  print((proc.time()-p)["elapsed"])
  #}  
  
  MCMC_posteriors <- matrix(nrow=n_iterations + 1, ncol=1)
  MCMC_posteriors[1] <- current_posterior
  
  Acceptances <- vector(length=n_iterations)
  time_vec <- vector(length=n_iterations)
  
  for(i in 1:n_iterations){
    if(i > 10000) test_pop <- 1e7
    p <- proc.time()
    cat("\n")
    param_interest <- (i-1) %% length(parameters) + 1
    print(i) 
    print(param_interest)
    
    if(proposal_type=="MH"){
      #adj_params <- length(parameters) - 10
      if((i < spline_burnin & i %% length(parameters) %in% c(1:10)) | i %% length(parameters) == ignored_spline_par){
        print("Rejected", quote=F)
        MCMC_output[i+1,] <- MCMC_output[i,]
        Acceptances[i] <- NA
        MCMC_posteriors[i+1] <- MCMC_posteriors[i]
        next
      }
      if(i >= spline_burnin){
        #adj_params <- length(parameters) - 10
        if(i %% length(parameters) %in% ignored_pars){
          print("Rejected", quote=F)
          MCMC_output[i+1,] <- MCMC_output[i,]
          Acceptances[i] <- NA
          MCMC_posteriors[i+1] <- MCMC_posteriors[i]
          next
        }
      }
    }
    
    if(i >= cov_start) proposal_type <- "Cov"
    
    if(i == cov_start) cov_matrix <- cov(MCMC_output[(i-20000):i,])*(2.38^2)/length(parameters)
      
    if((proposal_type == "Cov") & (i %% 10000 == 0) & (i > 20000)) cov_matrix <- cov(MCMC_output[(i-10000):i,])*(2.38^2)/length(parameters)
    #else adj_params <- length(parameters)
    
    current_parameters <- MCMC_output[i,]
    #names(current_parameters) <- param_names
    
    current_posterior <- MCMC_posteriors[i]
    
    proposed_parameters <- proposal(parameters=current_parameters,iteration=i,MCMC_sds=MCMC_sds,proposal_type=proposal_type, cov_matrix=cov_matrix, ignored_pars=ignored_pars, ignored_spline_par=ignored_spline_par, test_form=test_form)
    print(c("proposed parameters = ", proposed_parameters), quote=F)
    
    saveRDS(proposed_parameters, file="parameters3.rds")
    
    proposed_posterior <- posterior(parameters=proposed_parameters, prior_mean, prior_sds, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg)
    
    #print(c("current posterior = ", current_posterior), quote=FALSE)
    print(c("proposed posterior = ", proposed_posterior), quote=FALSE)
    
    likelihood_ratio <- ifelse(proposed_posterior == -Inf || is.nan(proposed_posterior), 0, exp(proposed_posterior-current_posterior))
    
    if(runif(1) < likelihood_ratio){
      print("Accepted", quote=F)
      MCMC_output[i+1,] <- proposed_parameters
      MCMC_posteriors[i+1] <- proposed_posterior
      Acceptances[i] <- 1
    }
    
    else{
      print("Rejected", quote=F)
      MCMC_output[i+1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      if(runif(1) < 0.01) MCMC_posteriors[i+1] <- posterior(current_parameters, prior_mean, prior_sds, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, day_agg)
      else MCMC_posteriors[i+1] <- MCMC_posteriors[i]
    }
    
    print(c("Last accepted = ", max(which(Acceptances==1))), quote=F)
    print(c("Acceptance rate =", round(sum(Acceptances[1:i], na.rm=T)/i,2)), quote=F)
    
    if(i %% 1000 == 0){
      Acceptances_vec <- vector(length=length(parameters))
      for(j in 1:length(Acceptances_vec)){
        Acceptances_param <- Acceptances[seq(j,i,length(parameters))]
        n_accept_param <- sum(Acceptances_param, na.rm=T)
        n_it_param <- length(Acceptances_param[is.na(Acceptances_param)==FALSE])
        Acceptances_vec[j] <- n_accept_param/n_it_param
      }
      names(Acceptances_vec) <- names(parameters)
      print(Acceptances_vec, quote=F)
      
      saveRDS(list(MCMC_output=MCMC_output, MCMC_posteriors=MCMC_posteriors, Acceptances=Acceptances, cov_matrix=cov_matrix, time_vec=time_vec, x=x), 
              file=paste0(strain,"/",Sys.Date(),"_",strain,"_",format(i,scientific=F),     "_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_growth_type,"_",vl_kinetics_type,"_",vl_truncation,"_",delay_dist,"_",tag,".rds"))
            unlink(paste0(strain,"/",Sys.Date(),"_",strain,"_",format(i-1000,scientific=F),"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_growth_type,"_",vl_kinetics_type,"_",vl_truncation,"_",delay_dist,"_",tag,".rds"))
            unlink(paste0(strain,"/",Sys.Date(),"_",strain,"_",format(i-2000,scientific=F),"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_growth_type,"_",vl_kinetics_type,"_",vl_truncation,"_",delay_dist,"_",tag,".rds"))
            unlink(paste0(strain,"/",Sys.Date(),"_",strain,"_",format(i-3000,scientific=F),"_",gene,"_",form,"_",proposal_type,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_growth_type,"_",vl_kinetics_type,"_",vl_truncation,"_",delay_dist,"_",tag,".rds"))
      
      unlink(paste0(strain,"/",Sys.Date()-1,"_",strain,"_",format(i-1000,scientific=F),"_",gene,"_",form,"_",proposal_type,"_",test_form,"_ignored_",gsub(" ","",toString(ignored_pars)),"_",vl_growth_type,"_",vl_kinetics_type,"_",vl_truncation,"_",delay_dist,"_",tag,".rds"))
    }
    print((proc.time()-p)["elapsed"])
    time_vec[i] <- ((proc.time()-p)["elapsed"])
  }
}

