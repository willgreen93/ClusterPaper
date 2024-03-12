p_array_func <- function(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5, name=F){
  if(length(dim(data_array))==3){
    max_test_day <- dim(data_array)[3]-1
    max_inf <- 2*(max_test_day+14)+1
  }
  if(length(dim(data_array))==2) max_inf <- 30
    
  inf_lims <- floor(max_inf/2)
  
  ### CHANGE THIS
  vl_vec <- as.numeric(substr(rownames(data_array), 4, nchar(rownames(data_array))))
  vl_range <- vl_vec[c(1,length(vl_vec))]
  vl_vec[1]="negative"
  
  theta_start <- which(names(parameters)=="theta1")
  thetas <- thetas_generator(parameters)[theta_start:(length(parameters)-1)]
  
  ## EXTRACT PARAMS
  inc <- parameters[c('inc1','inc2','inc3')]
  
  test_param <- c(parameters[grepl("test",names(parameters))],test6=1)
  
  #test_delay <- c(test_param, test6=max(0,1-sum(test_param)))
  
  test_delay <- setNames(vector(length=length(test_param)), names(test_param))
  test_delay[1] = test_param[1]
  for (i in 2:length(test_delay)) test_delay[i] <- prod(1-test_param[1:(i-1)])*test_param[i]
  
  pseeds = unlist(dqrng::generateSeedVectors(ncores,1)) 
  
  if(form=="incidence"){
    days <- 0:max_inf
    infecteds <- infecteds_generator(parameters, knots, population, max_day, form) # generate infecteds
    if(length(dim(data_array))==3) inc_period <- setNames(fGarch::psnorm(0:symp_delay_lim, mean=inc[1], sd=inc[2], xi=exp(inc[3]))-fGarch::psnorm(0:symp_delay_lim-1, mean=inc[1], sd=inc[2], xi=exp(inc[3])),0:symp_delay_lim)
    
    p_ct_tau_mat <- ct_func2_cpp_inc(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    dimnames(p_ct_tau_mat) <- list(0:(max_inf-1), rev(vl_vec))
    start_j = 35
    #print(t(round(p_ct_tau_mat),2)
  }
  
  if(form=="peak"){
    days <- -symp_delay_lim:symp_delay_lim
    infecteds <- infecteds_generator(parameters, knots, population, max_day, form) # generate infecteds
    if(length(dim(data_array)==3)) peak_to_symp <- setNames(fGarch::psnorm(days, mean=inc[1], sd=inc[2], xi=exp(inc[3]))-fGarch::psnorm(days-1, mean=inc[1], sd=inc[2], xi=exp(inc[3])),days)
  
    p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, rev(vl_vec))
    
    start_j = 21
    round(t(p_ct_tau_mat),7)
  }
    
  if(form %in% c("thresh","thresh_peak", "thresh_tdist")){
    days <- 0:max_inf
    infecteds <- infecteds_generator(parameters, knots, population, max_day, form) # generate infecteds
    
    p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, rev(vl_vec))
    
    if(form=="thresh") p_ct_tau_mat <- ct_func2_cpp_thresh(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], t=parameters[c(9,10,11)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    if(form=="thresh_peak") p_ct_tau_mat <- ct_func2_cpp_thresh_peak(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], t=parameters[c(9,10,11)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    if(form=="thresh_tdist") p_ct_tau_mat <- ct_func2_cpp_thresh_tdist(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], t=parameters[c(9,10,11)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg=vg, vk, vt=1, stoch=stoch)
    
    dimnames(p_ct_tau_mat) <- list(0:(max_inf-1), rev(vl_vec))
    start_j = 35
    round(t(p_ct_tau_mat),3)
  } 
  
  if(length(dim(data_array))==2) return(p_ct_tau_mat)
  
  else{
    ## CALCULATE p_ct_tau
    p_array <- array(NA, dim=c(length(vl_vec)-1, max_day, max_test_day+1)) 
    
    # generate an array where the row pertains to ct value, the column to calendar day and the slice to days since onset
    if(form %in% c("incidence", "peak")){
      for(j in start_j:max_day) for(k in 0:(dim(p_array)[3]-1)){
        if(form=="incidence") p_array[,j,k+1] <- t(c(infecteds[(j-k):(j-k-symp_delay_lim)]*inc_period[1:(symp_delay_lim+1)]*test_delay[k+1]))%*%p_ct_tau_mat[(k+1):(k+symp_delay_lim+1),(ncol(p_ct_tau_mat)-1):1]
        if(form=="peak")      p_array[,j,k+1] <- t(c(infecteds[(j-symp_delay_lim-k):(j+symp_delay_lim-k)]*rev(peak_to_symp)*test_delay[k+1]))%*%p_ct_tau_mat[(inf_lims+1+symp_delay_lim+k):(inf_lims+1-symp_delay_lim+k),(ncol(p_ct_tau_mat)-1):1]
      }
    }
    
    if(form %in% c("thresh", "thresh_peak", "thresh_tdist")){
      for(j in start_j:max_day){
        p_array[,j,] <- t(c(infecteds[j:(j-max_test_day)]*test_delay) * p_ct_tau_mat[1:(max_test_day+1),(ncol(p_ct_tau_mat)-1):1])
      }
    }
    
    dimnames(p_array) = list(as.numeric(vl_vec[-1]),1:max_day,0:max_test_day)
    
    if(form=="incidence") incubation <- inc_period
    if(form=="peak") incubation <- peak_to_symp
    if(form %in% c("thresh", "thresh_peak", "thresh_tdist")) incubation <- NA
    
    return(list(p_array, p_ct_tau_mat, incubation))  
  }
}

infecteds_generator <- function(parameters, knots, population, max_day, form){
  #thetas <- thetas_generator(parameters)
  thetas <- parameters[grepl("theta", names(parameters))]
  infecteds <- population*exp(stats::splinefun(x=c(0,knots*max_day), y=thetas)(1:(max_day)))*parameters['multiplier']
  if(form=="peak") infecteds <- c(rep(0,14),infecteds)
  
  names(infecteds) = 1:length(infecteds)
  
  return(infecteds)
}

thetas_generator <- function(parameters){
  start <- which(names(parameters) == "theta1")
  end <- length(parameters)-1
  
  #parameters[start:end] <- parameters[start:end]*parameters[length(parameters)]
  #parameters[length(parameters)] <- 1
  
  return(parameters)
}

## prior function
prior <- function(parameters, prior_mean, prior_cov_final, lower_bound, form, ignored_pars){
  spline_pars <- names(parameters)[which(grepl("theta",names(parameters)))]
  #parameters[spline_pars] <- parameters[spline_pars]*parameters['multiplier']
  
  vl_pars <- parameters[1:8][-ignored_pars]
  vl_pars_names <- names(vl_pars)
  
  last_par <- length(parameters)-1
  
  prior_val_vl <- TruncatedNormal::dtmvnorm(x=vl_pars, mu=prior_mean[vl_pars_names], sigma=prior_cov_final[vl_pars_names,vl_pars_names], lb=lower_bound[vl_pars_names], log=T)
  prior_val_non_vl <- log(truncnorm::dtruncnorm(x=parameters[9:last_par], a=lower_bound[9:last_par], mean=prior_mean[9:last_par], sd=sqrt(diag(prior_cov_final)[9:last_par])))
  prior_val <- sum(prior_val_vl, prior_val_non_vl)
  
  return(prior_val)
}

## likelihood function
likelihood_function3 <- function(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim){
  rdisp <- parameters['rdisp']
  p_array <- p_array_func(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5)[[1]]
  
  print(sum(p_array, na.rm=T))
  p_array <- weekly_aggregator(p_array)
  index_start <- which(is.na(p_array[1,,1])==F)[1]
  likelihood <- sum(dnbinom(x=data_array[-1,index_start:ncol(p_array),], size=1/rdisp, mu=p_array[,index_start:ncol(p_array),], log=T))
  likelihood
  return(likelihood)
}

likelihood_REACT <- function(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim){
  p_array <- p_array_func(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5, name=T)
  
  infecteds <- round(infecteds_generator(parameters, knots, population, max_day, form))
  
  N_matrix <- matrix(NA, nrow=nrow(data_array), ncol=ncol(data_array))
  rownames(N_matrix) <- c(rev(rownames(data_array)[-1]),"negative") 
  colnames(N_matrix) <- colnames(data_array)
  
  for(j in 31:ncol(N_matrix)) N_matrix[,j] <- t(p_array) %*% infecteds[j:(j-29)]
  
  N_matrix[nrow(N_matrix),] <- population-colSums(N_matrix[1:(nrow(N_matrix)-1),])
  p_matrix <- N_matrix/population
  
  likelihood <- sum(data_array[nrow(data_array):1,(31:max_day)] * log(p_matrix[,31:max_day]))
  
  return(likelihood)
}

weekly_aggregator <- function(daily_array){
  dimnames2 <- seq(3.5, dim(daily_array)[2]-3.5, 7)
  
  if(dim(daily_array)[2] %% 7 != 0) print("non integer number of weeks", quote=F)
  
  if(length(dim(daily_array))==2){
    weekly_array <- array(dim=c(dim(daily_array)[1], ceiling(dim(daily_array)[2]/7)))
    for(j in 1:ncol(weekly_array)) weekly_array[,j] <- apply(daily_array[,(7*(j-1)+1):(7*j)],c(1),sum)
    dimnames(weekly_array) <- list(dimnames(daily_array)[[1]], dimnames2)
  } 
  
  if(length(dim(daily_array))==3){
    weekly_array <- array(dim=c(dim(daily_array)[1], ceiling(dim(daily_array)[2]/7), dim(daily_array)[3]))
    for(j in 1:ncol(weekly_array)) weekly_array[,j,] <- apply(daily_array[,(7*(j-1)+1):(7*j),],c(1,3),sum)
    dimnames(weekly_array) <- list(dimnames(daily_array)[[1]], dimnames2, dimnames(daily_array)[[3]])
  }
  
  return(weekly_array)
}

output_list <- function(strain, form, corr, flat, recentness){
  files <- list.files(paste0(strain))
  form_filter = grepl(paste0("N_",form,"_Cov"), files) | grepl(paste0("N_",form,"_MH"), files)
  #corr_filter = grepl(paste0("corr",corr), files)
  flat_filter = grepl("flat", files)
  vg_filter = grepl(paste0("vg=",corr), files)
  mult_filter = grepl("mult", files)
  check_filter <- !grepl("checking", files)
  #vg_filter = !grepl("vg", files)
  
  if(flat==0) flat_filter <- !flat_filter
  
  filtered_files <- files[which(form_filter*flat_filter*vg_filter*mult_filter*check_filter==1)]
  
  return(filtered_files)
}

output_finder <- function(strain, form, corr, flat, recentness=1){
  filtered_files <- output_list(strain, form, corr, flat, recentness)
  
  #recentness <- (length(filtered_files))
  
  if(length(filtered_files)==0){
    print("No file") 
    return("No file")
  } 
  info <- file.info(paste0(strain,"/",filtered_files))
  
  times <- info$ctime
  file_name <- paste0(strain, "/", filtered_files[which(times==sort(times, decreasing=T)[recentness])])
  print(file_name, quote=F)
  return(readRDS(file_name))
}
