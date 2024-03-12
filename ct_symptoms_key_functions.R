p_array_func <- function(parameters, knots, vg, vk, vt, delay_dist, max_day, population, data_array, test_pop, ncores, form, test_form, symp_delay_lim, stoch=0.5){
  max_test_day <- dim(data_array)[3]-1
  max_inf <- 2*(max_test_day+14)+1
  inf_lims <- floor(max_inf/2)
  symp_delay_buckets <- symp_delay_lim*2 + 1
  
  ### CHANGE THIS
  vl_vec <- as.numeric(substr(rownames(data_array), 4, nchar(rownames(data_array))))
  vl_range <- vl_vec[c(1,length(vl_vec))]
  vl_vec[1]="negative"
  
  theta_start <- which(names(parameters)=="theta1")
  thetas <- thetas_generator(parameters)[theta_start:(length(parameters)-1)]
  
  ## EXTRACT PARAMS
  inc <- parameters[c('inc1','inc2','inc3')]
  
  test_param <- parameters[grepl("test",names(parameters))]
  
  ## CALCULATE p_tau_s
  if(form=="incidence"){
    days <- 0:max_inf
    infecteds <- infecteds_generator(thetas, knots, population, max_day, form) # generate infecteds
    inc_period <- setNames(psnorm(0:symp_delay_lim, mean=inc[1], sd=inc[2], xi=exp(inc[3]))-psnorm(0:symp_delay_lim-1, mean=inc[1], sd=inc[2], xi=exp(inc[3])),0:symp_delay_lim)
  }
  
  if(form=="peak"){
    days <- -symp_delay_lim:symp_delay_lim
    infecteds <- infecteds_generator(thetas, knots, population, max_day, form, inf_lims=inf_lims) # generate infecteds
    peak_to_symp <- setNames(psnorm(days, mean=inc[1], sd=inc[2], xi=exp(inc[3]))-psnorm(days-1, mean=inc[1], sd=inc[2], xi=exp(inc[3])),days)
    #norm_peak_to_symp <- (peak_to_symp/sum(peak_to_symp),days)
  }

  if(test_form == "empirical") test_delay <- c(test_param, test6=max(0,1-sum(test_param)))
  
  pseeds = unlist(generateSeedVectors(ncores,1)) 
  
  #if(form =="incidence") p_ct_tau_mat <- ct_func2_cpp2(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], p=c(max_inf, ct_range[2]+1, ct_range[1], test_pop), nthreads=ncores, pseeds=pseeds) 
  # p_ct_tau_mat <- ct_func2_cpp_vlmax(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], p=c(max_inf-1, ct_range[2]+1, ct_range[1], test_pop), nthreads=ncores, pseeds=pseeds)
  if(form=="peak"){
    p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(max_inf-1, vl_range[1], vl_range[2], test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt)
    dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, rev(vl_vec))
    round(t(p_ct_tau_mat),2)
  } 
  if(form=="incidence"){
    #sourceCpp("will3p_inc.cpp")
    p_ct_tau_mat <- ct_func2_cpp_inc(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(max_inf-1, vl_range[2]+1, vl_range[1], test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt, stoch)
    dimnames(p_ct_tau_mat) <- list(0:(max_inf-1), c(ct_range[1]:(ct_range[2]),"negative"))
    #round(t(p_ct_tau_mat)[1:31,1:15],2)
  } 
  
  ## CALCULATE p_ct_tau
  p_array <- array(NA, dim=c(length(vl_vec)-1, max_day, max_test_day+1)) 
  
  # generate an array where the row pertains to ct value, the column to calendar day and the slice to days since onset
  for(j in 36:max_day) for(k in 0:(dim(p_array)[3]-1)){
    if(form=="incidence"){
      p_array[,j,k+1] <- t(c(infecteds[(j-k):(j-k-symp_delay_lim)]*inc_period[1:(symp_delay_lim+1)]*test_delay[k+1]))%*%p_ct_tau_mat[(k+1):(k+symp_delay_lim+1),1:ct_diff]
    } 
    if(form=="peak"){
      p_array[,j,k+1] <- t(c(infecteds[(j-symp_delay_lim-k):(j+symp_delay_lim-k)]*rev(peak_to_symp)*test_delay[k+1]))%*%p_ct_tau_mat[(inf_lims+1+symp_delay_lim+k):(inf_lims+1-symp_delay_lim+k),1:(ncol(p_ct_tau_mat)-1)]
    }
  }
  return(p_array)
}

infecteds_generator <- function(thetas, knots, population, max_day, form, inf_lims=NA){
  if(form=="incidence") infecteds <- population*exp(setNames(stats::splinefun(x=c(0,knots*max_day), y=thetas)(1:(max_day)),1:max_day))
  if(form=="peak") infecteds <- population*exp(setNames(stats::splinefun(x=c(0,knots*max_day,max_day+inf_lims), y=thetas)(1:(max_day+inf_lims)), 1:(max_day+inf_lims)))
  
  return(infecteds)
}

thetas_generator <- function(parameters){
  start <- which(names(parameters) == "theta1")
  end <- length(parameters)-1
  
  parameters[start:end] <- parameters[start:end]*parameters[length(parameters)]
  parameters[length(parameters)] <- 1
  
  return(parameters)
}
