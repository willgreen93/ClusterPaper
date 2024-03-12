library(dplyr)
library(ggplot2)
library(matrixStats)
source("key_functions_simplified.R")

parameters <- c(a_bar=log(2), a_sigma=1, b_bar=1, b_sigma=1, vl_max_bar=8, vl_max_sigma=1, 
                l_bar=0, l_sigma=0, theta1=-9, theta2=-9, theta3=-9, theta4=-9, theta5=-9, 
                theta6=-9, theta7=-8, theta8=-10, theta9=-9, theta10=-9, multiplier=1)

input_generator_react <- function(strain, n_iterations=5e5, parameters=NA, test_pop=1e7, proposal_type="MH", cov_matrix=NA, ncores=16, ignored_pars, cov_start=100000){
  data_array_raw <- readRDS(paste0("data/REACT_data_matrix_",strain,".rds"))
  
  max_day <- dim(data_array_raw)[2]
  
  ### sort knots
  round_start_days <- intersect(which(colSums(data_array_raw)!=0), which(colSums(data_array_raw)==0)+1)-14
  round_end_days <- c(intersect(which(colSums(data_array_raw)==0), which(colSums(data_array_raw)!=0)+1), max_day)
  round_mid_points1 <- floor((2*round_start_days+round_end_days)/3)
  round_mid_points2 <- floor((round_start_days+2*round_end_days)/3)
  nround <- length(round_start_days)
  
  knots <- data.frame(round=seq(1:nround), start=round_start_days, mid1=round_mid_points1, mid2=round_mid_points2, end=round_end_days) %>%
    mutate(close1=as.numeric(lead(start,1))-as.numeric(end),
           close2=as.numeric(lead(start,1))-as.numeric(mid2)) %>%
    mutate(end=ifelse(close1<=7, NA, end),
           mid2=ifelse(close2<=7 & is.na(close2)==F, NA, mid2)) %>%
    dplyr::select(-round, -close1, -close2) %>%
    tidyr::gather(metric, day, 1:ncol(.)) %>%
    dplyr::select(day) %>% na.omit() %>% pull() %>% sort()

  p <- ggplot(data.frame(days=1:max_day, samples=colSums(data_array_raw)), aes(x=days, y=samples)) +
    geom_point() +
    geom_vline(data=data.frame(knots=c(0,knots,max_day)), aes(xintercept=knots), color="red", linetype="dotted") +
    theme_bw() +
    theme(panel.grid=element_blank())
    
  print(p)
  
  ### sort priors
  prior_output <- output_finder(strain=strain, form="incidence", corr=4, flat=0, 2)
  n_it <- which(is.na(prior_output$MCMC_posteriors))[1]-1
  
  vl_prior_means <- colMeans(prior_output$MCMC_output[(n_it-20000):n_it,1:8])
  vl_prior_sds <- setNames(colSds(prior_output$MCMC_output[(n_it-20000):n_it,1:8]), names(vl_prior_means))
  
  theta_names <- paste0("theta", seq(1:(length(knots)+2)))
  spline_prior_means <- setNames(rep(-7.5, length(theta_names)), theta_names)
  spline_prior_sds <- abs(spline_prior_means)
  
  prior_mean <- c(vl_prior_means, spline_prior_means, multiplier=1)
  prior_sds <- c(vl_prior_sds, spline_prior_sds, multiplier=10)
  n_parameters <- length(prior_mean)
  
  prior_cov_final <- matrix(0, nrow=n_parameters, ncol=n_parameters)
  diag(prior_cov_final) <- prior_sds^2
  colnames(prior_cov_final) <- names(prior_mean)
  rownames(prior_cov_final) <- names(prior_mean)
  
  MCMC_sds <- abs(prior_mean/100) 
  
  if(is.na(parameters)==T){
    parameters <- prior_mean
    print("setting starting parameters as prior means", quote=F)
  } 
  
  parameters[ignored_pars] <- 0
  MCMC_sds[ignored_pars] <- 0
  
  print(paste0("the proposal type is ", proposal_type), quote=F)
  
  lower_bound_vl <- prior_output$x$lower_bound[1:8]
  lower_bound_spline <- setNames(rep(-Inf, length(which(grepl("theta",names(prior_mean))))),theta_names)
  
  lower_bound <- c(lower_bound_vl, lower_bound_spline)
  
  ignored_spline_par <- which(names(parameters)==paste0("theta",round(length(theta_names)/2)))
  
  x <- list(n_iterations=n_iterations, 
            parameters=parameters, 
            knots=c(knots/max_day,1), 
            prior_mean=prior_mean, 
            prior_cov_final=prior_cov_final,
            population=5.5e7,
            data_array_raw=data_array_raw, 
            test_pop=test_pop,
            MCMC_sds=MCMC_sds,
            proposal_type=proposal_type,
            cov_matrix=cov_matrix,
            ncores=ncores,
            strain=strain, 
            ignored_pars=ignored_pars,
            vk=1,
            form=prior_output$x$form,
            vg=prior_output$x$vg,
            lower_bound=lower_bound,
            max_day=max_day,
            n_parameters=n_parameters,
            ignored_spline_par=ignored_spline_par,
            cov_start=cov_start)
  
  return(x)
}

x <- input_generator_react(strain="all", n_iterations=5e5, parameters=NA, test_pop=1e7, proposal_type="MH", 
                     cov_matrix=NA, ncores=16, ignored_pars=c(7,8), cov_start=100000)

M4 <- task_create_expr(MCMC_function(x), resources=resources)
task_log_show(M4)


knots <- seq(0.1,1,length.out=9)
vk <- 1
vg <- 3
population <- 5.5e7
data_array_raw <- readRDS("data/REACT_data_array_all.rds")[,-c(1:3)]
max_day <- dim(data_array)[2]
test_pop <- 1e7
ncores <- 8
form <- "incidence"
symp_delay_lim <- 28

proposal_type <- "MH"
strain <- "all"
ignored_spline_par <- 23
n_parameters <- length(parameters)
MCMC_sds <- abs(prior_mean)/100

ignored_pars <- c(7,8)
cov_start=100000
#vl_kiinetics_type <- 1


#likelihood_REACT <- function(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim){
  p_array <- p_array_func(parameters, knots, vk, vg, max_day, population, data_array, test_pop, ncores, form, symp_delay_lim, stoch=0.5, name=T)
  
  infecteds <- round(infecteds_generator(parameters, knots, population, max_day, form))
  
  N_matrix <- matrix(NA, nrow=nrow(data_array), ncol=ncol(data_array))
  rownames(N_matrix) <- c(rev(rownames(data_array)[-1]),"negative") 
  colnames(N_matrix) <- colnames(data_array)
  
  for(j in 31:ncol(N_matrix)){
    N_matrix[,j] <- t(p_array) %*% infecteds[j:(j-29)]
  }
  
  N_matrix[nrow(N_matrix),] <- population-colSums(N_matrix[1:(nrow(N_matrix)-1),])
  p_matrix <- N_matrix/population
  
  likelihood <- sum(data_array[,(31:max_day)] * log(p_matrix[,31:max_day]))
  
  return(likelihood)
}

