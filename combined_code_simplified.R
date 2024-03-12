library(dqrng)
library(truncnorm)
library(Rcpp)
library(mvtnorm)
library(stats)
library(ggplot2)
library(dplyr)
library(MultinomialCI)
library(matrixStats)
library(ggpubr)
library(zoo)
library(fGarch)
library(ggpmisc)
library(psych)
library(MASS)
library(png)
library(ess)
library(ggh4x)
library(tmvtnorm)
library(TruncatedNormal)
library(hipercow)
library(accelerometry)

source("key_functions_simplified.R")

setwd("C:/Users/wg4618/Documents/ClusterPaper")

Sys.setenv(PATH = paste("C:/rtools40/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")

#vl_data_all = readRDS("data/data_array_vl_all_N.rds")
hipercow_init()
hipercow_configure(driver = "windows")
hipercow_configuration()
hipercow_provision()
hipercow_provision_check(show_unchanged=TRUE)
hipercow_provision_compare()

hipercow_environment_create(sources = c("cluster_code_simplified_react.R"))

input_generator <- function(n_iterations=5e6, strain="all", test_pop=1e7, ncores=8, proposal_type="MH", cov_matrix=NULL, population=5.5e7, form="peak", 
                            parameters=NA, cov_start=150000, max_test_day=6, gene="N", tag="", corr, flat){
  
  data_array <- readRDS(paste0("data/data_array_vl_",strain,"_N",".rds"))[,,1:(max_test_day+1)]
  if(strain != "all") data_array <- data_array[,-1,]
  if("vl_0.5" %in% rownames(data_array)) data_array <- data_array[-which(rownames(data_array)=="vl_0.5"),,]
  
  
  array_dim_names <- dimnames(data_array)
  
  n_vls <- dim(data_array)[1]
  vls <- as.numeric(substr(array_dim_names[[1]], 4, 4+nchar(array_dim_names[[1]])-3))
  
  n_days <- dim(data_array)[2]
  dates <- array_dim_names[[2]]
  days <- 1:n_days
  
  vg = corr
  
  n_test_days <- dim(data_array)[3]
  test_days <- as.numeric(substr(array_dim_names[[3]], 5, 5))
  
  if(!strain %in% c("all", "omicron")) knots <- seq(0,n_days,14)[-1]/n_days
  if(strain=="all")                    knots <- c(seq(0,490,14)[-1], seq(497,n_days,7))/n_days
  if(strain=="omicron")                knots <- c(seq(0,n_days-7,7))[-1]/n_days
  
  theta_names <- paste0("theta", seq(1:(length(knots)+1)))
  #theta_names_peak <- paste0("theta", seq(1:(length(knots)+2)))
  
  ## vl parameters
  vl_kinetics_type <- paste0(ifelse(flat==1, "flat-", ""),ifelse(corr %in% c(0,3), "corr0", ifelse(corr %in% c(1,4), "corr1", ifelse(corr %in% c(2,5), "corr2", ifelse(corr==6, "corr6", NA)))))
  
  priors_file <- paste0("Q:/ClusterPaper/longitudinal_fits/longitudinal_posteriors_",vl_kinetics_type,"_final.rds")
  print(priors_file)
  priors <- readRDS(priors_file)
  vl_prior_means <- priors$mean
  vl_prior_sds <- priors$sd
  
  priors_cov_file <- paste0("Q:/ClusterPaper/longitudinal_fits/longitudinal_cov_",vl_kinetics_type,"_final.rds")
  print(priors_cov_file)
  prior_cov <- readRDS(priors_cov_file)
  
  test_prior_means_raw <- setNames((apply(data_array, 3, sum)/sum(data_array)), paste0("test",seq(0,max_test_day)))
  test_prior_means <- setNames(vector(length=length(test_prior_means_raw)-1), names(test_prior_means_raw[-length(test_prior_means_raw)]))
  
  test_prior_means[1] <- test_prior_means_raw[1]
  for(i in 2:length(test_prior_means)) test_prior_means[i] = test_prior_means_raw[i]/prod(1-test_prior_means[1:(i-1)])
  
  test_prior_sds <-   setNames(rep(0.5, max_test_day), paste0("test",seq(0,max_test_day-1)))
  
  #spline_prior_means <- setNames(c(log(apply(data_array, 2, sum)[knots*n_days]/population),-9), theta_names)
  spline_prior_means <- setNames(rep(-9,length(theta_names)),theta_names)
  spline_prior_sds <- setNames(c(2, 2, rep(9,length(theta_names)-4), 1, 1),theta_names)
  
  spline_multiplier_mean <- c(multiplier=1)
  spline_multiplier_sd <- c(multiplier=10)
  
  ## incubation period and peak to symp delay parameter
  if(form == "incidence"){
    incu_prior_means <- c(inc1=4, inc2=2, inc3=0)
    incu_prior_sds <-   c(inc1=4, inc2=2, inc3=2)
  }
  
  if(form=="peak"){
    incu_prior_means <- c(inc1=0, inc2=2, inc3=0)
    incu_prior_sds <-   c(inc1=4, inc2=2, inc3=2)
  }
  
  if(form=="thresh"){
    incu_prior_means <- c(inc1=4, inc2=2, inc3=2)
    incu_prior_sds <-   c(inc1=4, inc2=2, inc3=2)
  }
  
  if(form=="thresh_peak"){
    incu_prior_means <- c(inc1=4, inc2=2, inc3=2)
    incu_prior_sds <-   c(inc1=4, inc2=2, inc3=2)
  }
  
  if(form=="thresh_tdist"){
    incu_prior_means <- c(inc1=1, inc2=0.5, inc3=4)
    incu_prior_sds   <- c(inc1=1, inc2=0.5, inc3=4) 
  }
  
  prior_mean <- c(vl_prior_means, incu_prior_means, test_prior_means, rdisp=1,  spline_prior_means, spline_multiplier_mean)
  prior_sds <-  c(vl_prior_sds,   incu_prior_sds,   test_prior_sds,   rdisp=10, spline_prior_sds,   spline_multiplier_sd)
  
  ignored_spline_par <- length(prior_mean) - ceiling(length(theta_names)/2)
  
  if(flat==0) ignored_pars <- c(7,8)
  else ignored_pars <- ignored_spline_par
  
  if(is.na(parameters)[1]){
    parameters <- c(prior_mean)
    parameters[-ignored_spline_par] <- rnorm(n=length(parameters)-1, mean=parameters[-ignored_spline_par], sd=abs(parameters[-ignored_spline_par])*0.05)
  }
  
  MCMC_sds <-abs(parameters)/200
  MCMC_sds['inc3'] <- 0.05
  
  max_inf <- 43
  if(form=="peak") symp_delay_lim <- 14
  if(form=="incidence") symp_delay_lim <- 28
  if(form %in% c("thresh","thresh_peak","thresh_tdist")) symp_delay_lim <- NA
  
  spline_pars = names(parameters)[which(grepl("theta",names(parameters)))]
  
  lower_bound_vl <- c(a_bar=-Inf, a_sigma=0.05, b_bar=-Inf, b_sigma=0.05, vl_max_bar=-Inf, vl_max_sigma=0, l_bar=-Inf, l_sigma=0)
  lower_bound_inc <- c(inc1=ifelse(form %in% c("peak", "thresh_tdist"), -Inf, 0), inc2=0, inc3=ifelse(form %in% c("incidence", "peak"), -Inf, 0))
  lower_bound_test <- c(test0=0, test1=0, test2=0, test3=0, test4=0, test5=0, rdisp=-Inf)
  lower_bound_splines <- setNames(rep(-Inf, length(spline_pars)), spline_pars)
  
  lower_bound <- c(lower_bound_vl, lower_bound_inc, lower_bound_test, lower_bound_splines)
  
  prior_cov_final <- matrix(0, nrow=length(parameters), ncol=length(parameters))
  colnames(prior_cov_final) <- names(parameters)
  rownames(prior_cov_final) <- names(parameters)
  
  prior_cov_final[1:8,1:8] <- prior_cov
  for(i in 9:length(parameters)) prior_cov_final[i,i] <- prior_sds[i]^2
  
  input <- list(n_iterations=1000000,
                parameters=parameters,
                vl_kinetics_type=vl_kinetics_type,
                flat=flat,
                corr=corr,
                knots=knots,
                data_array=data_array,
                test_pop=1e7,
                MCMC_sds=MCMC_sds,
                prior_mean=prior_mean,
                prior_sds=prior_sds,
                lower_bound=lower_bound,
                prior_cov_final=prior_cov_final,
                ncores=ncores,
                proposal_type=proposal_type,
                cov_matrix = cov_matrix,
                population = 5.6e7,
                form=form,
                tag=tag,
                ignored_spline_par = ignored_spline_par,
                strain=strain,
                max_inf=max_inf,
                symp_delay_lim=symp_delay_lim,
                max_day=n_days,
                cov_start=cov_start,
                gene=gene,
                n_vls = n_vls,
                vls=vls,
                n_days=n_days, 
                days=days,
                dates=dates,
                n_test_days=n_test_days,
                test_days=test_days,
                n_parameters=length(parameters),
                vg=vg,
                ignored_pars=ignored_pars)
  
  return(input)
}

#####
input_generator_react <- function(strain, n_iterations=5e5, parameters=NA, test_pop=1e7, proposal_type="MH", cov_matrix=NA, ncores=32, ignored_pars, cov_start=100000){
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
  prior_output <- output_finder(strain=strain, form="incidence", corr=4, flat=0)
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

x0 <- input_generator(corr=0, flat=0, form="thresh")
x1 <- input_generator_react(strain="all", n_iterations=5e5, parameters=NA, test_pop=1e7, proposal_type="MH", cov_matrix=NA, ncores=32, ignored_pars=c(7,8), cov_start=100000)

Mn <- task_create_expr(MCMC_function(x1), resources=resources)  
task_log_show(Mn)

output_continuing_function_react <- function(){
  files_list <- list.files(path="all/react/")
  times <- file.info(paste0("all/react/",files_list))$ctime
  output <- readRDS(paste0("all/react/", files_list[which(times==max(times))]))
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  input <- output$x
  input$proposal_type="Cov"
  input$cov_matrix = cov(output$MCMC_output[(n_it-20000):n_it,])/length(output$x$parameters)
  
  input$parameters <- output$MCMC_output[n_it,]
  
  return(input)
}

x_react <- output_continuing_function_react()
x_react$ncores <- 16
resources_react <- hipercow_resources(cores = 16)
M0 <- task_create_expr(MCMC_function(x_react), resources=resources_react)
task_log_show(M0)


#####

pars <- expand.grid(flat = c(0), corr = c(3, 4, 5), form = c("peak", "thresh", "thresh_peak", "thresh_tdist"), 
                    strain=c("wt", "alpha", "delta", "omicron")) %>% mutate(cov_start = c(50000), form=as.character(form), strain=as.character(strain)) 


burnin_calculator <- function(output){
  posteriors <- output$MCMC_posterior
  n_it <- which(is.na(posteriors))[1]-1
  window <- 10000
  rolling_posteriors <- c(rep(NA,window),accelerometry::movingaves(posteriors[1:n_it], window))
  
  burnin <- which(rolling_posteriors>rolling_posteriors[n_it])[1]
  
  return(burnin)
}

pars$chain_length <- 0
pars$chain_length[c(10,15)] <- c(50000,10000)
pars$cores <- 7
pars$cores[c(8:10,13:15)] <- 15

pars_filtered <- pars %>% filter(cores==15)

pars6$chain_length <- 100000
pars6$chain_length[8] <- 75000
pars6$ncores <- 8
pars6$ncores[c(8, 10, 11, 12)] <- 32

output_continuing_function <- function(pars, i, chain=2, ncores=8){
  form=pars$form[i]
  corr=pars$corr[i]
  flat=pars$flat[i]
  strain=pars$strain[i]
  
  output <- output_finder(strain, form, corr, flat, recentness=1)
  
  #input_new_priors <- input_generator(form=form, corr=corr, flat=flat, strain=strain)
  
  input <- output$x
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  burnin <- pars$chain_length[i]
  #else burnin <- n_it-pars$chain_length[i]
  #print(c(burnin, n_it), quote=F)
  
  input$parameters <- output$MCMC_output[n_it,]
  
  input$cov_matrix <- cov(output$MCMC_output[(n_it-burnin):n_it,])*4/length(input$parameters)

  input$chain <- chain
  
  #input$lower_bound <- input_new_priors$lower_bound
  #input$prior_mean <- input_new_priors$prior_mean
  #input$prior_cov_final <- input_new_priors$prior_cov_final
  #input$prior_sds <- input_new_priors$prior_sds
    
  if("vl_0.5" %in% rownames(input$data_array)) input$data_array <- input$data_array[rownames(input$data_array) != "vl_0.5",,]
  
  input$ncores <- pars$ncores[i]
  
  return(input)
}

resources <- hipercow_resources(cores = 8)

for(i in 1:nrow(pars6)){
  x <- output_continuing_function(pars6, i, chain=3, ncores=NA)
  resources <- hipercow_resources(cores = x$ncores)
  if(x$ncores==8)  M1 <- task_create_expr(MCMC_function(x), resources=resources)
  #if(x$ncores==32) M2 <- task_create_expr(MCMC_function(x), resources=resources)
}

#x1 <- output_continuing_function(pars_filtered, 1)
#x2 <- output_continuing_function(pars_filtered, 2)
#x3 <- output_continuing_function(pars_filtered, 3)

#M1 <- task_create_expr(MCMC_function(x1), resources=resources)
#M2 <- task_create_expr(MCMC_function(x2), resources=resources)
#M3 <- task_create_expr(MCMC_function(x3), resources=resources)

task_log_show(M1)
task_log_show(M2)
task_log_show(M3)

pars6 <- expand.grid(flat = c(0), corr = c(6), form = c("peak", "thresh", "thresh_peak", "thresh_tdist"), 
                    strain=c("wt", "alpha", "delta", "omicron")) %>% mutate(cov_start = c(50000), form=as.character(form), strain=as.character(strain)) 


pars_vg_6_starter <- function(pars6, i){
  most_similar_output <- output_finder(pars6$strain[i], pars6$form[i], corr=4, pars6$flat[i], 1)
  n_it <- which(is.na(most_similar_output$MCMC_posteriors))[1]-2
  
  input <- input_generator(strain=pars6$strain[i], form=pars6$form[i], corr=6, flat=0, proposal_type="Cov")
  input$cov_matrix <- most_similar_output$cov_matrix
  input$proposal_type <- "Cov"
  input$parameters[-c(1:11)] <- most_similar_output$MCMC_output[n_it, -c(1:11)]
  input$ncores <- 8
  input$tag <- "bounded_a_b_"
  
  x <- input
  return(x)
  
}

######

resources4 <- hipercow_resources(cores = 8)
x1 <- pars_vg_6_starter(pars6, 1)
M1 <- task_create_expr(MCMC_function(x1), resources=resources4)
task_log_show(M1)

x2 <- pars_vg_6_starter(pars=pars6, i=2)
M2 <- task_create_expr(MCMC_function(x2), resources=resources4)

x3 <- pars_vg_6_starter(pars6, 3)
M3 <- task_create_expr(MCMC_function(x3), resources=resources4)

x4 <- pars_vg_6_starter(pars6, 4)
M4 <- task_create_expr(MCMC_function(x4), resources=resources4)

x5 <- pars_vg_6_starter(pars6, 5)
M5 <- task_create_expr(MCMC_function(x5), resources=resources4)

x6 <- pars_vg_6_starter(pars6, 6)
M6 <- task_create_expr(MCMC_function(x6), resources=resources4)

x7 <- pars_vg_6_starter(pars6, 7)
M7 <- task_create_expr(MCMC_function(x7), resources=resources4)

x8 <- pars_vg_6_starter(pars6, 8)
M8 <- task_create_expr(MCMC_function(x8), resources=resources4)

x9 <- pars_vg_6_starter(pars6, 9)
M9 <- task_create_expr(MCMC_function(x9), resources=resources4)

x10 <- pars_vg_6_starter(pars6, 10)
M10 <- task_create_expr(MCMC_function(x10), resources=resources4)

x11 <- pars_vg_6_starter(pars6, 11)
M11 <- task_create_expr(MCMC_function(x11), resources=resources4)

x12 <- pars_vg_6_starter(pars6, 12)
M12 <- task_create_expr(MCMC_function(x12), resources=resources4)

x13 <- pars_vg_6_starter(pars6, 13)
M13 <- task_create_expr(MCMC_function(x13), resources=resources4)

x14 <- pars_vg_6_starter(pars6, 14)
M14 <- task_create_expr(MCMC_function(x14), resources=resources4)

x15 <- pars_vg_6_starter(pars6, 15)
M15 <- task_create_expr(MCMC_function(x15), resources=resources4)

x16 <- pars_vg_6_starter(pars6, 16)
M16 <- task_create_expr(MCMC_function(x16), resources=resources4)

#####


chain_conv_check <- function(pars, i, ncores){
  input <- output_continuing_function(pars, i, ncores) 
  
  input$parameters[1:18] <- input_generator(strain=pars$strain[i],
                                      corr=pars$corr[i], 
                                      flat=pars$flat[i], 
                                      form=pars$form[i])$parameters[1:18]
  
  input$parameters[19:length(input$parameters)] <- 
    input$parameters[19:length(input$parameters)]*
    rnorm(n=length(input$parameters[19:length(input$parameters)]), mean=1, sd=0.04)
  
  input$tag <- "checking"
  
  return(input)
}

for(i in 1:nrow(pars_select)){
  x <- chain_conv_check(pars_select, i, ncores)
  M1 <- task_create_expr(MCMC_function(x), resources=resources)
}

task_log_show(M1)

pars_strain <- expand.grid(flat = c(0),
                    corr = c(3, 4, 5),
                    form = c("thresh_tdist"), strain=c("wt", "alpha", "delta", "omicron")) %>%
  mutate(cov_start = c(0), form=as.character(form), strain=as.character(strain)) 

thresh_tdist_start <- function(pars_strain, i){
  input <- input_generator(strain=pars_strain$strain[i], proposal_type="Cov", form=pars_strain$form[i], corr=pars_strain$corr[i], flat=pars_strain$flat[i])
  
  similar_output <- output_finder(strain=pars_strain$strain[i], form="thresh", corr=pars_strain$corr[i], flat=pars_strain$flat[i])
  n_it <- which(is.na(similar_output$MCMC_posteriors))[1]-2
  
  input$parameters <- similar_output$MCMC_output[n_it,]
  input$cov_matrix <- similar_output$cov_matrix
  input$ncores <- 16
  
  return(input)
}

resources16 <- hipercow_resources(cores = 16)

for(i in 1:nrow(pars_strain)){
  x5 <- thresh_tdist_start(pars_strain, i)
  M5 <- task_create_expr(MCMC_function(x5), resources=resources16)
}

task_log_show(M5)


strain_starting_function <- function(pars_strain, i){
  input <- input_generator(strain=pars_strain$strain[i], ncores=8, proposal_type="MH", form=pars_strain$form[i], 
                           parameters=NA, corr=pars_strain$corr[i], flat=pars_strain$flat[i], cov_start=pars_strain$cov_start[i])
  
  input_knot_dates <- as.Date(as.Date(min(input$dates))+dim(input$data_array)[2]*c(0,input$knots))
  input_params <- input$parameters
  input_theta_params <- which(grepl("theta", names(input_params)))
  
  all <- output_finder("all", form=pars_strain$form[i], corr=pars_strain$corr[i], flat=pars_strain$flat[i])
  n_it <- which(is.na(all$MCMC_posteriors))[1]-2
  
  all_knot_dates <- as.Date(min(all$x$dates))+dim(all$x$data_array)[2]*c(0,all$x$knots)
  all_params <- all$MCMC_output[n_it,]
  all_thetas <- all_params[grepl("theta", names(all_params))]
  all_df <- data.frame(dates=all_knot_dates, thetas=all_thetas)
  
  combined_dates <- as.Date(c(input_knot_dates, all_knot_dates))
  
  combined_df <- data.frame(dates = as.Date(min(combined_dates):max(combined_dates))) %>%
    left_join(., all_df, by="dates") %>% 
    mutate(value=approx(dates, thetas, xout = dates, method = "linear", rule = 2)$y)
  
  input$parameters[1:(min(input_theta_params)-1)] <- all_params[1:(min(input_theta_params)-1)] 
  input$parameters[input_theta_params] <- combined_df %>% filter(dates %in% input_knot_dates) %>% dplyr::select(value) %>% pull()
  
  input$cov_start <- pars_strain$cov_start[i]
  
  return(input)
}

for(i in 1:nrow(pars)){
  x <- strain_starting_function(pars_strain, i)
  M1 <- task_create_expr(MCMC_function(x), resources=resources)
}
task_log_show(M5)


#####
resources <- hipercow_resources(cores = 16)
pars <- pars %>% filter(form %in% c("thresh", "thresh_peak"), strain != "all")

x1 <- output_continuing_function(pars, i=1)
M1 <- task_create_expr(MCMC_function(x1), resources=resources)
task_log_show(M1)

x2 <- output_continuing_function(pars, i=2)
M2 <- task_create_expr(MCMC_function(x2), resources=resources)
task_log_show(M2)

x3 <- output_continuing_function(pars, i=3)
M3 <- task_create_expr(MCMC_function(x3), resources=resources)
task_log_show(M3)

x4 <- output_continuing_function(pars, i=4)
M4 <- task_create_expr(MCMC_function(x4), resources=resources)
task_log_show(M4)

x5 <- output_continuing_function(pars, i=5)
M5 <- task_create_expr(MCMC_function(x5), resources=resources)
task_log_show(M5)

x6 <- output_continuing_function(pars, i=6)
M6 <- task_create_expr(MCMC_function(x6), resources=resources)
task_log_show(M6)

x7 <- output_continuing_function(pars, i=7)
M7 <- task_create_expr(MCMC_function(x7), resources=resources)
task_log_show(M7)

x8 <- output_continuing_function(pars, i=8)
M8 <- task_create_expr(MCMC_function(x8), resources=resources)
task_log_show(M8)

x9 <- output_continuing_function(pars, i=9)
M9 <- task_create_expr(MCMC_function(x9), resources=resources)
task_log_show(M9)

x10 <- output_continuing_function(pars, i=10)
M10 <- task_create_expr(MCMC_function(x10), resources=resources)
task_log_show(M10)

x11 <- output_continuing_function(pars, i=11)
M11 <- task_create_expr(MCMC_function(x11), resources=resources)
task_log_show(M11)

x12 <- output_continuing_function(pars, i=12)
M12 <- task_create_expr(MCMC_function(x12), resources=resources)
task_log_show(M12)

x13 <- output_continuing_function(pars, i=13)
M13 <- task_create_expr(MCMC_function(x13), resources=resources)
task_log_show(M13)

x14 <- output_continuing_function(pars, i=14)
M14 <- task_create_expr(MCMC_function(x14), resources=resources)
task_log_show(M14)

x15 <- output_continuing_function(pars, i=15)
M15 <- task_create_expr(MCMC_function(x15), resources=resources)
task_log_show(M15)

x16 <- output_continuing_function(pars, i=16)
M16 <- task_create_expr(MCMC_function(x16), resources=resources)
task_log_show(M16)

x17 <- output_continuing_function(pars, i=17)
M17 <- task_create_expr(MCMC_function(x17), resources=resources)
task_log_show(M17)

x18 <- output_continuing_function(pars, i=18)
M18 <- task_create_expr(MCMC_function(x18), resources=resources)
task_log_show(M18)

x19 <- output_continuing_function(pars, i=19)
M19 <- task_create_expr(MCMC_function(x19), resources=resources)
task_log_show(M19)

x20 <- output_continuing_function(pars, i=20)
M20 <- task_create_expr(MCMC_function(x20), resources=resources)
task_log_show(M20)

x21 <- output_continuing_function(pars, i=21)
M21 <- task_create_expr(MCMC_function(x21), resources=resources)
task_log_show(M21)

x22 <- output_continuing_function(pars, i=22)
M22 <- task_create_expr(MCMC_function(x22), resources=resources)
task_log_show(M22)

x23 <- output_continuing_function(pars, i=23)
M23 <- task_create_expr(MCMC_function(x23), resources=resources)
task_log_show(M23)

x24 <- output_continuing_function(pars, i=24)
M24 <- task_create_expr(MCMC_function(x24), resources=resources)
task_log_show(M24)

x2 <- output_continuing_function(pars, 14)
x2$parameters[c(2,4)] <- 0.5
x2$vg <- 0
x2$corr <- 0
M2 <- task_create_expr(MCMC_function(x2), resources=resources)
task_log_show(M2)

x26 <- output_continuing_function(pars, 38)
x26$parameters[c(2,4)] <- 0.5
x26$vg <- 0
x26$corr <- 0
M26 <- task_create_expr(MCMC_function(x26), resources=resources)
task_log_show(M26)

x30 <- output_continuing_function(pars, 6)
x30$ignored_pars <- x30$ignored_spline_par
x30$parameters[c("l_bar", "l_sigma")] = c(-5,1)
x30$cov_matrix[c(7,8),c(7,8)] <- output_finder("all", "incidence", corr=1, flat=1)$x$cov_matrix[c(7,8),c(7,8)]
M30 <- task_create_expr(MCMC_function(x30), resources=resources)
task_log_show(M30)

x42 <- output_continuing_function(pars, 18)
x42$ignored_pars <- x42$ignored_spline_par
x42$parameters[c("l_bar", "l_sigma")] = c(-5,1)
x42$cov_matrix[c(7,8),c(7,8)] <- output_finder("all", "incidence", corr=4, flat=1)$x$cov_matrix[c(7,8),c(7,8)]
M42 <- task_create_expr(MCMC_function(x42), resources=resources)
task_log_show(M42)
#####


output_combining_function <- function(strain, form, corr, flat){
  outputs <- file.info(paste0(strain,"/",output_list(strain, form, corr, flat, recentness))) %>% arrange(ctime)
  output_names <- rownames(outputs)
  
  output <- readRDS(output_names[1])
  next_row <- which(is.na(output$MCMC_output[,1]))[1]
  output$chain <- vector(length=length(output$MCMC_posteriors))
  output$chain[1:next_row] <- 1
  
  for(i in 2:length(output_names)){
    print(i)
    if(length(output_names)==1) break
    output_next <- readRDS(output_names[i])
    n_it <- which(is.na(output_next$MCMC_output[,1]))[1]-1
    print(n_it)
    new_indexes <- next_row:(next_row+n_it-1)
    output$MCMC_output[new_indexes,] <- output_next$MCMC_output[1:n_it,]
    output$MCMC_likelihoods[new_indexes] <- output_next$MCMC_likelihoods[1:n_it,]
    output$MCMC_posteriors[new_indexes] <- output_next$MCMC_posteriors[1:n_it,]
    output$Acceptances[new_indexes] <- output_next$Acceptances[1:n_it]
    output$MCMC_likelihoods[new_indexes] <- output_next$MCMC_likelihoods[1:n_it,]
    output$time_vec[new_indexes] <- output_next$time_vec[1:n_it]
    output$chain[new_indexes] <- max(output$chain, na.rm=T)+1
    #output
    next_row <- next_row + n_it
  }
  
  output$cov_matrix <- output$cov_matrix
  output$x$lower_bound["vl_max_sigma"] <- 0
  
  n_tot_raw <- which(is.na(output$MCMC_posteriors))[1]
  n_tot <- format(n_tot_raw - n_tot_raw %% 100, scientific=F)
  
  saveRDS(list(MCMC_output=output$MCMC_output, MCMC_posteriors=output$MCMC_posteriors, MCMC_likelihoods=output$MCMC_likelihoods, Acceptances=output$Acceptances, cov_matrix=output$cov_matrix, time_vec=output$time_vec, x=output$x, chain=output$chain), 
          file=paste0(strain,"/",output$x$strain,"_",output$x$gene,"_",output$x$form,"_",output$x$proposal_type,"_ignored_",gsub(" ","",toString(output$x$ignored_pars)),"_",output$x$vl_kinetics_type,"_",output$x$tag,"vg=",output$x$vg,"_",n_tot,"_",Sys.Date(),"_mult.rds"))
  
  cat("\n")
}

pars_new <- pars6[10:12,]# %>% filter(strain=="omicron")
for(i in 1:nrow(pars_new)) output_combining_function(strain=pars_new$strain[i], form=pars_new$form[i], corr=pars_new$corr[i], flat=pars_new$flat[i])

out <- output_finder(strain="all", form="thresh_peak", corr=2, flat=1)

comparison_loop_1D <- function(output, n_samples, chain_length=NA){
  
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  
  if(is.na(chain_length)) start = 1
  else start <- n_it - chain_length
  
  d_sums <- list()
  storage <- list()
  plot_dfs <- list()
  plot_list <- list()
  variables <- c("viral_load", "day", "test_day")
  
  for (d in 1:3){
    d_sums[[d]] <- apply(output$x$data_array[-1,,], d, sum, na.rm=T)
    storage[[d]] <- matrix(nrow=n_samples, ncol=dim(output$x$data_array[-1,,])[d])
  } 
  
  iterations <- round(seq(start, n_it, length.out=n_samples))
  print(iterations)
  
  for(i in 1:length(iterations)){
    parameters <- output$MCMC_output[iterations[i],]
    #parameters[c("test2", "test3", "test4")] = parameters[c("test2", "test3", "test4")]*0.8 
    p_array_raw <- p_array_func(parameters, knots=output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, population=output$x$population, 
                            data_array=output$x$data_array, test_pop=1e7, ncores=12, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim, 
                            stoch=0.5, name=T)[[1]]
    p_array_raw[is.na(p_array_raw)] <- 0
    
    p_array <-  array(rnbinom(n=length(p_array_raw), size=1/parameters['rdisp'], mu=p_array_raw), dim=dim(p_array_raw))
    dimnames(p_array) <- dimnames(p_array_raw)
    
    for(d in 1:3) storage[[d]][i,] <- apply(p_array, d, sum, na.rm=T)
  }
  
  for(d in 1:3){
    plot_dfs[[d]] <- data.frame(x=as.numeric(dimnames(p_array)[[d]]), 
                                p=colQuantiles(storage[[d]], probs=c(0.025, 0.5, 0.975)), 
                                d=d_sums[[d]]) %>% magrittr::set_colnames(c(variables[d], "lower", "median", "upper", "count")) 
    
    if(d==2) plot_dfs[[d]] <- plot_dfs[[2]] %>% filter(count>100)
    
    plot_list[[d]] <- ggplot(plot_dfs[[d]], aes(x=.data[[variables[d]]])) +
      geom_point(aes(y=count)) +
      geom_line(aes(y=median)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
      theme_bw() + 
      scale_y_log10()
    
    #print(plot_list[[d]])
    print(d)
  } 
  
  return(plot_list)
}

comparison_loop_2D <- function(output, n_samples, chain_length=NA){
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  start <- n_it - chain_length
  iterations <- round(seq(start, n_it, length.out=n_samples))
  
  storage <- array(dim=c(dim(output$x$data_array[-1,,])[1], dim(output$x$data_array[-1,,])[3], n_samples))
  
  for(i in 1:length(iterations)){
    parameters <- output$MCMC_output[iterations[i],]
    #parameters[c("test2", "test3", "test4")] = parameters[c("test2", "test3", "test4")]*0.8 
    
    p_array_raw <- p_array_func(parameters, output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, population=output$x$population, 
                            data_array=output$x$data_array, test_pop=1e7, ncores=12, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim, 
                            stoch=0.5, name=T)[[1]]
    p_array_raw[is.na(p_array_raw)] <- 0
    
    p_array <-  array(rnbinom(n=length(p_array_raw), size=1/parameters['rdisp'], mu=p_array_raw), dim=dim(p_array_raw))
    dimnames(p_array) <- dimnames(p_array_raw)
    
    storage[,,i] <- round(apply(p_array, c(1,3), sum, na.rm=T))
    dimnames(storage) = list(dimnames(p_array)[[1]],dimnames(p_array)[[3]],1:n_samples)
  }
  
  cut_data = output$x$data_array[-1,,]
  dimnames(cut_data) = dimnames(p_array)
  
  d_sum <- apply(cut_data, c(1,3), sum, na.rm=T) %>% as.table() %>% as.data.frame() %>% 
    magrittr::set_colnames(c("viral_load", "test_day", "number")) 
  
  plotting_df <- apply(storage, c(1,2), quantile, probs=c(0.025, 0.5, 0.975)) %>% as.table() %>% as.data.frame() %>% 
    tidyr::spread(Var1, Freq) %>%
    magrittr::set_colnames(c("viral_load", "test_day", "lower", "median", "upper")) %>%
    left_join(d_sum, by=c("viral_load","test_day")) %>%
    mutate(viral_load=as.numeric(as.character(viral_load)), test_day=as.numeric(as.character(test_day)))
  
  output_plot <- ggplot(plotting_df, aes(x=viral_load)) +
    geom_point(aes(y=number)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3) +
    facet_wrap(~test_day) +
    theme_bw() + 
    scale_y_log10()
  
  return(output_plot)
  
}

output_plotter <- function(output, chain_length, n_samples=10){
  #output <- output_finder(strain, form, corr, flat)
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  
  MCMC_output <- output$MCMC_output[1:n_it,]
  MCMC_posteriors <- output$MCMC_posteriors[1:n_it]
  data_array <- output$x$data_array
  knots <- output$x$knots
  population <- output$x$population
  max_day <- output$x$max_day
  form <- output$x$form
  symp_delay_lim <- output$x$symp_delay_lim
  
  parameters <- MCMC_output[n_it,]
  thetas <- thetas_generator(parameters)[which(names(parameters)=="theta1"):(length(parameters)-1)]
  
  infecteds <- infecteds_generator(thetas, knots, population, max_day, form)
  #plot(infecteds, type="l")
  
  comp1 <- comparison_loop_1D(output, n_samples, chain_length=chain_length)
  comp2 <- comparison_loop_2D(output, n_samples, chain_length=chain_length)
  
  plot <- ggpubr::ggarrange(comp1[[1]], comp1[[2]], comp1[[3]], comp2, nrow=2, ncol=2)
  print(plot)
  ggsave(paste0(output$x$strain,"/images/fitting/",output$x$form,"_corr=", output$x$corr,"_flat=", output$x$flat,"_",output$x$tag,
                      "_vg=", output$x$vg,"_mult3.png"), plot, width=7, height=7)
}

pars_omicron <- pars6

for(i in 1:nrow(pars_omicron)){
  output <- output_finder(strain=pars_omicron$strain[i], pars_omicron$form[i], pars_omicron$corr[i], pars_omicron$flat[i])
  if(length(output)==1) next
  output_plotter(output=output, chain_length=20000, n_samples=10)
} 

param_output_post_func <- function(output){
  output <- output_finder(strain="wt", form="peak", corr=3, flat=0)
  n_it <-  which(is.na(output$MCMC_posteriors))[1]-2
  parameters_vg2 <- output$MCMC_output[n_it,]
  
  plot(output$MCMC_posteriors[seq(1000,235000,100)], type="l")
  
  n_it0 <- 100000
  n_it1 <- 200000
  #n_it2 <- which(is.na(output2$MCMC_posteriors))[1]-2
  
  output$MCMC_output[n_it,] <- output$MCMC_output[n_it0,]
  
  output_plotter(output, n_it, n_samples=1)
  
  post <- posterior(parameters=output$MCMC_output[n_it,], prior_mean=output$x$prior_mean, prior_cov_final=output$x$prior_cov_final, 
                    lower_bound=output$x$lower_bound, knots=output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, 
                    population=output$x$population, data_array=weekly_aggregator(output$x$data_array), 
                    test_pop=output$x$test_pop, ncores=output$x$ncores, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim, 
                    ignored_pars=output$x$ignored_pars)
  
  print(post)
  
  vl <- vls_calculator2(strain=output$x$strain, output$x$form, output$x$corr, output$x$flat, n_it1, 1)
  
  ggplot(vl, aes(x=t, y=median)) + geom_line() +
    theme_bw()
}

MCMC_plotter <- function(output){
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  thinned_iter <- round(seq(1,n_it,length.out=10000))
  
  combined <- cbind(iteration=thinned_iter, output$MCMC_output[thinned_iter,], posteriors=output$MCMC_posteriors[thinned_iter]) %>%
    as.data.frame() %>% mutate(posteriors2 = ifelse(iteration < n_it/2, NA, posteriors)) 
  
  variable_set <- c("viral", "theta", "other")
  
  viral_params <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma", "inc1", "inc2", "inc3")
  other_params <- c("test0", "test1", "test2", "test3", "test4", "test5", "rdisp", "posteriors", "posteriors2")
  theta_params <- colnames(combined)[grepl("theta", colnames(combined))]
  
  plot_list <- list()
  
  for(i in 1:3){
    filtered_combined <- combined[,c("iteration", get(paste0(variable_set[i],"_params")))] %>%
      tidyr::gather(param, value, 2:ncol(.)) %>% mutate(param=factor(param, levels=colnames(combined))) 
    
    plot_list[[i]] <- ggplot(filtered_combined, aes(x=iteration, y=value)) +
      geom_line() +
      facet_wrap(~param, scales="free_y") +
      theme_bw()
  }
  
  return(plot_list)
  
}

a <- MCMC_plotter(output_finder("delta", "incidence", corr=4, flat=0))
a[[1]]
a[[2]]
a[[3]]

MCMC_comparator <- function(pars, strain_interest=NA, chain_length, thinning){
  out <- data.frame(number=numeric(), row=numeric(), burnin=numeric(), strain=character(), form=character(), corr=numeric(), flat=numeric(), posterior=numeric(), chain=numeric())
  next_row <- 1
  if(is.na(strain_interest)==F) pars <- pars %>% filter(strain==strain_interest)
  
  list <- list()
  
  for(i in 1:nrow(pars)){
    print(i)
    output <- output_finder(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
    if(length(output)==1) next
    n_it <- which(is.na(output$MCMC_posteriors))[1]-2
    if(is.null(output$chain)) output$chain=rep(1,n_it)
    if(is.na(chain_length)) start <- 1
    else if(chain_length < 1) start <- (1-chain_length)*n_it
    else start <- n_it - chain_length
    out_length <- round((n_it-start)*thinning)
    burnin <- burnin_calculator(output)
    out[next_row:(next_row+out_length-1),] = data.frame(number=i,
                                                        row=round(seq(start,n_it,length.out=out_length)),
                                                        burnin=burnin,
                                                        strain=pars$strain[i],
                                                        form=pars$form[i],
                                                        corr=pars$corr[i],
                                                        flat=pars$flat[i],
                                                        posterior=output$MCMC_posteriors[round(seq(start,n_it,length.out=out_length))],
                                                        chain=output$chain[round(seq(start,n_it,length.out=out_length))])
    
    next_row=next_row+out_length
  }
  
  out$flat = ifelse(out$flat==0, "tri", "trap")
  out$flat = factor(out$flat, levels=c("tri", "trap"))
  out$chain <- factor(out$chain)
  
  p <- ggplot(out, aes(x=row, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    facet_nested(strain ~ corr + form, scales="free_x", labeller = label_wrap_gen(multi_line=FALSE), independent="x") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-4)) +
    labs(x="Number of iterations (tens thousands)")
  
  p2 <- ggplot(out, aes(x=row, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    facet_nested(strain ~ corr + form, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-4)) +
    labs(x="Number of iterations (tens thousands)", y="posterior (tens)")
  
  p3 <- ggplot(out %>% filter(row > burnin) %>% mutate(row_norm=row-burnin), aes(x=row_norm, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    facet_nested(strain ~ form + corr, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-3)) +
    labs(x="Number of iterations (thousands)", y="posterior (tens)")
  
  ggsave("all/images/MCMC1.png", p, width=10, height=5)
  ggsave("all/images/MCMC2.png", p2, width=10, height=5)
  
  return(list(p, p2, p3))
}

MCMC_facet <- MCMC_comparator(pars=pars6, strain_interest=NA, chain_length=0.5, thinning=0.25)
MCMC_facet[[1]]
MCMC_facet[[2]]
MCMC_facet[[3]]

ess_calculator <- function(strain, form, corr, flat, start){
  output <- output_finder(strain=strain, form=form, corr=corr, flat=flat)
  if(length(output)==1) return(c(na=NA))
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  if(n_it<start) return(c(na=NA))
  return(sns::ess(output$MCMC_output[start:n_it,]))
}

ess_comparator <- function(pars, start){
  for(i in 1:nrow(pars)){
    print(i)
    ess_hold <- round(ess_calculator(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], start=start))
    if(i==1) ess <- data.frame(param=names(ess_hold), form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], ess=ess_hold) 
    else ess <- rbind(ess, data.frame(param=names(ess_hold), form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], ess=ess_hold)) 
  }
  rownames(ess) <- NULL
  return(ess %>% tidyr::spread(param, ess))
}

ess <- ess_comparator(pars, start=5000)

AIC_calculator <- function(output, start){
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  n_pars <- length(output$x$parameters)-length(unique(c(output$x$ignored_pars, output$x$ignored_spline_par)))
  
  likelihood_vec <- vector(length = n_it - start)
  for(i in 1:length(likelihood_vec)) likelihood_vec[i] = output$MCMC_posteriors[i] - prior(output$MCMC_output[i,], output$x$prior_mean, output$x$prior_sds, output$x$form, output$x$ignored_pars)
  
  maximum_likelihood <- max(likelihood_vec)
  AIC <- 2*n_it-2*maximum_likelihood
  
  return(AIC)
}

AIC_comparator <- function(pars, start){
  out <- data.frame(form=character(), corr=logical(), flat=numeric(), AIC=numeric())
  
  for(i in 1:nrow(pars)){
    output <- output_finder(strain="all", form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
    AIC <- AIC_calculator(output, start)
    out[i,1:3] <- pars[i,1:3]
    out[i,4] <- AIC
  }
  cat("\n")
  
  out$AIC <- out$AIC-min(out$AIC)
  out$corr = ifelse(out$corr==TRUE, "corr", "no_corr")
  out$flat = ifelse(out$flat==0, "tri", "trap")
  
  p <- ggplot(out, aes(x=form, y=AIC, fill=corr)) +
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~flat, scales="free") +
    theme_bw()
  
  print(p)
  
  return(out)
}

DIC_calculator <- function(output, start){
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  likelihoods_vec <- output$MCMC_likelihoods[start:n_it]
  deviance_vec <- -2*likelihoods_vec
  
  expectation_of_deviance <- mean(deviance_vec)
  deviance_of_expectation <- mean(max(unique(deviance_vec)[100]))
  
  DIC = 2*expectation_of_deviance-deviance_of_expectation
  
  return(DIC)
}

DIC_comparator <- function(pars, start){
  out <- data.frame(form=character(), corr=logical(), flat=numeric(), DIC=numeric())
  
  for(i in 1:nrow(pars)){
    output <- output_finder(strain="all", form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
    DIC <- DIC_calculator(output, start)
    out[i,1:3] <- pars[i,1:3]
    out[i,4] <- DIC
  }
  cat("\n")
  
  out$DIC <- out$DIC-min(out$DIC)
  out$corr = ifelse(out$corr==TRUE, "corr", "no_corr")
  out$flat = ifelse(out$flat==0, "tri", "trap")
  
  p <- ggplot(out, aes(x=form, y=DIC, fill=corr)) +
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~flat, scales="free") +
    theme_bw()
  
  print(p)
  
  return(out)
}

DICs <- DIC_comparator(pars=pars, start=90000)
DICs

metrics_calculator <- function(output){
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  acceptances <- sum(output$Acceptances)/n_it
  time_it <- sum(output$time_vec)/n_it
  return(c(n_it, acceptances, time_it))
}

metrics_comparator <- function(pars){
  for(i in 1:nrow(pars)){
    print(i)
    output <- output_finder(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
    metrics <- metrics_calculator(output)
    
    start <- burnin_calculator(output)
    if((metrics[[1]] - start) < 10000) start <- metrics[[1]] - 10000
    if(length(output)==1){
      pars[i,c("cov_start", "n_it", "Acceptances", "time_it", "min_ess", "min_ess_par", "DIC")] <- NA
      next
    } 
    ess <- ess_calculator(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], start=start)
    
    pars$n_it[i] <- metrics[1]
    pars$warm_up[i] <- start
    pars$chain_length[i] <- metrics[[1]] - start
    pars$Acceptances[i] <- round(metrics[2],2)
    pars$time_it[i] <- round(metrics[3],2)
    pars$min_ess[i] <- round(min(ess[ess!=0]),1) 
    pars$min_ess_par[i] <- ifelse(is.na(ess[1]), NA, names(which(ess==min(ess[ess!=0]))))
    pars$DIC[i] <- DIC_calculator(output, start)
    pars$posterior[i] <- output$MCMC_posteriors[metrics[1]]
    
    #pars$a_lower[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"a_bar"])-median(output$MCMC_output[start:metrics[1],"a_sigma"])),2)
    #pars$a_median[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"a_bar"])),2)
    #pars$a_upper[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"a_bar"])+median(output$MCMC_output[start:metrics[1],"a_sigma"])),2)
    
    #pars$a[i] <- paste0(pars$a_median[i], " (", pars$a_lower[i], " - ", pars$a_upper[i], ")")
    #
    #pars$b_lower[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"b_bar"])-median(output$MCMC_output[start:metrics[1],"b_sigma"])),2)
    #pars$b_median[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"b_bar"])),2)
    #pars$b_upper[i] <- round(exp(median(output$MCMC_output[start:metrics[1],"b_bar"])+median(output$MCMC_output[start:metrics[1],"b_sigma"])),2)
    #
    #pars$b[i] <- paste0(pars$b_median[i], " (", pars$b_lower[i], " - ", pars$b_upper[i], ")")
    
    #pars$vlmax_lower[i] <- round(median(output$MCMC_output[start:metrics[1],"vl_max_bar"])-median(output$MCMC_output[start:metrics[1],"vl_max_sigma"]),2)
    #pars$vlmax_median[i] <- round(median(output$MCMC_output[start:metrics[1],"vl_max_bar"]),2)
    #pars$vlmax_upper[i] <- round(median(output$MCMC_output[start:metrics[1],"vl_max_bar"])+median(output$MCMC_output[start:metrics[1],"vl_max_sigma"]),2)
    
    #pars$inc1[i] <- round(median(output$MCMC_output[start:metrics[1],"inc1"]),2)
    #pars$inc2[i] <- round(median(output$MCMC_output[start:metrics[1],"inc2"]),2)
    #pars$inc3[i] <- round(median(output$MCMC_output[start:metrics[1],"inc3"]),2)
    
    #pars$vlmax[i] <- paste0(pars$vlmax_median[i], " (", pars$vlmax_lower[i], " - ", pars$vlmax_upper[i], ")")
    
  }
  
  pars_new <- pars %>% group_by(strain, corr) %>% filter(form!="incidence") %>% mutate(DIC=round(DIC-min(DIC, na.rm=T))) %>% mutate(posterior_norm=round(posterior-max(posterior, na.rm=T),0))
  
  #pars_new %>% dplyr::select(strain, form, corr, min_ess, DIC, posterior_norm, a, b, vlmax, inc1, inc2, inc3) %>% arrange(strain, DIC) %>% print(n=50)
  
  cat("\n")
  
  p1 <- ggplot(pars_new[,c("strain", "form", "corr", "DIC")] %>% mutate(corr=factor(corr)), aes(x=corr, fill=form, y=DIC)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(~strain, scales="free_y") +
    theme_bw()
  
  p1
  
  p2 <- ggplot(pars_new, aes(x=form, y=posterior_norm)) +
    geom_point() +
    facet_grid(strain~corr, scales="free_y") +
    theme_bw()
  
  p2
  
  ggsave("all/images/model_DIC.png", p1, width=10, height=5)
  ggsave("all/images/model_posteriors.png", p2, width=10, height=5)
  
  
  return(list(DIC=p1, posteriors=p2, metrics=pars_new))
}

metrics <- metrics_comparator(pars=pars6)
metrics[[2]]
metrics[[3]] %>% group_by(strain) %>% arrange(strain, posterior_norm, desc=T)
metrics[[3]] %>% arrange(form, corr, flat)

vl_func <- function(a1, b1, l1, t1, vlmax1, vi, corr, flat){
  a_par <- exp(a1)
  tmax1 <- (vlmax1-vi)/a_par
  
  if(corr %in% c(4,5)) b_par <- exp(b1)/exp(a1)
  else b_par <- exp(b1)
  
  if(corr == 5) vlmax_par = vlmax1*exp(a1)
  else vlmax_par = vlmax1
  
  if(corr == 6){
    a_par <- 0.5+exp(a1)
    b_par <- 0.25+exp(b1)
  } 
  
  
  l_par <- ifelse(flat==0, 0, exp(l1))
  
  vl = ifelse(t1< -l_par/2, vlmax_par+a_par*(t1+l_par/2), ifelse(t1<l_par/2, vlmax_par, vlmax_par-b_par*(t1-l_par/2)))
    
  return(pmax(vl,vi))
}

vls_calculator <- function(strain, form, corr, flat, chain_length, n_samples=33, n_indiv=33){
  times <- seq(-14,28,0.1)
  vl_storage <- array(dim=c(length(times), n_indiv, n_samples))
  t_peak <- array(dim=c(n_indiv, n_samples))
  symp_rel_peak <- array(dim=c(n_indiv, n_samples))
  symp_rel_det <- array(dim=c(n_indiv, n_samples))
  tot_infection_time <- array(dim=c(n_indiv, n_samples))
  
  output_raw <- output_finder(strain=strain, form=form, corr=corr, flat=flat)
  if(length(output_raw) == 1) return("skip")
  output <- output_raw$MCMC_output
  
  n_it <- which(is.na(output))[1]-2
  start <- n_it - chain_length
  iterations <- seq(start,n_it,length.out=n_samples)
  
  for(i in 1:n_samples){
    vl_max_samples <- rnorm(n=n_indiv, output[iterations[i],"vl_max_bar"], output[iterations[i],"vl_max_sigma"])
    
    a_samples_vl <- rnorm(n=n_indiv, mean=output[iterations[i],"a_bar"], sd=output[iterations[i],"a_sigma"])
    b_samples_vl <- rnorm(n=n_indiv, mean=output[iterations[i],"b_bar"], sd=output[iterations[i],"b_sigma"])
    
    a_samples <- 0.5+exp(a_samples_vl)
    b_samples <- 0.25+exp(b_samples_vl)
    
    for(m in 1:n_indiv){
      vl_storage[,m,i] <- vl_func(a1=a_samples_vl[m], b1=b_samples_vl[m], l1=l_samples[m], t1=times, vlmax1=vl_max_samples[m], vi=1, corr=output_raw$x$corr, flat=output_raw$x$flat)
    }
    
    if(output_raw$x$form=="peak" & output_raw$x$vg==6){
      days = -output_raw$x$symp_delay_lim:output_raw$x$symp_delay_lim
      
      inc_dist = setNames(fGarch::psnorm(days,   mean=output[iterations[i],"inc1"], sd=output[iterations[i],"inc2"], xi=exp(output[iterations[i],"inc3"]))-
                          fGarch::psnorm(days-1, mean=output[iterations[i],"inc1"], sd=output[iterations[i],"inc2"], xi=exp(output[iterations[i],"inc3"])),days)
      
      time_of_peak = (vl_max_samples-vi)/a_samples
      inc_samples = sample(days, size=n_samples, prob=inc_dist, replace=T)
      t_symptoms = time_of_peak+inc_samples
      accepted <- rep(T,length(inc_samples))
    } 
    if(output_raw$x$form=="thresh" & output_raw$x$vg==6){
      inc_samples <- rnorm(n=n_indiv, mean=output[iterations[i],"inc1"], sd=output[iterations[i],"inc2"])
      
      vl_thresh <- inc_samples
      accepted <- vl_thresh < vl_max_samples
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+output[iterations[i],"inc3"]
    } 
    if(output_raw$x$form=="thresh_peak" & output_raw$x$vg==6){
      inc_samples <- rtruncnorm(n=n_indiv, a=0, mean=output[iterations[i],"inc1"], sd=output[iterations[i],"inc2"])
      
      vl_thresh <- vl_max_samples-inc_samples
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+output[iterations[i],"inc3"]
      accepted <- rep(T,length(inc_samples))
    } 
    if(output_raw$x$form=="thresh_tdist" & output_raw$x$vg==6){
      inc_samples <- exp(rnorm(n=n_indiv, mean=output[iterations[i],"inc1"], sd=output[iterations[i],"inc2"]))
      
      vl_thresh <- output[iterations[i],"inc3"]
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+inc_samples
      accepted <- rep(T,length(inc_samples))
    } 
    
    t_peak[,i] <- (vl_max_samples-vi)/a_samples
    
    symp_rel_det[,i] <- t_symptoms
    symp_rel_peak[,i] <- t_symptoms-t_peak[,i]
    tot_infection_time[,i] <- (vl_max_samples-vi)*(a_samples+b_samples)/(a_samples*b_samples)
    
    symp_rel_det[!accepted,i] <- NA
    symp_rel_peak[!accepted,i] <- NA
    tot_infection_time[!accepted,i] <- NA
  }
  
  dimnames(vl_storage)[[1]] <- round(times,2)
  
  median_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, median))
  lower_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, quantile, probs=0.025))
  upper_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, quantile, probs=0.975))
  
  plotting_df <- data.frame(t=rep(times,3),
                            form=form, 
                            strain=strain,
                            corr=corr, 
                            flat=flat,
                            type=rep(c("median","lower","upper"), c(length(times),length(times), length(times))),
                            vls=rbind(rowQuantiles(median_each_slice, probs=c(0.025, 0.5, 0.975)),
                                      rowQuantiles(lower_each_slice, probs=c(0.025, 0.5, 0.975)),
                                      rowQuantiles(upper_each_slice, probs=c(0.025, 0.5, 0.975)))) %>%
    magrittr::set_colnames(c("t", "form", "strain", "corr", "flat", "type", "lower", "median", "upper"))
  
  p <- ggplot(plotting_df, aes(x=t, color=type, fill=type)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color="white") +
    theme_bw()
  
  plotting_df2 <- data.frame(t=rep(times,3),
                            form=form, 
                            strain=strain,
                            corr=corr, 
                            flat=flat,
                            vls=unname(t(apply(vl_storage, MARGIN=c(1), quantile, probs=c(0.025, 0.5, 0.975))))) %>%
    magrittr::set_colnames(c("t", "form", "strain", "corr", "flat", "lower", "median", "upper"))
  
  p2 <- ggplot(plotting_df2, aes(x=t)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color="white") +
    theme_bw()
  
  #print(p2)
  
  plotting_df3 <- reshape2::melt(vl_storage, varnames=c("t", "person", "sample")) %>%
    mutate(ID=paste0(person,sample)) %>%
    filter(!lag(value,1) < 1.000001 |
           !lead(value,1) < 1.00000001 |
            t %in% c(0,-14,28)) %>%
    group_by(ID) %>%
    filter(value < 1.000001 | value==max(value) | t %in% c(-14,28)) %>%
    mutate(form=form, 
           strain=strain,
           corr=corr, 
           flat=flat)
    
  p3 <- ggplot(plotting_df3, aes(x=t, y=value, group=ID)) +
    geom_line(color="blue", alpha=0.2) +
    theme_bw() +
    theme(legend.position="none") +
    geom_line(data=plotting_df2, aes(x=t, y=median, group=1), size=1)
    
  key_times <- data.frame(form=form, 
                             strain=strain, 
                             corr=corr, 
                             flat=flat,
                             det_to_symp = c(symp_rel_det),
                             peak_to_symp = c(symp_rel_peak), 
                             det_to_peak = c(t_peak),
                             total_infection_time = c(tot_infection_time)) %>% 
    tidyr::gather(param, value, 5:8) 
  
  p4 <- ggplot(key_times, aes(x=param, y=value)) +
    geom_boxplot() +
    #facet_wrap() +
    theme_bw()
  
  #ggsave(paste0(strain,"/vl_trajectories/",strain,"_",form,"_corr=",corr,"_flat=",flat), p)
  
  return(list(plotting_df, plotting_df2, plotting_df3, key_times))
}

hist_plotter <- function(strain, form, corr, flat, chain_length, n_samples){
  #times <- seq(-14,28,0.1)
  #vl_storage <- array(dim=c(length(times), n_samples, n_samples))
  
  output <- output_finder(strain=strain, form=form, corr=corr, flat=flat)
  vi <- as.numeric(substr(rownames(output$x$data_array)[1],4,4))
  
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  start <- n_it-chain_length
  
  median_params <- setNames(colMedians(output$MCMC_output[(n_it-start):n_it,]), colnames(output$MCMC_output))
  
  vl_max_samples = rnorm(n=n_samples, median_params["vl_max_bar"], median_params["vl_max_sigma"])
  a_samples = rnorm(n=n_samples, mean=median_params["a_bar"], sd=median_params["a_sigma"])
  b_samples = rnorm(n=n_samples, mean=median_params["b_bar"], sd=median_params["b_sigma"])
  
  if(output$x$form %in% c("peak", "incidence")){
    if(output$x$form == "incidence") days = 0:output$x$symp_delay_lim
    if(output$x$form == "peak") days = -output$x$symp_delay_lim:output$x$symp_delay_lim
    
    inc_dist = setNames(fGarch::psnorm(days,   mean=median_params["inc1"], sd=median_params["inc2"], xi=exp(median_params["inc3"]))-
                        fGarch::psnorm(days-1, mean=median_params["inc1"], sd=median_params["inc2"], xi=exp(median_params["inc3"])),days)
    
    time_of_peak = (vl_max_samples-vi)/exp(a_samples)
    time_of_symptoms = sample(days, size=n_samples, prob=inc_dist, replace=T)
    time_of_symptom_rel_peak = time_of_symptoms - time_of_peak
    vl_at_symptoms = vl_func(a1=a_samples, b1=b_samples, l1=0, t1=time_of_symptom_rel_peak, vlmax1=vl_max_samples, vi=vi, corr=output$x$corr, flat=output$x$flat)
  }
  
  if(output$x$form == "peak"){
    vt_samples = rnorm(n=n_samples, median_params["inc1"], sd=median_params["inc2"])
    lag = median_params["inc3"]
    time_of_peak = (vl_max_samples-vi)/exp(a_samples)
    time_of_threshold = (vt_samples-vi)/exp(a_samples)
    time_of_symptoms = time_of_threshold+lag
    time_of_symptom_rel_peak = time_of_symptoms - time_of_peak
    vl_at_symptoms = vl_func(a1=a_samples, b1=b_samples, l1=0, t1=time_of_symptom_rel_peak, vlmax1=vl_max_samples, vi=vi, corr=output$x$corr, flat=output$x$flat)
  }
  
  if(output$x$form == "thresh"){
    vt_samples = rnorm(n=n_samples, median_params["inc1"], sd=median_params["inc2"])
    lag = median_params["inc3"]
    time_of_peak = (vl_max_samples-vi)/exp(a_samples)
    time_of_threshold = (vt_samples-vi)/exp(a_samples)
    time_of_symptoms = time_of_threshold+lag
    time_of_symptom_rel_peak = time_of_symptoms - time_of_peak
    vl_at_symptoms = vl_func(a1=a_samples, b1=b_samples, l1=0, t1=time_of_symptom_rel_peak, vlmax1=vl_max_samples, vi=vi, corr=output$x$corr, flat=output$x$flat)
  } 
  
  if(output$x$form == "thresh_peak"){
    vt_peak_samples = rtruncnorm(n=n_samples, a=0, median_params["inc1"], sd=median_params["inc2"])
    vt_samples = vl_max_samples - vt_peak_samples
    lag = median_params["inc3"]
    time_of_peak = (vl_max_samples-vi)/exp(a_samples)
    time_of_threshold = time_of_peak-vt_peak_samples/exp(a_samples)
    time_of_symptoms = time_of_threshold+lag
    time_of_symptom_rel_peak = time_of_symptoms - time_of_peak
    vl_at_symptoms = vl_func(a1=a_samples, b1=b_samples, l1=0, t1=time_of_symptom_rel_peak, vlmax1=vl_max_samples, vi=vi, corr=output$x$corr, flat=output$x$flat)
  } 
    
  if(output$x$form == "thresh_tdist"){
    vt_samples = median_params["inc3"]
    lag = exp(rnorm(n=n_samples, median_params["inc1"], sd=median_params["inc2"]))
    time_of_peak = (vl_max_samples-vi)/exp(a_samples)
    time_of_threshold = (vt_samples-vi)/exp(a_samples)
    time_of_symptoms = time_of_threshold+lag
    time_of_symptom_rel_peak = time_of_symptoms - time_of_peak
    vl_at_symptoms = vl_func(a1=a_samples, b1=b_samples, l1=0, t1=time_of_symptom_rel_peak, vlmax1=vl_max_samples, vi=vi, corr=output$x$corr, flat=output$x$flat)
  } 
   
  combined <- data.frame(a=exp(a_samples),
                         b=exp(b_samples),
                         vlmax=vl_max_samples,
                         time_of_peak=time_of_peak,
                         time_of_symptoms=time_of_symptoms,
                         vl_at_symptoms=vl_at_symptoms) %>%
    mutate(vthresh=ifelse(form %in% c("thresh", "thresh_peak", "thresh_tdist"), unname(vt_samples), NA)) 
  
  if(form %in% c("thresh", "thresh_peak", "thresh_tdist")) combined <- combined %>% 
    mutate(vthresh = unname(vt_samples)) %>%
    filter(vthresh < vlmax)
  
  p <- ggplot(combined %>% tidyr::gather("param", "value", 1:ncol(.)), aes(x=value)) +
    geom_histogram(fill="blue") +
    facet_wrap(~param, scales="free") +
    theme_bw()
  
  p
  
  ggsave(paste0(output$x$strain,"/images/histplots/median_histplot_",output$x$strain,"_",output$x$form,"_corr_",output$x$corr,"_vg_",output$x$vg,".png"), p)
  
  #for(m in 1:n_samples){
    #  vl_storage[,m,i] <- vl_func(a1=a_samples[m], b1=b_samples[m], l1=l_samples[m], t1=times, vlmax1=vl_max_samples[m], vi=0.5, corr=output$x$corr, flat=output$x$flat)
    #}
  #}
  
  #median_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, median))
  #lower_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, quantile, probs=0.025))
  #upper_each_slice <- apply(vl_storage, MARGIN = 3, FUN = function(slice) apply(slice, MARGIN = 1, quantile, probs=0.975))
  
  #plotting_df <- data.frame(t=rep(times,3),
  #                          strain=strain,
  #                          form=form, 
  #                          corr=corr, 
  #                          flat=flat,
  #                          type=rep(c("median","lower","upper"), c(length(times),length(times), length(times))),
  #                          vls=rbind(rowQuantiles(median_each_slice, probs=c(0.025, 0.5, 0.975)),
  #                                    rowQuantiles(lower_each_slice, probs=c(0.025, 0.5, 0.975)),
  #                                    rowQuantiles(upper_each_slice, probs=c(0.025, 0.5, 0.975)))) %>%
  #  magrittr::set_colnames(c("t", "strain", "form", "corr", "flat", "type", "lower", "median", "upper"))
  
  #p <- ggplot(plotting_df, aes(x=t, color=type, fill=type)) +
  #  geom_line(aes(y=median)) +
  #  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color="white") +
  #  theme_bw()
  
  #print(p)
  
  return(p)
}

for(i in 1:nrow(pars)) hist_plotter(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], chain_length=3000, n_samples=1000)

vls_comparator <- function(pars, chain_length){
  out <- data.frame(t=numeric(), form=character(), strain=character(), corr=logical(), flat=character(), lower=numeric(), median=numeric(), upper=numeric())
  out2 <- data.frame(t=numeric(), person=numeric(), sample=numeric(), value=numeric(), ID=character(), form=character(), strain=character(), corr=numeric(), flat=numeric())
  out_symp <- data.frame(form=character(), strain=character(), corr=numeric(), flat=numeric(), param=character(), value=numeric())
  
  next_line <- 1
  next_line2 <- 1
  
  for(i in 1:nrow(pars)){
    vls <- vls_calculator(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], chain_length=chain_length, n_samples=10)
    if(length(vls) != 1){
      out[((i-1)*nrow(vls[[2]])+1):(i*nrow(vls[[2]])),] <- vls[[2]]
      out2[next_line:(next_line+nrow(vls[[3]])-1),] <- vls[[3]]
      out_symp[next_line2:(next_line2+nrow(vls[[4]])-1),] <- vls[[4]]
      
      next_line <- next_line+nrow(vls[[3]])-1
      next_line2 <- next_line2+nrow(vls[[4]])-1
    } 
  }
  
  p <- ggplot(out %>% filter(corr!=5), aes(x=t)) +
    geom_line(aes(y=median)) +
    facet_nested(form ~ strain + corr, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color="white") +
    theme_bw()
  
  p2 <- ggplot(out2 %>% filter(is.na(strain)==FALSE), aes(x=t, y=value, group=ID)) +
        geom_line(color="blue", alpha=0.1) +
        theme_bw() +
        theme(legend.position="none") +
        facet_nested(form ~ strain + corr, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
        geom_line(data=out, aes(x=t, y=median, group=1), size=1)
  
  p3 <- ggplot(out_symp, aes(x=form, y=value)) +
    geom_boxplot(outlier.shape=NA) +
    facet_nested(param ~ strain + corr, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    theme_bw()
  
  
  ggsave("images/viral_load_comparison_aggregage.png", p)
  ggsave("images/viral_load_comparison_samples.png", p2)
  ggsave("images/symptom_comparisons.png", p3)
  
  print(p)
  
  return(out)
}

vl_comparison <- vls_comparator(pars=pars6, chain_length=20000)

vl_comparison[[3]]

output <- output_finder("all", "peak", flat=1, corr=2)
start=40000

p_vl_day <- function(output, chain_length, n_samples){
  out <- data.frame(day=character(), vl=character())
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  start <- n_it - chain_length
  iterations <- round(seq(start,n_it,length.out=n_samples))
  
  for(i in 1:length(iterations)){
    parameters <- output$MCMC_output[i,]
    p_ct_tau_mat <- p_array_func(parameters, knots=output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, population=output$x$population, data_array=output$x$data_array, test_pop=output$x$test_pop, ncores=output$x$ncores, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim, stoch=0.5, name=F)[[2]]
    next_line <- 1
    
    for(j in 1:nrow(p_ct_tau_mat)){
      vls <- rep(colnames(p_ct_tau_mat), p_ct_tau_mat[j,]*1000)
      out[next_line:(next_line+length(vls)-1),] = data.frame(day=rownames(p_ct_tau_mat)[j], vl=vls)
      next_line <- next_line+length(vls)
    }
  }
  out_plot <- out %>% mutate(day=factor(day, levels=rownames(p_ct_tau_mat)), 
                             vl=ifelse(vl=="negative", "0.5", vl)) %>%
    mutate(vl=as.numeric(vl))
  
  p <- ggplot(out_plot, aes(x=day, y=vl)) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw()

  print(p)
  
  return(out_plot)
}

p_vl_day(output, 5000, n_samples=10)

p_vl_day_comparator <- function(pars, chain_length, n_samples){
  p_vl_day_df <- data.frame(day=factor(), vl=numeric(), strain=character(), form=character(), corr=numeric(), flat=numeric())
  
  for(i in 1:nrow(pars)){
    print(i)
    output <- output_finder(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
    if(length(output)==1) next
    p_vl_day_i <- p_vl_day(output, chain_length, n_samples) 
    if(i==1) p_vl_day_df <- p_vl_day_i %>% mutate(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i]) 
    else p_vl_day_df <- rbind(p_vl_day_df, p_vl_day_i %>% mutate(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i]))
  }
  
  p <- ggplot(p_vl_day_df %>% mutate(day=factor(day, levels=-20:40)), aes(x=day, y=vl)) +
    geom_boxplot(outlier.shape=NA) +
    facet_nested(form ~ corr + strain, scales="free_x", labeller = label_wrap_gen(multi_line=FALSE), independent="x") +
    theme_bw() +
    scale_x_discrete(breaks = as.character(seq(-20,40,10)))
  
  print(p)
  
  p2 <- ggplot(p_vl_day_df %>% mutate(day=factor(day, levels=-20:40), corr=factor(corr, levels=c(0,1,2))), aes(x=day, y=vl)) +
    geom_boxplot(outlier.shape=NA) +
    facet_nested(form ~ corr + trunc + flat, scales="free_x", labeller = label_wrap_gen(multi_line=FALSE), independent="x") +
    theme_bw() +
    scale_x_discrete(breaks = as.character(seq(-20,40,10)))
  
  print(p)
  
  ggsave(paste0("images/boxplot.png"), plot=p, width=12, height=7)
  ggsave(paste0("/images/boxplot2.png"), plot=p2, width=12, height=7)

}

p_vl_day_comparator(pars=pars6[c(1,5,9,13),], chain_length=2000, n_samples=10)

inc_period_plotter <- function(output, n_samples){
  form <- output$x$form
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  start <- burnin_calculator(output)
  
  if(form %in% c("incidence", "peak")){
    incubation <- p_array_func(parameters=output$MCMC_output[180000,], knots=output$x$knots, vk=1, vg=output$x$vg,
                               max_day=output$x$max_day, population=output$x$population, data_array=output$x$data_array, 
                               test_pop=output$x$test_pop, ncores=4, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim, stoch=0.5)[[3]]
    
    incubation_df <- data.frame(days=as.numeric(names(incubation)), prob=incubation) %>% mutate(cum_prob=cumsum(prob), prob=prob/max(prob))
    
    p <- ggplot(incubation_df, aes(x=days)) +
      geom_line(aes(y=prob)) +
      geom_line(aes(y=cum_prob)) +
      theme_bw()
  }
  
  if(form %in% c("thresh", "thresh_peak")){
    samples <- data.frame(samples=rnorm(n=10000, mean=output$MCMC_output[n_it,c("inc1")], sd=output$MCMC_output[n_it,c("inc2")]))
    
    p <- ggplot(samples, aes(x=samples)) +
      geom_histogram(aes(color="white")) +
      theme_bw() +
      geom_text(x=max(samples)-4, y=nrow(samples)/10, label=paste0("delay from threshold = ", round(output$MCMC_output[n_it, c("inc3")],1))) +
      theme(legend.position = "none")
    
  }
  print(p)
}

inc_period_plotter(output=output_finder(strain="wt", form="incidence", corr=5, flat=0), n_samples=10)
inc_period_plotter(output=output_finder(strain="wt", form="peak", corr=5, flat=0), n_samples=10)
inc_period_plotter(output=output_finder(strain="wt", form="thresh", corr=5, flat=0), n_samples=10)
inc_period_plotter(output=output_finder(strain="wt", form="thresh_peak", corr=5, flat=0), n_samples=10)

#####
output <- readRDS("all/2024-01-17_all_105000_N_incidence_MH_ignored_NA_non-differentiable___.rds")

N_array <- p_array_func(parameters=out2$MCMC_output[610000,], knots=out2$x$knots, vk=1, max_day=out2$x$max_day, population=out2$x$population, 
                        data_array, test_pop=1e7, ncores=24, form=out2$x$form, symp_delay_lim=out2$x$symp_delay_lim, stoch=0.5)
#####

#drat:::add("mrc-ide")

#options(
#  didehpc.username = "wg4618",   # change to your username 
#  didehpc.home = "Q:/")    # can be any letter 

#didehpc:::didehpc_config(credentials=list(username="wg4618"), home = "Q:/", r_version="4.2.1")
#
#context::context_log_start()
#
#root <- "context"

#ctx <- context::context_save(root, packages=c("dplyr", "ggplot2", "zoo", "matrixStats", "Rcpp", "dqrng", "truncnorm", "mvtnorm", "TruncatedNormal", "BH", "sitmo", "fGarch"),
#                             sources=c("cluster_code_simplified.R"))

#obj1 <- didehpc::queue_didehpc(ctx, config = didehpc:::didehpc_config(cores=8, credentials=list(username="wg4618"), home = "Q:/", r_version="4.1.3", cluster="wpia-hn", template="AllNodes"))
#obj2 <- didehpc::queue_didehpc(ctx, config = didehpc:::didehpc_config(cores=20, credentials=list(username="wg4618"), home = "Q:/", r_version="4.1.3", cluster="fi--didemrchnb", template="20Core"))

output_starting_func <- function(pars, i, trunc=T){
  input <- input_generator(parameters=NA, cov_start=pars$cov_start[i], max_test_day=6, gene="N", tag="", corr=pars$corr[i], flat=pars$flat[i], form=pars$form[i], trunc=trunc)
  most_similar_output <- output_finder(strain="all", form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
  
  if(length(most_similar_output)==1){
    most_similar_output <- output_finder(strain="all", form=pars$form[i], corr=0, flat=0)
    n_it <- which(is.na(most_similar_output$MCMC_posteriors))[1]-2
    input$parameters[12:length(input$parameters)] <- most_similar_output$MCMC_output[n_it,(12:length(input$parameters))]
    #input$cov_start <- 200000
  } 
  else{
    n_it <- which(is.na(most_similar_output$MCMC_posteriors))[1]-2
    input$parameters <- most_similar_output$MCMC_output[n_it,]
  }
  
  input$cov_start <- max(input$cov_start-n_it,50001)
  input$proposal_type <- "Cov"
  input$cov_matrix <- most_similar_output$cov_matrix
  
  input$ncores <- 16
  
  return(input)
}

switch_to_mult_function <- function(pars, i){
  form=pars$form[i]
  corr=pars$corr[i]
  flat=pars$flat[i]
  
  output <- output_finder(strain="all", form, corr, flat)
  input <- output$x
  n_it <- which(is.na(output$MCMC_posteriors))[1]-1
  
  input$parameters <- output$MCMC_output[n_it,]
  theta_pars <- which(grepl("theta",names(input$parameters))==TRUE)
  input$parameters[theta_pars] <- input$parameters[theta_pars]*input$parameters['multiplier']
  input$parameters['multiplier'] <- 1
  input$cov_matrix["multiplier", "multiplier"] <- input$cov_matrix["multiplier", "multiplier"]*10
  
  return(input)
}
