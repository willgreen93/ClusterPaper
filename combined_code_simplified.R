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
#library(hipercow)
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

burnin_calculator <- function(output){
  posteriors <- output$MCMC_posterior
  n_it <- which(is.na(posteriors))[1]-1
  window <- 10000
  rolling_posteriors <- c(rep(NA, window), rollmean(posteriors[1:n_it], k = window, align = "right", fill = NA))
  
  burnin <- which(rolling_posteriors>rolling_posteriors[n_it])[1]
  
  return(burnin)
}

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

pars6 <- expand.grid(flat = c(0), corr = c(6), form = c("peak", "thresh", "thresh_peak", "thresh_tdist"), 
                    strain=c("wt", "alpha", "delta", "omicron")) %>% mutate(cov_start = c(50000), form=as.character(form), strain=as.character(strain)) 


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

pars_new <- pars6
for(i in 1:nrow(pars_new)) output_combining_function(strain=pars_new$strain[i], form=pars_new$form[i], corr=pars_new$corr[i], flat=pars_new$flat[i])

out <- output_finder(strain="alpha", form="thresh_peak", corr=6, flat=0)

comparison_loop_1D <- function(output, n_samples, start){
  
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  
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
    print(i)
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
                                p=matrixStats::colQuantiles(storage[[d]], probs=c(0.025, 0.5, 0.975)), 
                                d=d_sums[[d]]) %>% magrittr::set_colnames(c(variables[d], "lower", "median", "upper", "count")) 
    
    if(d==2) plot_dfs[[d]] <- plot_dfs[[2]] %>% dplyr::filter(count>100)
    
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

comparison_loop_2D <- function(output, n_samples, start){
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  
  iterations <- round(seq(start, n_it, length.out=n_samples))
  
  storage <- array(dim=c(dim(output$x$data_array[-1,,])[1], dim(output$x$data_array[-1,,])[3], n_samples))
  
  for(i in 1:length(iterations)){
    print(i)
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

output_plotter <- function(output, chain_length, n_samples=100){
  #output <- output_finder(strain, form, corr, flat)
  n_it <- which(is.na(output$MCMC_output[,1])==TRUE)[1]-2
  if(chain_length<1) start = (1-chain_length)*n_it
  else start = n_it - chain_length
  
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
  
  comp1 <- comparison_loop_1D(output, n_samples, start=start)
  comp2 <- comparison_loop_2D(output, n_samples, start=start)
  
  plot <- ggpubr::ggarrange(comp1[[1]], comp1[[2]], comp1[[3]], comp2, nrow=2, ncol=2)
  print(plot)
  ggsave(paste0(output$x$strain,"/images/fitting/",output$x$strain,"_",output$x$form,"_corr=", output$x$corr,"_flat=", output$x$flat,"_",output$x$tag,
                      "_vg=", output$x$vg,"_mult3.png"), plot, width=7, height=7)
}

pars_omicron <- pars6 %>% dplyr::filter(form=="thresh_peak")

for(i in 1:nrow(pars_omicron)){
  output <- output_finder(strain=pars_omicron$strain[i], pars_omicron$form[i], pars_omicron$corr[i], pars_omicron$flat[i])
  if(length(output)==1) next
  output_plotter(output=output, chain_length=pars_omicron$chain_length[i], n_samples=100)
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

MCMC_comparator <- function(pars, strain_interest=NA, chain_length=NA, thinning){
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
    if(is.na(chain_length)) start=(1-pars$chain_length[i])*n_it
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
  
  p <- ggplot(out %>% ungroup(), aes(x=row, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    ggh4x::facet_nested(strain ~ form, scales="free_x", labeller = label_wrap_gen(multi_line=FALSE), independent="x") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-4)) +
    labs(x="Number of iterations (tens thousands)")
  
  p2 <- ggplot(out, aes(x=row, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    ggh4x::facet_nested(strain ~ corr + form, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-4)) +
    labs(x="Number of iterations (tens thousands)", y="posterior (tens)")
  
  p3 <- ggplot(out %>% dplyr::filter(row > burnin) %>% mutate(row_norm=row-burnin), aes(x=row_norm, y=posterior/10)) +
    geom_line(aes(color=chain)) +
    ggh4x::facet_nested(strain ~ form, scales="free", labeller = label_wrap_gen(multi_line=FALSE), independent="all") +
    theme_bw() +
    scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-3)) +
    labs(x="Number of iterations (thousands)", y="posterior (tens)")
  
  ggsave("images/MCMC1.png", p, width=10, height=5)
  ggsave("images/MCMC2.png", p2, width=10, height=5)
  
  return(list(p, p2, p3))
}

pars6$chain_length <- c(0.9,0.9,0.8,0.7,0.9,0.8,0.8,0.8,0.9,0.15,0.5,0.4,0.9,0.9,0.9,0.8)

MCMC_facet <- MCMC_comparator(pars=pars6, strain_interest=NA, chain_length=NA, thinning=0.1)
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
  
  parameter_medians <- setNames(colMedians(output$MCMC_output[start:n_it,]), colnames(output$MCMC_output))
  parameter_means <- colMeans(output$MCMC_output[start:n_it,])
  
  expectation_of_deviance <- mean(deviance_vec)
  deviance_of_expectation <- mean(max(unique(deviance_vec)[100]))
  deviance_of_param_medians <- -2*likelihood_function3(parameters=parameter_medians, knots=output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, population=output$x$population, data_array=weekly_aggregator(output$x$data_array), test_pop=output$x$test_pop, ncores=output$x$ncores, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim)
  deviance_of_param_medians
  deviance_of_param_means <- -2*likelihood_function3(parameters=parameter_means, knots=output$x$knots, vk=1, vg=output$x$vg, max_day=output$x$max_day, population=output$x$population, data_array=weekly_aggregator(output$x$data_array), test_pop=output$x$test_pop, ncores=output$x$ncores, form=output$x$form, symp_delay_lim=output$x$symp_delay_lim)
  deviance_of_param_means
  
  DIC1 = 2*expectation_of_deviance-deviance_of_expectation
  DIC2 = 2*expectation_of_deviance-deviance_of_param_medians
  DIC3 = 2*expectation_of_deviance-deviance_of_param_means
  DIC4 = 0.5*var(likelihoods_vec)+mean(deviance_vec)
  
  return(list(DIC1=DIC1, DIC2=DIC2, DIC3=DIC3, DIC4=DIC4))
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
    
    n_it <- which(is.na(output$MCMC_posteriors))[1]-2
    start <- (1-pars$chain_length[i])*n_it #burnin_calculator(output)
    #if((metrics[[1]] - start) < 10000) start <- metrics[[1]] - 10000
    if(length(output)==1){
      pars[i,c("cov_start", "n_it", "Acceptances", "time_it", "min_ess", "min_ess_par", "DIC")] <- NA
      next
    } 
    ess <- ess_calculator(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], start=start)
    
    pars$n_it[i] <- n_it
    pars$warm_up[i] <- start
    pars$chain_length_abs[i] <- metrics[[1]] - start
    pars$Acceptances[i] <- round(metrics[2],2)
    pars$time_it[i] <- round(metrics[3],2)
    pars$min_ess[i] <- round(min(ess[ess!=0]),1) 
    pars$min_ess_par[i] <- ifelse(is.na(ess[1]), NA, names(which(ess==min(ess[ess!=0]))))
    
    DICs <- DIC_calculator(output, start)
    pars$DIC1[i] <- DICs[[1]]
    pars$DIC2[i] <- DICs[[2]]
    pars$DIC3[i] <- DICs[[3]]
    pars$DIC4[i] <- DICs[[4]]
    pars$posterior[i] <- output$MCMC_posteriors[metrics[1]]
    
  }
  
  pars_new <- pars %>% group_by(strain, corr) %>% mutate(posterior_norm=round(posterior-max(posterior, na.rm=T),0),
                                                         DIC1_norm = DIC1-min(DIC1),
                                                         DIC2_norm = DIC2-min(DIC2),
                                                         DIC3_norm = DIC3-min(DIC3),
                                                         DIC4_norm = DIC4-min(DIC4)) %>% arrange(strain, posterior_norm)
  pars_new %>% arrange(strain, DIC1)
  #pars_new %>% dplyr::select(strain, form, corr, min_ess, DIC, posterior_norm, a, b, vlmax, inc1, inc2, inc3) %>% arrange(strain, DIC) %>% print(n=50)
  
  cat("\n")
  
  p1 <- ggplot(pars_new[,c("strain", "form", "corr", "DIC1_norm")], aes(x=form, y=DIC1_norm)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(~strain, scales="free_y") +
    theme_bw()
  
  p1
  
  p1.1 <- ggplot(pars_new[,c("strain", "form", "corr", "DIC3_norm")], aes(x=form, y=DIC3_norm)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(~strain, scales="free_y") +
    theme_bw()
  
  p1.1
  
  p2 <- ggplot(pars_new, aes(x=form, y=DIC1)) +
    geom_point() +
    facet_grid(strain~corr, scales="free_y") +
    theme_bw()
  
  p2
  
  pars_model <- pars_new[,c("form", "strain", "DIC1_norm", "DIC2_norm", "DIC3_norm", "DIC4_norm")] %>% tidyr::gather(param, value,3:6)
  
  ggplot(pars_model, aes(x=form, y=value, fill=param, group=param)) +
    geom_col(position="dodge") +
    facet_wrap(~strain, scales="free") + theme_bw()
  
  ggsave("images/model_DIC.png", p1, width=10, height=5)
  ggsave("images/model_posteriors.png", p2, width=10, height=5)
  
  
  return(list(DIC=p1, posteriors=p2, metrics=pars_new))
}

metrics <- metrics_comparator(pars=pars6)
metrics[[2]]
metrics[[3]] %>% group_by(strain) %>% arrange(strain, posterior_norm, desc=T)
metrics[[3]] %>% arrange(form, corr, flat)
metrics[[4]]

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

vls_calculator <- function(pars, i, n_samples=50, n_indiv=50){
  strain <- pars$strain[i]
  form <- pars$form[i]
  corr <- pars$corr[i]
  flat <- pars$flat[i]
  chain_length <- pars$chain_length[i]
  
  times <- seq(-14,28,0.1)
  vl_storage <- array(dim=c(length(times), n_indiv, n_samples))
  t_peak <- array(dim=c(n_indiv, n_samples))
  t_peak_end <- array(dim=c(n_indiv, n_samples))
  symp_rel_peak <- array(dim=c(n_indiv, n_samples))
  symp_rel_det <- array(dim=c(n_indiv, n_samples))
  symp_to_end <- array(dim=c(n_indiv, n_samples))
  
  tot_infection_time <- array(dim=c(n_indiv, n_samples))
  
  growth_matrix <- array(dim=c(n_indiv, n_samples))
  decay_matrix <- array(dim=c(n_indiv, n_samples))
  vl_matrix <- array(dim=c(n_indiv, n_samples))
  
  output_raw <- output_finder(strain=strain, form=form, corr=corr, flat=flat)
  if(length(output_raw) == 1) return("skip")
  output <- output_raw$MCMC_output
  
  n_it <- which(is.na(output))[1]-2
  if(chain_length < 1) start <- (1-pars$chain_length[i])*n_it
  else start <- n_it - chain_length
  iterations <- seq(start,n_it,length.out=n_samples)
  
  for(i in 1:n_samples){
    if(n_samples==1) pars_run <- setNames(colMedians(output[start:n_it,]), colnames(output))
    else pars_run <- output[iterations[i],]
    
    vl_max_samples <- rnorm(n=n_indiv, pars_run["vl_max_bar"], pars_run["vl_max_sigma"])
    
    a_samples_vl <- rnorm(n=n_indiv, mean=pars_run["a_bar"], sd=pars_run["a_sigma"])
    b_samples_vl <- rnorm(n=n_indiv, mean=pars_run["b_bar"], sd=pars_run["b_sigma"])
    
    a_samples <- 0.5+exp(a_samples_vl)
    b_samples <- 0.25+exp(b_samples_vl)
    
    for(m in 1:n_indiv){
      vl_storage[,m,i] <- vl_func(a1=a_samples_vl[m], b1=b_samples_vl[m], l1=l_samples[m], t1=times, vlmax1=vl_max_samples[m], vi=1, corr=output_raw$x$corr, flat=output_raw$x$flat)
    }
    
    if(output_raw$x$form=="peak" & output_raw$x$vg==6){
      days = -output_raw$x$symp_delay_lim:output_raw$x$symp_delay_lim
      
      inc_dist = setNames(fGarch::psnorm(days,   mean=pars_run["inc1"], sd=pars_run["inc2"], xi=exp(pars_run["inc3"]))-
                          fGarch::psnorm(days-1, mean=pars_run["inc1"], sd=pars_run["inc2"], xi=exp(pars_run["inc3"])),days)
      
      time_of_peak = (vl_max_samples-vi)/a_samples
      inc_samples = sample(days, size=n_indiv, prob=inc_dist, replace=T)
      t_symptoms = time_of_peak+inc_samples
      accepted <- rep(T,length(inc_samples))
    } 
    if(output_raw$x$form=="thresh" & output_raw$x$vg==6){
      inc_samples <- rnorm(n=n_indiv, mean=pars_run["inc1"], sd=pars_run["inc2"])
      
      vl_thresh <- inc_samples
      accepted <- vl_thresh < vl_max_samples
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+pars_run["inc3"]
    } 
    if(output_raw$x$form=="thresh_peak" & output_raw$x$vg==6){
      inc_samples <- rtruncnorm(n=n_indiv, a=0, mean=pars_run["inc1"], sd=pars_run["inc2"])
      
      vl_thresh <- vl_max_samples-inc_samples
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+pars_run["inc3"]
      accepted <- rep(T,length(inc_samples))
    } 
    if(output_raw$x$form=="thresh_tdist" & output_raw$x$vg==6){
      inc_samples <- exp(rnorm(n=n_indiv, mean=pars_run["inc1"], sd=pars_run["inc2"]))
      
      vl_thresh <- pars_run["inc3"]
      
      t_thresh <- (vl_thresh-vi)/a_samples
      t_symptoms <- t_thresh+inc_samples
      accepted <- rep(T,length(inc_samples))
    } 
    
    t_peak[,i] <- (vl_max_samples-vi)/a_samples
    t_peak_end[,i] <- (vl_max_samples-vi)/b_samples
    
    symp_rel_det[,i] <- t_symptoms
    symp_rel_peak[,i] <- t_symptoms-t_peak[,i]
    tot_infection_time[,i] <- (vl_max_samples-vi)*(a_samples+b_samples)/(a_samples*b_samples)
    symp_to_end[,i] <- tot_infection_time[,i]-t_symptoms
    
    symp_rel_det[!accepted,i] <- NA
    symp_rel_peak[!accepted,i] <- NA
    tot_infection_time[!accepted,i] <- NA
    symp_to_end[!accepted,i] <- NA
    
    growth_matrix[,i] <- a_samples
    decay_matrix[,i] <- b_samples
    vl_matrix[,i] <- vl_max_samples
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
                            vls=rbind(matrixStats::rowQuantiles(median_each_slice, probs=c(0.025, 0.5, 0.975)),
                                      matrixStats::rowQuantiles(lower_each_slice, probs=c(0.025, 0.5, 0.975)),
                                      matrixStats::rowQuantiles(upper_each_slice, probs=c(0.025, 0.5, 0.975)))) %>%
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
  
  plotting_df4 <- reshape2::melt(vl_storage, varnames=c("t", "person", "sample")) %>%
    dplyr::filter(t%%1 == 0) %>%
    mutate(form=form, 
           strain=strain,
           corr=corr, 
           flat=flat) 
  
  plotting_df3 <- reshape2::melt(vl_storage, varnames=c("t", "person", "sample")) %>%
    mutate(ID=paste0(person,sample)) %>% 
    mutate(lead1 =lead(value), lag1 = dplyr::lag(value)) %>%
    group_by(ID) %>%
    dplyr::filter((lead1 != 1 & lag1 == 1) | (lead1 == 1 & lag1 != 1) | t %in% c(-14,28) | value==max(value)) %>%
    mutate(form=form, 
           strain=strain,
           corr=corr, 
           flat=flat) %>%
    dplyr::select(-c(lead1, lag1))
    
  p3 <- ggplot(plotting_df3, aes(x=t, y=value, group=ID)) +
    geom_line(color="blue", alpha=0.2) +
    theme_bw() +
    theme(legend.position="none") +
    geom_line(data=plotting_df2, aes(x=t, y=median, group=1), size=1)
    
  key_vars <- data.frame(form=form, 
                         strain=strain, 
                         corr=corr, 
                         flat=flat,
                         det_to_symp = c(symp_rel_det),
                         symp_to_end = c(symp_to_end),
                         peak_to_symp = c(symp_rel_peak), 
                         det_to_peak = c(t_peak),
                         peak_to_end = c(t_peak_end),
                         total_infection_time = c(tot_infection_time),
                         growth=c(growth_matrix),
                         decay=c(decay_matrix),
                         peak=c(vl_matrix)
                         ) %>% 
    tidyr::gather(param, value, 5:13) 
  
  key_vars_mult <- data.frame(matrixStats::colQuantiles(as.matrix(data.frame(det_to_symp = matrixStats::colQuantiles(symp_rel_det, probs=c(0.025, 0.5, 0.975)),
                              symp_to_end = matrixStats::colQuantiles(symp_to_end, probs=c(0.025, 0.5, 0.975)),
                              peak_to_symp = matrixStats::colQuantiles(symp_rel_peak, probs=c(0.025, 0.5, 0.975)), 
                              det_to_peak = matrixStats::colQuantiles(t_peak, probs=c(0.025, 0.5, 0.975)),
                              peak_to_end = matrixStats::colQuantiles(t_peak_end, probs=c(0.025, 0.5, 0.975)),
                              total_infection_time = matrixStats::colQuantiles(tot_infection_time, probs=c(0.025, 0.5, 0.975)),
                              growth=matrixStats::colQuantiles(growth_matrix, probs=c(0.025, 0.5, 0.975)),
                              decay=matrixStats::colQuantiles(decay_matrix, probs=c(0.025, 0.5, 0.975)),
                              peak=matrixStats::colQuantiles(vl_matrix, probs=c(0.025, 0.5, 0.975)))), probs=c(0.025, 0.5, 0.975))) %>%
    magrittr::set_colnames(c("lower", "median", "upper")) %>%
    tibble::rownames_to_column(var="row_name") %>%
    mutate(param=stringr::str_extract(row_name, "^[^.]+")) %>%
    mutate(type=ifelse(grepl("2.5", row_name), "lower", ifelse(grepl("50", row_name), "median", "upper"))) %>%
    mutate(strain=strain, form=form, corr=corr, flat=flat) %>%
    dplyr::select(strain, form, corr, flat, param, type, lower, median, upper)
  
  p4 <- ggplot(key_vars, aes(x=param, y=value)) +
    geom_boxplot() +
    #facet_wrap() +
    theme_bw()
  
  #ggsave(paste0(strain,"/vl_trajectories/",strain,"_",form,"_corr=",corr,"_flat=",flat), p)
  
  return(list(plotting_df, plotting_df2, plotting_df3, plotting_df4, key_vars, key_vars_mult))
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

mapping_func <- function(out){
  mapped_data <- out %>%
    mutate(param_plotting = case_when(
      param == "det_to_symp" ~ "Detection to symptom onset",
      param == "symp_to_end" ~ "Symptom onset to end of infection",
      param == "peak_to_symp" ~ "Time of peak to symptom onset",
      param == "det_to_peak" ~ "Detection to peak viral load",
      param == "peak_to_end" ~ "Peak viral load to end of infection",
      param == "total_infection_time" ~ "Total infection time",
      TRUE ~ as.character(param)  # If none of the above conditions match, keep original value
    ))  %>%
    mutate(strain_plotting = case_when(
      strain == "omicron" ~ "Omicron",
      strain == "alpha" ~ "Alpha",
      strain == "delta" ~ "Delta",
      strain == "wt" ~ "Wild type",
      TRUE ~ as.character(param)  # If none of the above conditions match, keep original value
    ))
  return(mapped_data)
}

symp_plotter <- function(out6, params_interest, legend="none", y_label=NULL){
  
  out6_adj <- mapping_func(out6) %>% mutate(strain_plotting=factor(strain_plotting, levels=c("Wild type", "Alpha", "Delta", "Omicron")))
  
  p6 <- ggplot(out6_adj %>% dplyr::filter(param %in% params_interest) %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(x=strain_plotting, color=type)) +
    geom_point(aes(y=median)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1) +
    theme_bw() +
    facet_wrap(~param_plotting, scales="free_x") +
    coord_flip() +
    theme(legend.position = legend) +
    labs(y=y_label, x=NULL) +
    theme(panel.grid=element_blank())
  
  out7_raw <- out6_adj %>% group_by(strain, param, form, corr, flat) %>% tidyr::gather(CI, value, 7:9) %>% tidyr::spread(type, value)
  out7_median <- out7_raw %>% dplyr::filter(CI=="median", param %in% params_interest)
  out7_lower <- out7_raw %>% dplyr::filter(CI=="lower", param %in% params_interest)
  out7_upper <- out7_raw %>% dplyr::filter(CI=="upper", param %in% params_interest)
  
  p7 <- ggplot(out7_median, aes(x=strain_plotting, color=strain_plotting, fill=strain_plotting)) +
    geom_point(aes(y=median, size=0.4)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5) +
    geom_crossbar(data=mapping_func(out6_adj) %>% dplyr::filter(param %in% params_interest), aes(ymin=lower, y=median, ymax=upper, alpha=0.1), width=0.5, color="transparent") +
    theme_bw() +
    facet_wrap(~param_plotting, scales="free_x") +
    coord_flip() +
    theme(legend.position = legend) +
    labs(y=y_label, x=NULL) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(color="white", fill="white", size=0.5, linetype="solid"),
          strip.text.x = element_text(size=12, hjust = 0.5, margin=margin(l=2, b=2), color="black"),
          panel.border=element_rect(color="white"),
          axis.title.x = element_text(size=12, color="black", margin=margin(t=10)),
          axis.title.y = element_text(size=12, color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black")) + theme(plot.margin = margin(0.25,0,0.25,0, "cm"))
  
  p7
  
  return(list(p6,p7))
}

for(i in 1:nrow(pars)) hist_plotter(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i], chain_length=3000, n_samples=1000)

vls_comparator <- function(pars, tag="", n_samples=50, n_indiv=50){
  out1 <- data.frame(t=numeric(), form=character(), strain=character(), corr=logical(), flat=character(), lower=numeric(), median=numeric(), upper=numeric())
  out2 <- data.frame(t=numeric(), person=numeric(), sample=numeric(), value=numeric(), ID=character(), form=character(), strain=character(), corr=numeric(), flat=numeric())
  out3 <- data.frame(t=numeric(), person=numeric(), sample=numeric(), value=numeric(), form=character(), strain=character(), corr=numeric(), flat=numeric())
  out_symp <- data.frame(form=character(), strain=character(), corr=numeric(), flat=numeric(), param=character(), value=numeric())
  out6 <- data.frame(strain=character(), form=character(), corr=numeric(), flat=numeric(), param=character(), type=character(), lower=numeric(), median=numeric(), upper=numeric())
  
  out1_med <- data.frame(t=numeric(), form=character(), strain=character(), corr=logical(), flat=character(), lower=numeric(), median=numeric(), upper=numeric())
  out2_med <- data.frame(t=numeric(), person=numeric(), sample=numeric(), value=numeric(), ID=character(), form=character(), strain=character(), corr=numeric(), flat=numeric())
  out3_med <- data.frame(t=numeric(), person=numeric(), sample=numeric(), value=numeric(), form=character(), strain=character(), corr=numeric(), flat=numeric())
  out_symp_med <- data.frame(form=character(), strain=character(), corr=numeric(), flat=numeric(), param=character(), value=numeric())
  out6_med <- data.frame(strain=character(), form=character(), corr=numeric(), flat=numeric(), param=character(), type=character(), lower=numeric(), median=numeric(), upper=numeric())
  
  next_line1 <- 1
  next_line2 <- 1
  next_line3 <- 1
  next_line4 <- 1
  next_line6 <- 1
  
  next_line1_med <- 1
  next_line2_med <- 1
  next_line3_med <- 1
  next_line4_med <- 1
  next_line6_med <- 1
  
  for(i in 1:nrow(pars)){
    vls <- vls_calculator(pars=pars, i=i, n_indiv=n_indiv, n_samples=n_samples)
    vls_med <- vls_calculator(pars=pars, i=i, n_indiv=1000, n_samples=1)
    if(length(vls) != 1){
      out1[next_line1:(next_line1+nrow(vls[[2]])-1),] <- vls[[2]]
      out2[next_line2:(next_line2+nrow(vls[[3]])-1),] <- vls[[3]]
      out3[next_line3:(next_line3+nrow(vls[[4]])-1),] <- vls[[4]]
      out_symp[next_line4:(next_line4+nrow(vls[[5]])-1),] <- vls[[5]]
      out6[next_line6:(next_line6+nrow(vls[[6]])-1),] <- vls[[6]]
      
      out1_med[next_line1_med:(next_line1_med+nrow(vls_med[[2]])-1),] <- vls_med[[2]]
      out2_med[next_line2_med:(next_line2_med+nrow(vls_med[[3]])-1),] <- vls_med[[3]]
      out3_med[next_line3_med:(next_line3_med+nrow(vls_med[[4]])-1),] <- vls_med[[4]]
      out_symp_med[next_line4_med:(next_line4_med+nrow(vls_med[[5]])-1),] <- vls_med[[5]]
      out6_med[next_line6_med:(next_line6_med+nrow(vls_med[[6]])-1),] <- vls_med[[6]]
      
      next_line1 <- next_line1+nrow(vls[[2]])
      next_line2 <- next_line2+nrow(vls[[3]])
      next_line3 <- next_line3+nrow(vls[[4]])
      next_line4 <- next_line4+nrow(vls[[5]])
      next_line6 <- next_line6+nrow(vls[[6]])
      
      next_line1_med <- next_line1_med+nrow(vls_med[[2]])
      next_line2_med <- next_line2_med+nrow(vls_med[[3]])
      next_line3_med <- next_line3_med+nrow(vls_med[[4]])
      next_line4_med <- next_line4_med+nrow(vls_med[[5]])
      next_line6_med <- next_line4_med+nrow(vls_med[[6]])
    } 
  }
  
  p <- ggplot(out1 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(x=t)) +
    geom_line(aes(y=median)) +
    facet_nested(form ~ strain, scales="fixed", labeller = label_wrap_gen(multi_line=FALSE)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color="white") +
    theme_bw() +
    scale_x_continuous(breaks=seq(-14,28,7)) +
    scale_y_continuous(breaks=seq(2,10,2)) +
    theme(panel.grid = element_blank())
  
  p
  
  p1 <- ggplot(out1_med %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(x=t)) +
    geom_ribbon(data=out1 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(ymin=lower, ymax=upper), alpha=0.1, fill="blue") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.1, fill="red") +
    geom_line(aes(y=median)) +
    facet_nested(form ~ strain, scales="fixed", labeller = label_wrap_gen(multi_line=FALSE)) +
    theme_bw() +
    scale_x_continuous(breaks=seq(-14,28,7)) +
    scale_y_continuous(breaks=seq(2,10,2)) +
    theme(panel.grid = element_blank())
  
  p1
  
  p2 <- ggplot(out2 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))) %>%
                 mutate(strain=recode(strain, "wt"="Wild type",
                                      "alpha"="Alpha", "delta"="Delta", "omicron"="Omicron")), aes(x=t, y=value, group=ID)) +
    geom_line(aes(color=strain), alpha=0.03) +
    theme_bw() +
    theme(legend.position="none") +
    facet_grid(~strain, scales="fixed", labeller = label_wrap_gen(multi_line=FALSE)) +
    geom_line(data=out1 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))) %>%
                mutate(strain=recode(strain, "wt"="Wild type",
                                     "alpha"="Alpha", "delta"="Delta", "omicron"="Omicron")), aes(x=t, y=median, group=1), size=1) +
    scale_x_continuous(breaks=seq(-14,28,7)) + 
    scale_y_continuous(breaks=seq(2,10,2), name = expression(paste("Viral Load (", log[10]," RNA copies mL"^{-1},")", sep = ""))) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(color="white", fill="white", size=0.5, linetype="solid"),
          strip.text.x = element_text(size=12, hjust = 0, margin=margin(l=2, b=2), color="black"),
          panel.border=element_blank(),
          axis.title.x = element_text(size=12, color="black", margin=margin(t=10)),
          axis.title.y = element_text(size=12, color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black")) +
    labs(x="Time from peak (days)") 
  
  p2
  
  p2.2 <- ggplot(out3 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron")), t=factor(t, levels=c(-14:28))), aes(x=t, y=value)) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw() +
    facet_nested(form ~ strain, scales="fixed", labeller = label_wrap_gen(multi_line=FALSE)) +
    scale_x_discrete(breaks=seq(-14,28,7)) + 
    scale_y_continuous(breaks=seq(2,10,2)) +
    theme(panel.grid = element_blank()) 
  
  p2.2
  
  out4 <- out3 %>% mutate(t=factor(t, levels=c(-14:28))) %>% group_by(strain, form, t) %>% summarise(lower=quantile(value, probs=c(0.025)),
                                                                                             median=quantile(value, probs=c(0.5)), 
                                                                                             upper=quantile(value, probs=0.975))
  
  p2.3 <- ggplot(out4 %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(x=t)) +
    geom_pointrange(aes(y=median, ymin=lower, ymax=upper)) +
    theme_bw() +
    facet_nested(form ~ strain, scales="fixed", labeller = label_wrap_gen(multi_line=FALSE)) +
    scale_x_discrete(breaks = levels(out4$t)[c(T, rep(F, 6))]) +
    scale_y_continuous(breaks=seq(2,10,2)) +
    theme(panel.grid = element_blank())
  
  p2.3
  
  removed5 <- out_symp %>% group_by(form, strain, corr, flat, param) 
  
  p5 <- ggplot(out6 %>% dplyr::filter(param %in% c("growth", "decay", "peak")) %>% mutate(strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))), aes(x=strain, color=type)) +
    geom_point(aes(y=median)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1) +
    theme_bw() +
    facet_wrap(~param, scales="free") 
  
  out6 %>% tidyr::pivot_wider(names_from="type", values_from=c("lower", "median", "upper")) %>%
    mutate(lower=paste0(round(median_lower,1), " (", round(lower_lower,1), "-", round(upper_lower,1), ")"),
           median=paste0(round(median_median,1), " (", round(lower_median,1), "-", round(upper_median,1), ")"),
           upper=paste0(round(median_upper,1), " (", round(lower_upper,1), "-", round(upper_upper,1), ")")) %>%
    dplyr::select(strain, param, lower, median, upper)
  
  det_peak_end <- ggarrange(symp_plotter(out6, c("det_to_peak","peak_to_end"), legend="none")[[2]])
  tot_inf <- ggarrange(symp_plotter(out6, c("total_infection_time"), legend="none")[[2]])
  peak_symp <- ggarrange(symp_plotter(out6, c("peak_to_symp"), legend="none")[[2]] + theme(plot.margin = margin(0,6,0,6, "cm")))
  det_symp_end <- ggarrange(symp_plotter(out6, c("det_to_symp","symp_to_end"), legend="none", y_label="days")[[2]])
  
  p7 <- ggarrange(det_peak_end, tot_inf, peak_symp, det_symp_end, nrow=4, heights=c(1,1,1,1), widths=c(1,1,0.5,1))
  
  vk_CrIs <- removed5 %>% reframe(value=quantile(value, probs=c(0.025, 0.5, 0.975))) %>%
    mutate(bound=rep(c("lower", "median", "upper"), nrow(.)/3)) %>%
    tidyr::pivot_wider(
      id_cols = c("form", "strain", "corr", "flat"),
      names_from = c("param", "bound"),
      values_from = c("value")
    )
  
  ggsave(paste0("images/viral_load_comparison_aggregate_",tag,".png"), p)
  ggsave(paste0("images/viral_load_comparison_samples_",tag,".png"), p2)
  ggsave(paste0("images/viral_load_comparison_samples_box_",tag,".png"), p2.2)
  ggsave(paste0("images/viral_load_comparison_samples_pointrange_",tag,".png"), p2.3)
  ggsave(paste0("images/delay_comparisons_2_",tag,".png"), p7, width=8, height=12)
  ggsave(paste0("images/vk_comparisons_",tag,".png"), p4)
  
  print(p)
  
  return(list(out1, vk_CrIs))
}

vl_comparison <- vls_comparator(pars=pars6 %>% dplyr::filter(form=="thresh_peak"), 
                                tag="thresh_peak", n_samples=50, n_indiv=50)

vl_comparison[[2]]

prior_vs_posterior <- function(pars, i, n_samples=1000){
  output <- output_finder(strain=pars$strain[i], form=pars$form[i], corr=pars$corr[i], flat=pars$flat[i])
  
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  start <- (1-pars$chain_length[i])*n_it
  
  samples_prior <- mvtnorm::rmvnorm(n=n_samples, mean=output$x$prior_mean, sigma=output$x$prior_cov_final)
  xa <- c(seq(-5,1.5,length.out=1001))
  xb <- c(seq(-5,0,length.out=1001))
  x2 <- c(seq(3,12,length.out=1001))
  
  prior_array_a <- matrix(nrow=n_samples, ncol=length(xa))
  post_array_a <- matrix(nrow=n_samples, ncol=length(xa))
  prior_array_b <- matrix(nrow=n_samples, ncol=length(xb))
  post_array_b <- matrix(nrow=n_samples, ncol=length(xb))
  prior_array_vl <- matrix(nrow=n_samples, ncol=length(x2))
  post_array_vl <- matrix(nrow=n_samples, ncol=length(x2))
  iterations <- round(seq(start, n_it, length.out=n_samples))
  
  for(i in 1:length(x)){
    prior_array_a[,i] <- dnorm(x=xa[i], mean=samples_prior[,1], sd=samples_prior[,2])
    post_array_a[,i] <- dnorm(x=xa[i], 
                               mean=output$MCMC_output[iterations,1], 
                               sd=output$MCMC_output[iterations,2])
    prior_array_b[,i] <- dnorm(x=xb[i], mean=samples_prior[,3], sd=samples_prior[,4])
    post_array_b[,i] <- dnorm(x=xb[i], 
                              mean=output$MCMC_output[iterations,3], 
                              sd=output$MCMC_output[iterations,4])
    
    prior_array_vl[,i] <- dnorm(x=x2[i], mean=samples_prior[,5], sd=samples_prior[,6])
    post_array_vl[,i] <- dnorm(x=x2[i], 
                              mean=output$MCMC_output[iterations,5], 
                              sd=output$MCMC_output[iterations,6])
  }
  
  df <- data.frame(strain=output$x$strain,
                   xa=xa,
                   xb=xb,
                   x2=x2,
                   matrixStats::colQuantiles(prior_array_a, probs=c(0.025, 0.5, 0.975)),
                   matrixStats::colQuantiles(post_array_a, probs=c(0.025, 0.5, 0.975)),
                   matrixStats::colQuantiles(prior_array_b, probs=c(0.025, 0.5, 0.975)),
                   matrixStats::colQuantiles(post_array_b, probs=c(0.025, 0.5, 0.975)),
                   matrixStats::colQuantiles(prior_array_vl, probs=c(0.025, 0.5, 0.975)),
                   matrixStats::colQuantiles(post_array_vl, probs=c(0.025, 0.5, 0.975))) %>%
    magrittr::set_colnames(c("strain", "xa", "xb", "x2", "a_lower_Prior", "a_median_Prior", "a_upper_Prior",
                             "a_lower_Posterior", "a_median_Posterior", "a_upper_Posterior",
                             "b_lower_Prior", "b_median_Prior", "b_upper_Prior",
                             "b_lower_Posterior", "b_median_Posterior", "b_upper_Posterior",
                             "vl_lower_Prior", "vl_median_Prior", "vl_upper_Prior",
                             "vl_lower_Posterior", "vl_median_Posterior", "vl_upper_Posterior")) %>%
    tidyr::gather(param, value, 5:ncol(.)) %>%
    mutate(vk_param = sub("_.*", "", param),
           type = sub("^[^_]*_(.*?)_.*$", "\\1", param),
           pp = sub("^.+_", "", param)) %>%
    dplyr::select(-param) %>%
    tidyr::pivot_wider(., names_from=type, values_from=value) %>%
    mutate(implied = ifelse(vk_param=="a", 0.5+exp(xa), ifelse(vk_param=="b", 0.25+exp(xb), x2)),
           vk=ifelse(vk_param=="a", "Growth rate", ifelse(vk_param=="b", "Decay rate", "Peak viral load"))) %>%
    mutate(vk=factor(vk, levels=c("Growth rate", "Decay rate", "Peak viral load")))
  
  ab_pp <- ggplot(df %>% dplyr::filter(vk != "Peak viral load"), aes(x=implied)) +
    geom_line(aes(y=median, color=pp)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=pp), alpha=0.2) +
    facet_wrap(~vk, scale="free") +
    theme_bw() +
    #lims(x=c(0.5, 5)) +
    theme(panel.grid=element_blank(),
          legend.title = element_blank(),
          legend.position = "none") +
          #strip.background = element_rect(color="white", fill="white", size=0.5, linetype="solid"),
          #strip.text.x = element_text(size=12, hjust = 0.5, margin=margin(l=2, b=2), color="black"),
          #panel.border=element_rect(color="black")) +
    labs(y="Probability density", x=expression("Viral load (" * log[10] * " RNA copies mL"^{-1} * "day"^{-1} * ")"))
    
  vl_pp <- ggplot(df %>% dplyr::filter(vk_param == "vl"), aes(x=implied)) +
    geom_line(aes(y=median, color=pp)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=pp), alpha=0.2) +
    facet_wrap(~vk, scale="free") +
    theme_bw() +
    #lims(x=c(0.5, 5)) +
    theme(panel.grid=element_blank(),
          legend.position = "none") +
    labs(y="Probability density", x=expression("Viral load (" * log[10] * " RNA copies mL"^{-1} * ")"))
  
  p_pp <- ggarrange(ab_pp, vl_pp, nrow=2)
  
  ggsave(paste0("images/",output$x$strain,"_vk_prior_posterior_",".png"), p_pp, width=8, height=6)
  
  
}

for(i in 1:nrow(pars)) prior_vs_posterior(pars, i)

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
