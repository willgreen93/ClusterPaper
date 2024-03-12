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

setwd("Q:/ClusterPaper")

Sys.setenv(PATH = paste("C:/rtools40/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")

source("ct_symptoms_key_functions.R")

input_generator <- function(n_iterations=5e6, strain="all", test_pop=1e7, ncores, proposal_type, vl_growth_type="lognorm", vl_kinetics_type="differentiable", vl_truncation="truncated",
                            delay_dist="skew_normal", cov_matrix=NULL, population=6.7e7, form="peak", tag, test_form="empirical", ignored_pars=NA, spline_burnin=NA, parameters=NA, 
                            day_agg="weekly", cov_start=200000, ct_range=c(10,39), max_test_day=6, gene="N"){
  
  #if(gene=="ORF1ab") data_array <- readRDS(paste0("data/data_array_vl_",strain,".rds"))[1:((ct_range[2]-ct_range[1])+1),,1:(max_test_day+1)]
  if(gene=="N") data_array <- readRDS(paste0("data/data_array_vl_",strain,"_N",".rds"))[,,1:(max_test_day+1)]
  
  max_day <- dim(data_array)[2]
  ndays <- dim(data_array)[2]
  if(ndays < 490) knots <- seq(0,ndays,14)[-1]/ndays
  if(ndays > 490) knots <- c(seq(0,490,14)[-1], seq(497,ndays,7))/ndays
  
  if(strain=="omicron") knots <- c(seq(0,ndays-7,7))[-1]/ndays
  theta_names_inc <- paste0("theta", seq(1:(length(knots)+1)))
  theta_names_peak <- paste0("theta", seq(1:(length(knots)+2)))
  
  ## vl parameters
  priors <- readRDS(paste0("Q:/StanPriors/hierarchical_2.13b_combined_",vl_growth_type,"_ab_",vl_kinetics_type,".rds"))
  vl_prior_means <- priors$mean
  vl_prior_sds <- priors$sd
  
  vl_prior_means['l_bar'] <- 2
  vl_prior_means['l_sigma'] <- 1
  
  vl_prior_sds['l_bar'] <- 4
  vl_prior_sds['l_sigma'] <- 1
  
  if(test_form=="empirical"){
    test_prior_means <- (apply(data_array, 3, sum)/sum(data_array))[-dim(data_array)[3]]
    #test_prior_means <- rep(1/(max_test_day+1), max_test_day)
    test_prior_sds <-   rep(0.5, max_test_day)
    
    names(test_prior_means) <- paste0("test",seq(0,max_test_day-1))  
    names(test_prior_sds) <- paste0("test",seq(0,max_test_day-1))  
    
    test_rdisp_prior_means <- c(test_prior_means, rdisp=1) # remember prior for inc is for the mean and sd
    test_rdisp_prior_sds <-   c(test_prior_sds, rdisp=10)
  }
  
  ## incubation period and peak to symp delay parameter
  if(form == "incidence"){
    incidence_inc_prior_means <- c(inc1=2, inc2=1, inc3=0)
    incidence_inc_prior_sds <-   c(inc1=2, inc2=1, inc3=2)
   
    incidence_spline_prior_means <- setNames(c(log(apply(data_array, 2, sum)[knots*ndays]/population),-Inf), theta_names_inc)
    incidence_spline_prior_means[incidence_spline_prior_means==-Inf] <- -14
    incidence_spline_prior_sds <- abs(incidence_spline_prior_means)
    
    incidence_spline_multiplier_mean <- c(multiplier=1)
    incidence_spline_multiplier_sd <- c(multiplier=10)
    
    prior_mean <- c(vl_prior_means, incidence_inc_prior_means, test_rdisp_prior_means, incidence_spline_prior_means, incidence_spline_multiplier_mean)
    prior_sds <- c(vl_prior_sds, incidence_inc_prior_sds, test_rdisp_prior_means, incidence_spline_prior_sds, incidence_spline_multiplier_sd)
    
    ignored_spline_par <- length(prior_mean) - ceiling(length(theta_names_peak)/2)
    
    if(is.na(parameters)[1]){
      parameters <- c(prior_mean)
      parameters[-ignored_spline_par] <- rnorm(n=length(parameters)-1, mean=parameters[-ignored_spline_par], sd=abs(parameters[-ignored_spline_par])*0.05)
    }  
  }
  
  if(form=="peak"){
    peak_symp_prior_means <- c(inc1=-2, inc2=1, inc3=0)
    peak_symp_prior_sds <- c(inc1=0.5, inc2=0.3, inc3=3)
    
    peak_spline_prior_means <- setNames(c(-Inf,log(apply(data_array, 2, sum)[knots*ndays]/population),-Inf), theta_names_peak)
    peak_spline_prior_means[peak_spline_prior_means==-Inf] <- -14
    peak_spline_prior_sds <- abs(peak_spline_prior_means)
    
    peak_spline_multiplier_mean <- c(multiplier=1)
    peak_spline_multiplier_sd <- c(multiplier=10)
    
    prior_mean <- c(vl_prior_means, peak_symp_prior_means, test_rdisp_prior_means, peak_spline_prior_means, peak_spline_multiplier_mean)
    prior_sds <- c(vl_prior_sds, peak_symp_prior_sds, test_rdisp_prior_means, peak_spline_prior_sds, peak_spline_multiplier_sd)
    
    ignored_spline_par <- length(prior_mean) - ceiling(length(theta_names_peak)/2)
    
    if(is.na(parameters)[1]){
      parameters <- c(prior_mean)
      parameters[-ignored_spline_par] <- rnorm(n=length(parameters)-1, mean=parameters[-ignored_spline_par], sd=abs(parameters[-ignored_spline_par])*0.05)
    } 
  }
  
  MCMC_sds <-abs(parameters)/200
  if(delay_dist=="skew_normal") MCMC_sds['inc3'] <- 0.05
  
  max_inf <- 43
  if(form=="peak") symp_delay_lim <- 14
  if(form=="incidence") symp_delay_lim <- 21
  
  #setup <- ignored_pars_func(parameters, ignored_pars, delay_dist)
  #ignored_pars <- setup$ignored_pars
  #parameters <- setup$parameters
  
  input <- list(n_iterations=1000000,
                parameters=parameters,
                vl_growth_type=vl_growth_type, 
                vl_kinetics_type=vl_kinetics_type,
                vl_truncation=vl_truncation,
                delay_dist=delay_dist,
                knots=knots,
                data_array=data_array,
                test_pop=1e7,
                MCMC_sds=MCMC_sds,
                prior_mean=prior_mean,
                prior_sds=prior_sds,
                ncores=ncores,
                proposal_type=proposal_type,
                cov_matrix = cov_matrix,
                population = 6.7e7,
                form=form,
                tag=tag,
                test_form=test_form,
                ignored_pars = ignored_pars,
                ignored_spline_par = length(parameters)-floor(length(theta_names_peak)/2),
                spline_burnin = spline_burnin,
                strain=strain,
                max_inf=max_inf,
                symp_delay_lim=symp_delay_lim,
                day_agg=day_agg,
                max_day=max_day,
                cov_start=cov_start,
                ct_range=ct_range,
                gene=gene,
                n_parameters=length(parameters))
  
  return(input)
}

ignored_pars_func <- function(parameters, ignored_pars, delay_dist){
  ignored_par1 <- which(names(parameters)=="l_sigma")
  parameters['l_sigma'] <- 0

  if(delay_dist=="normal"){
    ignored_pars2 <- which(names(parameters)=="inc3")
    parameters['inc3'] <- 0
    if(is.na(ignored_pars)) ignored_pars_final <- c(ignored_par1, ignored_pars2)
    else ignored_pars_final <- c(ignored_par1, ignored_pars2, ignored_pars)
  }
  else ignored_pars_final <- ifelse(is.na(ignored_pars),c(ignored_par1),c(ignored_par1, ignored_pars))
  
  return(list(parameters=parameters,ignored_pars=ignored_pars_final))
}

input_cov_generator <- function(file_name, ncores=16){
  file <- readRDS(paste0(file_name,".rds"))
  
  start <- which(is.na(file$MCMC_posteriors))[1]-2
  
  nroll <- 10000
  rollmean <- rollmean(file$MCMC_posteriors[1:start], nroll)
  
  #plot(c(file$MCMC_posteriors[1:start]), type="l")
  #lines(c(rep(NA,nroll),rollmean[1:(start-nroll)]), col="red")
  
  #greater <- which(rollmean[1:(start-nroll)] > c(file$MCMC_posteriors[(nroll+1):start]))
  #burnin <- greater[greater>file$x$spline_burnin][10000]
  
  burnin <- which(rollmean > rollmean[start-nroll])[1]
  
  print(paste0("burnin = ", burnin), quote=F)
  print(paste0("iterations = ", start-burnin), quote=F)
  
  plot(file$MCMC_posteriors[burnin:start], type="l")
  
  input <- list(n_iterations=1000000,
                parameters=file$MCMC_output[start,],
                knots=file$x$knots,
                data_array=file$x$data_array,
                test_pop=file$x$test_pop,
                MCMC_sds=file$x$MCMC_sds,
                prior_mean=file$x$prior_mean,
                prior_sds=file$x$prior_sds,
                ncores=ifelse(is.null(ncores),file$x$ncores,ncores),
                proposal_type="Cov",
                cov_matrix = cov(file$MCMC_output[burnin:start,])*(2.38^2)/length(file$x$parameters),
                population = 6.7e7,
                form=file$x$form,
                tag=file$x$tag,
                test_form=file$x$test_form,
                ignored_pars = NA,
                ignored_spline_par = NA,
                spline_burnin = NA,
                strain=file$x$strain,
                max_inf=file$x$max_inf,
                symp_delay_lim=file$x$symp_delay_lim,
                day_agg=file$x$day_agg,
                max_day=file$x$max_day)
  
}

drat:::add("mrc-ide")

options(
  didehpc.username = "wg4618",   # change to your username 
  didehpc.home = "Q:/")    # can be any letter 

didehpc:::didehpc_config(credentials=list(username="wg4618"), home = "Q:/", r_version="4.2.1")

context::context_log_start()

root <- "context"
setwd("Q:/ClusterPaper")

ctx <- context::context_save(root, packages=c("dplyr", "ggplot2", "zoo", "matrixStats", "Rcpp", "dqrng", "truncnorm", "mvtnorm", "BH", "sitmo", "fGarch"),
                             sources=c("ct_symptoms_MCMC_cluster_form.R"))


obj3 <- didehpc::queue_didehpc(ctx, config = didehpc:::didehpc_config(cores=24, credentials=list(username="wg4618"), home = "Q:/", r_version="4.1.3", cluster="fi--didemrchnb", template="24Core"))

input <- input_generator(n_iterations=5e6, strain="all", test_pop=1e7, ncores=24, proposal_type="MH", vl_growth_type="lognorm", vl_kinetics_type="non-differentiable", vl_truncation="truncated",
                        delay_dist="skew_normal", cov_matrix=NULL, population=6.7e7, form="peak", tag="test1", test_form="empirical", ignored_pars=NA, spline_burnin=1000, parameters=NA, 
                        day_agg="weekly", cov_start=200000, ct_range=c(10,39), max_test_day=6, gene="N")

M3 <- obj3$enqueue(MCMC_function(x=input))
M3$log()


baseR.rollmean <- function(dat, window) {
  n <- length(dat)
  y <- dat[window:n] - dat[c(1, 1:(n-window))]
  y[1] <- sum(dat[1:window])
  return(cumsum(y) / window)
}

output_combine_function <- function(file_name1, file_name2){
  file1 <- readRDS(paste0(file_name1,".rds"))
  file2 <- readRDS(paste0(file_name2,".rds"))
  
  n_it1 <- which(is.na(file1$MCMC_posteriors))[1]-2
  n_it2 <- which(is.na(file2$MCMC_posteriors))[1]-2
  
  MCMC_posteriors <- c(file1$MCMC_posteriors[1:n_it1],file2$MCMC_posteriors[1:n_it2])
  MCMC_output <- rbind(file1$MCMC_output[1:n_it1,], file2$MCMC_output[1:n_it2,])
  
  plot(MCMC_posteriors, type="l")
  
  list <- list(MCMC_output = MCMC_output, 
               MCMC_posteriors = MCMC_posteriors,
               n_it1 = n_it1, 
               n_it2 = n_it2)
  
  return(list)
  
}

chain_convergence_checker <- function(strain, tag1, tag2, param_interest=NULL, warmup=NULL){
  file1 <- output_extractor(strain, tag1)
  file2 <- output_extractor(strain, tag2)
  
  n_it1 <- which(is.na(file1$MCMC_posteriors))[1]-2
  n_it2 <- pmin(which(is.na(file2$MCMC_posteriors))[1]-2,file2$x$n_iterations, na.rm=T)

  df <- data.frame(iter=0:file1$x$n_iterations, chain1=file1$MCMC_posteriors, chain2=file2$MCMC_posteriors) %>% filter(iter<max(n_it1, n_it2), iter>100, iter%% 100 == 0) %>%
    tidyr::gather(chain, value, 2:3)
  
  ggplot(data=df, aes(x=iter, y=value)) +
    geom_line(aes(color=chain)) +
    theme_bw() +
    ggtitle(label=paste0(file1$x$strain))
  
  df2 <- rbind(data.frame(chain=1,iter=1:n_it1,file1$MCMC_output[1:n_it1,], posterior=file1$MCMC_posteriors[1:n_it1]),
               data.frame(chain=2,iter=1:n_it2,file2$MCMC_output[1:n_it2,], posterior=file2$MCMC_posteriors[1:n_it2])) %>%
    filter(iter %% 100 == 0) %>%
    tidyr::gather(param, value, 3:(ncol(.))) %>%
    mutate(chain=factor(chain)) %>%
    mutate(param=factor(param, levels=c(colnames(file1$MCMC_output),"posterior"))) 
  
  if(is.null(param_interest)==F) df2 <- df2 %>% filter(param==param_interest)
  if(is.null(warmup)==F) df2 <- df2 %>% filter(iter > warmup)
  
  p2 <- ggplot(data=df2, aes(x=iter, y=value, color=chain)) +
    geom_line() +
    facet_wrap(~param, scales="free_y") +
    theme_bw() 
  
  return(p2)
}

chain_convergence_checker(strain="alpha", tag1="2023-07-19_alpha_408000_N_incidence_Cov_ignored_NA_lognorm_non-differentiable_truncated_skew_normal_",
                          tag2="2023-07-11_alpha_1000000_N_incidence_Cov_ignored_NA_lognorm_non-differentiable_truncated_skew_normal_", 
                          param_interest=NA, warmup=NULL)

output_extractor <- function(strain, tag){
  file_list <- paste0(strain,"/",list.files(paste0(strain,"/")))
  tagged_files <- file_list[which(grepl(paste0(tag),file_list))]
  
  file_info <- file.info(tagged_files)
  recent_tagged_file <- tagged_files[which(file_info$mtime == max(file_info$mtime))]
  
  print(recent_tagged_file, quote=F)
  
  output <- readRDS(paste0(recent_tagged_file))
  return(output)
}

burnin_calculator <- function(output, nroll, n_it){
  rollmean <- baseR.rollmean(dat=output$MCMC_posteriors[1:n_it], window=nroll)
  burnin <- which(rollmean > rollmean[n_it-nroll])[1]
  return(burnin)
}

all_spline_extractor <- function(input, all_tag){
  output <- output_extractor("all", all_tag)
  n_it <- which(is.na(output$MCMC_posteriors))[1] - 2
  all_knots <- c(NA,output$x$knots*output$x$max_day,NA)
  strain <- input$strain
  
  warmup <- burnin_calculator(output, 10000)
  
  theta_start <- which(colnames(output$MCMC_output)=="theta1")
  all_thetas <- colMeans(output$MCMC_output[warmup:n_it,])[theta_start:(ncol(output$MCMC_output)-1)]
  
  strain_dates <- readRDS("strain_dates.rds")
  
  strain_knots <- input$knots+as.numeric(strain_dates$start[strain]-strain_dates$start['all'])
  
  closest_knot <- sapply(strain_knots, function(x) which.min(abs(all_knots - x)))
  thetas <- setNames(all_thetas[closest_knot], paste0("theta",seq(2,(length(strain_knots)+1))))
  
  input$parameters[names(thetas)] <- thetas
  input$parameters['multiplier'] <- mean(output$MCMC_output[(warmup:n_it),'multiplier'])
  
  return(input)
  
}

plotting_func_loop <- function(tag_list, strain, n_samples=1000){
  plotting_func(strain=strain, 
                tag="N_incidence_Cov_ignored_NA_lognorm_non-differentiable_non_truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_incidence_Cov_ignored_NA_lognorm_non-differentiable_truncated",
                warmup_raw=NULL, n_samples=1000)
  
  #plotting_func(strain=strain, 
  #              tag="N_incidence_Cov_ignored_8_lognorm_non-differentiable_non_truncated",
  #              warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_incidence_Cov_ignored_8_lognorm_non-differentiable_truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_incidence_Cov_ignored_NA_lognorm_differentiable_truncated",
                warmup_raw=NULL, n_samples=10)
  
  plotting_func(strain=strain, 
                tag="626000_N_incidence_Cov_ignored_7,8_lognorm_non-differentiable_truncated_skew_normal_",
                warmup_raw=NULL, n_samples=1000)
  
  
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_NA_lognorm_non-differentiable_non-truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_NA_lognorm_non-differentiable_truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_8_lognorm_non-differentiable_non-truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_8_lognorm_non-differentiable_truncated",
                warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="2023-08-11_all_1000000_N_peak_Cov_ignored_7,8_lognorm_non-differentiable_truncated_skew_normal_",
                warmup_raw=NULL, n_samples=n_samples)
  
  #plotting_func(strain=strain, 
  #              tag="N_incidence_Cov_ignored_7,8_lognorm_non-differentiable_truncated",
  #              file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
  
  #plotting_func(strain=strain, 
  #              tag="N_incidence_Cov_ignored_7,8_lognorm_non-differentiable_non_truncated",
  #              file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
}

plotting_func_loop(strain="all", n_samples=100)

plotting_func_loop_strains <- function(strain_vector=c("alpha", "delta", "omicron", "wt"), n_samples){
  for(i in strain_vector){
    plotting_func(strain=i, 
                  tag="1000000_N_incidence_Cov_ignored_NA_lognorm_non-differentiable_truncated_skew_normal",
                  warmup_raw=NULL, n_samples=n_samples)
  }
}

plotting_func_loop_strains(n_samples=1000, strain_vector=c("wt", "alpha", "delta", "omicron"))

i <- "wt"

plotting_func_peak <- function(strain, n_samples){
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_NA_lognorm_non-differentiable_non-truncated_skew_normal",
                file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_NA_lognorm_non-differentiable_truncated_skew_normal",
                file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_8_lognorm_non-differentiable_truncated_skew_normal",
                file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
  
  plotting_func(strain=strain, 
                tag="N_peak_Cov_ignored_8_lognorm_non-differentiable_non-truncated_skew_normal",
                file_name2=NULL, warmup_raw=NULL, n_samples=n_samples)
  
}

output_comparison_table <- function(tag, strain, vl_int_grad=c(13.698,0.328), n_samples){
  output <- output_extractor(strain, tag)
  
  colnames(output$MCMC_output)[c(5,6)] <- c("ct_min_bar", "ct_min_sigma")
  n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
  warmup <- burnin_calculator(output, nroll=10000, n_it)
  symp_delay_lim <- output$x$symp_delay_lim
  
  mean_sd_store <- matrix(nrow=n_samples, ncol=2)
  
  out1 <- data.frame(params=c("a_log_bar", "a_sigma", "b_log_bar", "b_sigma", "ct_min_bar", "ct_min_sigma",
                              "l_bar", "l_sigma", "inc1", "inc2", "inc3", "test0", "test1", "test2", "test3", "test4",
                              "test5"),
                     mean=round(unname(colMeans(output$MCMC_output[warmup:n_it,1:17])),4),
                     sd=round(unname(colSds(output$MCMC_output[warmup:n_it,1:17])),4)) 
  table1 <- ggtexttable(out1, rows = NULL, theme=ttheme("blank"))
  
  transpars <- matrix(nrow=n_samples, ncol=25)
  colnames(transpars) <- c("a", "a_lower", "a_upper",
                           "b", "b_lower", "b_upper",
                           "ctmin", "ctmin_lower", "ctmin_upper",
                           "l", "llower", "lupper",
                           "ttp", "ttp_lower", "ttp_upper",
                           "ttd", "ttd_lower", "ttd_upper",
                           "tit", "tit_lower", "tit_upper",
                           "inc", "inc_lower", "inc_upper","twenty_percent")
  
  area_samples_ct <- matrix(nrow=n_samples, ncol=1000)
  twenty_percent <- vector(length=n_samples)
    
  ## SETGET STUFF IN VL SPACE
  for(m in 1:n_samples){
    print(c("m = ", m), quote=F)
    iter <- floor(seq(warmup,n_it-1,length.out=n_samples)[m])
    ct_min_samples <- rnorm(n=1000, mean=output$MCMC_output[iter,'ct_min_bar'], sd=output$MCMC_output[iter,'ct_min_sigma'])
    log_vlmax_samples <- vl_int_grad[1]-ct_min_samples*vl_int_grad[2]
    
    lower_symp <- ifelse(output$x$form=="incidence",0,-symp_delay_lim)
    
    inc_prob <- psnorm(seq(lower_symp,symp_delay_lim,0.01), mean=output$MCMC_output[iter,'inc1'], sd=output$MCMC_output[iter,'inc2'], xi=exp(output$MCMC_output[iter,'inc3']))-
      psnorm(seq(lower_symp,symp_delay_lim,0.01)-1, mean=output$MCMC_output[iter,'inc1'], sd=output$MCMC_output[iter,'inc2'], xi=exp(output$MCMC_output[iter,'inc3']))
    
    inc_samples <- rep(seq(lower_symp,symp_delay_lim, 0.01), round(inc_prob*1000))              
    
    lower_a <- -Inf
    lower_b <- -Inf
    lower_l <- 0
    
    if(output$x$vl_truncation=="truncated"){
      lower_a <- log((40-ct_min_samples)/14)
      lower_b <- log((40-ct_min_samples)/28)
    } 
    if(output$x$vl_kinetics_type=="differentiable") lower_l <- -Inf
    
    a_samples_ct <- exp(rtruncnorm(n=1000, a=lower_a, mean=output$MCMC_output[iter,'a_bar'], sd=output$MCMC_output[iter,'a_sigma']))
    b_samples_ct <- exp(rtruncnorm(n=1000, a=lower_b, mean=output$MCMC_output[iter,'b_bar'], sd=output$MCMC_output[iter,'b_sigma']))
    
    a_samples_vl <- a_samples_ct*vl_int_grad[2]
    b_samples_vl <- b_samples_ct*vl_int_grad[2]
    
    if((is.na(output$x$ignored_pars))[1]==FALSE) l_samples <- output$MCMC_output[iter,'l_bar']
    else l_samples <- rtruncnorm(n=1000, a=lower_l, mean=output$MCMC_output[iter,'l_bar'], sd=output$MCMC_output[iter,'l_sigma'])
    
    if(output$x$vl_kinetics_type=="differentiable") log_vlmax_samples <- vl_int_grad[1]-(ct_min_samples-log((a_samples_ct+b_samples_ct)/(a_samples_ct+b_samples_ct+exp(l_samples))))*vl_int_grad[2]
    
    if(output$x$vl_kinetics_type=="differentiable"){
      ttp_samples <- (1/a_samples_ct) * log(((a_samples_ct+b_samples_ct)/b_samples_ct)*exp(40-ct_min_samples)-exp(l_samples)/b_samples_ct)
      ttd_samples <- (1/b_samples_ct) * log(((a_samples_ct+b_samples_ct)/a_samples_ct)*exp(40-ct_min_samples)-exp(l_samples)/a_samples_ct) 
      tit_samples <- ttp_samples + ttd_samples 
    }
    
    if(output$x$vl_kinetics_type!="differentiable"){
      ttp_samples = (40-ct_min_samples)/a_samples_ct + 0.5*l_samples
      ttd_samples = (40-ct_min_samples)/b_samples_ct + 0.5*l_samples
      tit_samples = ttp_samples+ttd_samples
    }
    
    gt <- (40-ct_min_samples)/a_samples_ct
    ft <- l_samples
    dt <- (40-ct_min_samples)/b_samples_ct
    
    area_samples_ct[m,] <- 
      sort((1/(a_samples_vl*log(10)))*(exp(a_samples_vl*gt*log(10))-1) +
      ft*exp(log_vlmax_samples*log(10)) +
      (1/(b_samples_vl*log(10)))*(exp(b_samples_vl*dt*log(10))-1))
    
    twenty_percent[m] <- sum(area_samples_ct[m,800:1000])/sum(area_samples_ct[m,])
    #p_symp <- vector(length=10000)
    #for(i in 1:length(p_symp)) p_symp[i] <- length(which(inc_samples < ttp_samples[i]))/length(inc_samples)
    #p_symp_pre_peak <- sum(rbinom(n=length(p_symp), size=1,prob=p_symp))/length(p_symp)
    
    transpars[m,] <- c(median(a_samples_vl), quantile(a_samples_vl, probs=0.025), quantile(a_samples_vl, probs=0.975), 
                       median(b_samples_vl), quantile(b_samples_vl, probs=0.025), quantile(b_samples_vl, probs=0.975), 
                       median(log_vlmax_samples), quantile(log_vlmax_samples, probs=0.025), quantile(log_vlmax_samples, probs=0.975),
                       median(l_samples), quantile(l_samples, probs=0.025), quantile(l_samples, probs=0.975),
                       median(ttp_samples), quantile(ttp_samples, probs=0.025), quantile(ttp_samples, probs=0.975),
                       median(ttd_samples), quantile(ttd_samples, probs=0.025), quantile(ttd_samples, probs=0.975),
                       median(tit_samples), quantile(tit_samples, probs=0.025), quantile(tit_samples, probs=0.975),
                       median(inc_samples), quantile(inc_samples, probs=0.025), quantile(inc_samples, probs=0.975),twenty_percent[m]) 
    
  }
  
  if(is.na(output$x$ignored_pars)[1]==F) transpars[,c('llower','lupper')] <- NA
  
  out2 <- data.frame(param=c("vl_growth_rate", "vl_decay_rate", "peak_vl", "flattening", "time_to_peak", "time_to_decay", "total_infection_time", "incubation", "twenty_percent"),
                     median =     round(colMedians(transpars[,seq(1,25,3)]),2),
                     mean_range=c(paste0(round(colQuantiles(transpars[,seq(1,25,3)], probs=0.025),2)," - ", round(colQuantiles(transpars[,seq(1,25,3)], probs=0.975),2))),
                     lower =    c(round(colMedians(transpars[,seq(2,24,3)]),2),NA),
                     lower_range=c(paste0(round(colQuantiles(transpars[,seq(2,24,3)], probs=0.025),2)," - ", round(colQuantiles(transpars[,seq(2,24,3)], probs=0.975),2)),NA),
                     upper =     c(round(colMedians(transpars[,seq(3,24,3)]),2),NA),
                     upper_range=c(paste0(round(colQuantiles(transpars[,seq(3,24,3)], probs=0.025),2)," - ", round(colQuantiles(transpars[,seq(3,24,3)], probs=0.975),2)),NA))
  
  table2 <- ggtexttable(out2, rows = NULL, theme=ttheme("blank"))
  
  outputs <- list(table1=table1, table2=table2)
  
}


plotting_func <- function(strain, tag="", warmup_raw=NULL, n_samples=10, vl_int_grad=c(13.698,0.328)){
  output <- output_extractor(strain, tag)
  output$x$ct_range <- c(10,39)
  
  vg <- ifelse(output$x$vl_growth_type == "lognorm", 1, 0)
  vk <- ifelse(output$x$vl_kinetics_type == "non-differentiable", 1, 0)
  vt <- ifelse(output$x$vl_truncation == "truncated", 1, 0)
  
  colnames(output$MCMC_output) <- names(output$x$prior_mean)
  colnames(output$MCMC_output)[c(5,6)] <- c("ct_min_bar", "ct_min_sigma")
  
  sourceCpp("will3p_vlmax_all_tstoch.cpp")
  sourceCpp("will3p_inc.cpp")
  
  source("ct_symptoms_key_functions.R")
  
  ncores <- 8
  
  n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
  warmup <- ifelse(is.null(warmup_raw), burnin_calculator(output, nroll=10000, n_it), warmup_raw)
  print(warmup)
  
  ct_range <- output$x$ct_range
  
  symp_delay_lim <- output$x$symp_delay_lim
  if(grepl("peak", output$x$form)) symp_days <- -symp_delay_lim:symp_delay_lim
  if(grepl("incidence", output$x$form)) symp_days <- 0:symp_delay_lim
  
  prior_vec <- vector(length=n_it)
  for(i in 1:n_it){
    #print(i)
    prior_vec[i] <- prior(output$MCMC_output[i,], output$x$prior_mean, output$x$prior_sds, output$x$form, output$x$delay_dist, vt)
  } 
  likelihood_vec <- output$MCMC_posteriors[1:n_it]-prior_vec
  
  dev <- -2*likelihood_vec[warmup:n_it]
  mean_dev <- mean(dev)
  
  iter_maximised_post <- which(output$MCMC_posteriors==max(output$MCMC_posteriors, na.rm=T))[1]
  pars_maximised_post <- output$MCMC_output[iter_maximised_post,]
  L_maximised_post <- vector(length=10)
  
  for(i in 1:length(L_maximised_post)){
    print(i)
    L_maximised_post[i] <- likelihood_function3(parameters=pars_maximised_post, output$x$knots, vg, vk, vt, output$x$delay_dist,
                                                output$x$max_day, output$x$population, weekly_aggregator(output$x$data_array),
                                                output$x$test_pop, ncores=8, output$x$form, output$x$test_form, output$x$symp_delay_lim, output$x$day_agg)
  } 
  
  post_maximised_post2 <- sort(unique(c(output$MCMC_posteriors)), decreasing=T)[1:20]
  L_maximised_post_vec2 <- vector(length=20)
  for(i in 1:length(L_maximised_post_vec2)) L_maximised_post_vec2[i] <- which(c(output$MCMC_posteriors)==post_maximised_post2[i])[1]
  L_maximised_post2 <- likelihood_vec[L_maximised_post_vec2]
  
  DIC1 <- 2*mean_dev+2*mean(L_maximised_post)
  DIC2 <- 2*mean_dev+2*mean(L_maximised_post2)
  
  DIC_table <- data.frame(metric=c("DIC", "DIC2"), value=round(c(DIC1, DIC2),1))
  
  ## traceplots
  #MCMC_output_plot <- data.frame(output$MCMC_output[1:n_it,],posterior=output$MCMC_posteriors[1:n_it,]) %>% dplyr::mutate(iteration=row_number()) %>% tidyr::gather("param", "value", 1:(ncol(.)-1)) %>% mutate(param=factor(param, levels=c(names(output$x$prior_mean), "posterior"))) %>% filter(iteration %% 100 == 0) %>%
  #  mutate(form=output$x$form, kinetics=output$x$vl_kinetics_type, truncation=output$x$vl_truncation)

  #traceplot <- ggplot(MCMC_output_plot %>% filter(iteration %% 100 == 0, iteration>warmup), aes(x=iteration, y=value)) +
  #  geom_line() +
  #  facet_wrap(~param, scales="free_y") +
  #  theme_bw() +
  #  theme(legend.position="none",
  #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #  ggtitle("MCMC traceplots from post-burnin") +
  #  scale_x_continuous(breaks=c(floor(warmup/100000)*100000,1e6))

  # Sample data
  data <- as.data.frame(matrix(rnorm(300), ncol = 3))
  
  # Create a pairs plot with contour lines
  sel_out = output$MCMC_output[warmup:n_it,c("a_bar", "b_bar", "ct_min_bar", "l_bar")] 
  
  library(GGally)
  
  pairs_plot = ggpairs(as.data.frame(sel_out[seq(1,nrow(sel_out),100),])) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Pairs plots of population mean viral kinetics parameters") 
  
  ggsave(filename="test2.png", pairs_plot, width=7.5, height=24/5)
  img <- readPNG("test2.png")
  
  pairs_plot2 = ggplot() + 
    background_image(img) 
  
  
  ## PRIOR vs POSTERIOR
  prior_df <- as.data.frame(output$MCMC_output[seq(warmup,n_it-1),1:11]) %>% tidyr::gather("param", "value", 1:ncol(.)) %>%
    mutate(truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type, ignored=paste(output$x$ignored_pars, sep="", collapse=""), form=output$x$form) %>%
    mutate(prior=rnorm(n=nrow(.), mean=output$x$prior_mean[param], sd=output$x$prior_sds[param])) 
  
  prior_plot <- ggplot(prior_df) + 
    geom_density(aes(x=value, color=truncation, linetype=kinetics)) + 
    scale_color_manual(values=c("red")) +
    geom_density(aes(x=prior), color="black", linetype="dashed") +
    facet_wrap(~param, scales="free") +
    theme_bw() +
    ggtitle("Prior vs posterior") +
    theme(legend.position="none")
  
  ct_dist_p_mat <- matrix(nrow=dim(output$x$data_array)[1], ncol=n_samples)
  npos_p_mat <- matrix(nrow=dim(output$x$data_array)[2]-35, ncol=n_samples)
  ntest_p_mat <- matrix(nrow=dim(output$x$data_array)[3], ncol=n_samples)
  
  onset_delays <- 0:(dim(output$x$data_array)[3]-1)
  ct_from_onset_mat <- matrix(NA, nrow=length(onset_delays), ncol=n_samples)
  ct_from_onset_median_mat <- matrix(NA, nrow=length(onset_delays), ncol=n_samples)
  
  time_vec <- vector(length=n_samples)
  ct_range <- output$x$ct_range
  max_inf <- output$x$max_inf
  
  inf_lims <- floor(max_inf/2)
  
  if(output$x$form=="peak") symp_days <- -output$x$symp_delay_lim:output$x$symp_delay_lim
  if(output$x$form=="incidence") symp_days <- 0:output$x$symp_delay_lim
  
  inc_mat <- matrix(nrow=length(symp_days), ncol=n_samples)
  cum_inc_mat <- matrix(nrow=length(symp_days), ncol=n_samples)
  
  vl_posterior_array <- array(dim=c(length(0:(max_inf-1)), n_samples, 5))
  pseeds <- unlist(generateSeedVectors(ncores,1)) 
  
  for(k in 1:n_samples){
    p <- proc.time()
    print(k)
    iter <- floor(seq(warmup,n_it-1,length.out=n_samples)[k])
    parameters <- output$MCMC_output[iter,]
    p_array_raw <- p_array_func(parameters, knots=output$x$knots, vg=ifelse(output$x$vl_growth_type=="lognorm", 1, 0), vk=ifelse(output$x$vl_kinetics_type=="non-differentiable", 1, 0), vt=ifelse(output$x$vl_truncation=="truncated", 1, 0), delay_dist=output$x$delay_dist, max_day=output$x$max_day, population=6.7e7, data_array = output$x$data_array, test_pop=1e7, ncores=ncores, form=output$x$form, test_form=output$x$test_form, symp_delay_lim=output$x$symp_delay_lim, stoch=0.5)
    p_array_raw[is.na(p_array_raw)] <- 0
    
    p_array <-  array(rnbinom(n=length(p_array_raw), size=1/output$MCMC_output[iter,'rdisp'], mu=p_array_raw), dim=dim(p_array_raw))[,-c(1:35),]
    dimnames(p_array) = dimnames(output$x$data_array[,-c(1:35),])
    
    pseeds <- unlist(generateSeedVectors(ncores,1)) 
    
    #vl_posterior[,k] <- p_ct_tau_mat %*% c(10:35)
    time_vec[k] <- (proc.time()-p)["elapsed"]
    
    # 1D fits
    ct_dist_p_mat[,k] <- apply(p_array,1, sum)
    npos_p_mat[,k] <- apply(p_array,2, sum)
    ntest_p_mat[,k] <- apply(p_array,3,sum)
    
    # 2d fit (ct dist by symp day)
    p_array_cdays <- apply(p_array,c(1,3),sum)
    
    if(k==1) df_p_raw <-                 p_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1))) %>% mutate(ct=10:(dim(output$x$data_array)[1]+9)) %>% tidyr::gather(symp_day, value_p, 1:(ncol(.)-1))%>% magrittr::set_colnames(c("ct", "symp_day", paste0("value",k)))
    else df_p_raw <- left_join(df_p_raw, p_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1))) %>% mutate(ct=10:(dim(output$x$data_array)[1]+9)) %>% tidyr::gather(symp_day, value_p, 1:(ncol(.)-1)) %>% magrittr::set_colnames(c("ct", "symp_day", paste0("value",k))), by=c("ct", "symp_day"))
    
    # 2d fit (mean ct by symp day)
    ct_symp_day_mat <- apply(p_array, c(3,1), sum)
    ct_from_onset_mat[,k] <- (ct_symp_day_mat %*% (ct_range[1]:(ct_range[2]))) / rowSums(ct_symp_day_mat) #mean
    
    if(output$x$form=="peak"){
      p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(output$x$max_inf-1, output$x$ct_range[2]+1, output$x$ct_range[1], output$x$test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt)
      dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, c(ct_range[1]:(ct_range[2]+1)))
      
      peak_to_symp <- setNames(psnorm(symp_days, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3']))-psnorm(symp_days-1, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3'])),symp_days)
      norm_peak_to_symp <- setNames(peak_to_symp/sum(peak_to_symp),symp_days)
      cum_norm_inf_to_symp <- setNames(cumsum(norm_peak_to_symp),symp_days)
      inc_mat[,k] <- norm_peak_to_symp
      cum_inc_mat[,k] <- cum_norm_inf_to_symp
    } 
    if(output$x$form=="incidence"){
      p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(output$x$max_inf-1, output$x$ct_range[2]+1, output$x$ct_range[1], output$x$test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt)
      dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, c(ct_range[1]:(ct_range[2]+1)))
      
      #p_ct_tau_mat2 <- ct_func2_cpp_inc(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(output$x$max_inf-1, output$x$ct_range[2]+1, output$x$ct_range[1], output$x$test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt, stoch=0.0001)
      #dimnames(p_ct_tau_mat2) <- list(0:(max_inf-1), c(ct_range[1]:(ct_range[2]+1)))
      
      inf_to_symp <- setNames(psnorm(symp_days, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3']))-psnorm(symp_days-1, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3'])),symp_days)
      norm_inf_to_symp <- setNames(inf_to_symp/sum(inf_to_symp),symp_days)
      cum_norm_inf_to_symp <- setNames(cumsum(norm_inf_to_symp),symp_days)
      inc_mat[,k] <- norm_inf_to_symp
      cum_inc_mat[,k] <- cum_norm_inf_to_symp
    }
    
    for(i in 1:nrow(p_ct_tau_mat)){
      vl_list <- vl_int_grad[1] - vl_int_grad[2]*rep(as.numeric(names(p_ct_tau_mat[i,])),round(p_ct_tau_mat[i,]*100000))
      vl_posterior_array[i,k,] <- quantile(vl_list, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    } 
    
  }
  
  ct_symp_day_data_mat <- apply(output$x$data_array, c(3,1), sum)
  ct_from_onset_data <- (ct_symp_day_data_mat %*% (ct_range[1]:(ct_range[2]))) / apply(output$x$data_array, 3, sum) 
  
  ct_from_onset_posterior <- data.frame(day=onset_delays,
                                        x=rowQuantiles(ct_from_onset_mat, probs=c(0.025, 0.5, 0.975)),
                                        data=ct_from_onset_data) %>%
    magrittr::set_colnames(c("day", "lower", "median", "upper", "data"))
  
  plot_ct_day <- ggplot(ct_from_onset_posterior, aes(x=day)) +
    geom_point(aes(y=data)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
    theme_bw() +
    ggtitle("Mean ct value by days after symptom onset")
  
  print(plot_ct_day)
  
  ct_dist_data <- apply(output$x$data_array, 1, sum)
  npos_data <- apply(output$x$data_array, 2, sum)[36:dim(output$x$data_array)[2]]
  ntest_data <- apply(output$x$data_array, 3, sum)
  
  ct_dist_plot <- data.frame(ct=ct_range[1]:ct_range[2],
                             rowQuantiles(ct_dist_p_mat, probs=c(0.025, 0.5, 0.975)),
                             data=ct_dist_data,
                             multinomialCI(ct_dist_data, alpha=0.95)*sum(ct_dist_data)) %>%
    magrittr::set_colnames(c("ct", "lower", "count", "upper", "data", "data_lower", "data_upper"))
  
  npos_dist_plot <- data.frame(day=36:dim(output$x$data_array)[2],
                               rowQuantiles(npos_p_mat, probs=c(0.025, 0.5, 0.975)),
                               data=npos_data) %>%
    magrittr::set_colnames(c("day", "lower", "positives", "upper", "data"))
  
  ntest_dist_plot <- data.frame(test_day=0:(dim(output$x$data_array)[3]-1),
                                rowQuantiles(ntest_p_mat, probs=c(0.025, 0.5, 0.975)),
                                data=ntest_data) %>%
    magrittr::set_colnames(c("test_day", "lower", "tests", "upper", "data"))
  
  df <- tibble(x = max(ct_dist_plot$ct), y = max(ct_dist_plot$count), tb = list(DIC_table))
  
  # using defaults
  plot_ct <- ggplot(ct_dist_plot, aes(x=ct)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=count)) +
    geom_pointrange(aes(y=data, ymin=data_lower, ymax=data_upper)) +
    #geom_point(aes(y=data)) +
    ggtitle("Count of ct values in data and model fit") +
    geom_table(data=df,aes(x=x, y=y, label=tb)) +
    theme_bw() 
  
  plot_npos <- ggplot(npos_dist_plot %>% tibble::rownames_to_column(var="date") %>% mutate(date=as.Date(date)), aes(x=date)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=positives)) +
    geom_point(aes(y=data)) +
    theme_bw() +
    ggtitle("Count of positives in data and model fit") +
    scale_x_date(date_breaks = ifelse(output$x$strain=="all", "3 months", "1 month"), date_labels =  "%b-%Y") 
  
  # prev_plot <- ggplot(prev_df %>% filter(upper < 0.08), aes(x=date)) +
  #   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
  #   geom_line(aes(y=median)) +
  #   geom_point(aes(y=data)) +
  #   geom_vline(data=data.frame(knots=knot_days), aes(xintercept=knots), color="red", linetype="dashed") +
  #   theme_bw() +
  #   ggtitle("Infection prevelence by day in data and model") +
  #   scale_x_date(date_breaks = date_breaks, date_labels =  "%b-%Y") 
  # 
  plot_ntest <- ggplot(ntest_dist_plot, aes(x=test_day)) +
    #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=tests)) +
    geom_point(aes(y=data)) +
    theme_bw() +
    ggtitle("Count of positives by day of test since onset")
  
  data_array2 <- output$x$data_array[,-c(1:35),]
  
  d_array_cdays <- matrix(nrow=dim(data_array2)[1], ncol=dim(data_array2)[3])
  
  for(i in 1:nrow(p_array_cdays)) for(j in 1:ncol(p_array_cdays)){
    d_array_cdays[i,j] <- sum(output$x$data_array[i,,j])
  }
  
  df_d <- (d = d_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1)))%>% mutate(ct=output$x$ct_range[1]:output$x$ct_range[2]) %>% tidyr::gather(symp_day, value_d, 1:(ncol(.)-1)))
  df_p <- df_p_raw[,c(1,2)] %>% mutate(lower=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.025)),
                                       median=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.5)),
                                       upper=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.975)),)
  
  ct_plot <- left_join(df_p, df_d, by=c("symp_day", "ct")) %>% mutate(symp_day=factor(symp_day, levels=onset_delays)) %>% mutate(symp_day=paste0("symp_day=",symp_day))
  
  plot_symp_ct <- ggplot(ct_plot, aes(x=ct)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_point(aes(y=value_d)) +
    facet_wrap(~symp_day, scales="free_y") +
    theme_bw() +
    labs(y="count") +
    ggtitle("Count of ct values by day of test since symptoms onset in data (points) and model fit (line)")
  
  plot3 <- ggplot(ct_plot %>% filter(ct<40) %>% group_by(symp_day) %>% mutate(symp_day=substr(symp_day, 10, nchar(symp_day))) %>% mutate(median=median/sum(median)), aes(x=ct)) +
    geom_line(aes(y=median, color=symp_day)) +
    #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    #geom_point(aes(y=value_d)) +
    #facet_wrap(~symp_day, ncol=1, scales="free_y") +
    theme_bw() +
    labs(y="count") +
    ggtitle("Count of ct values by day of test since symptoms onset in model fit")
  
  ## Overall counts
  #if(output$x$form=="peak") days <- -inf_lims:inf_lims
  #if(output$x$form=="incidence") days <- 0:(max_inf-1)
  
  days <- -inf_lims:inf_lims
  
  vl_posterior_df <- data.frame(day=days,
                                     apply(vl_posterior_array, 3, rowMedians)) %>%
    magrittr::set_colnames(c("day", "lower", "lower_quart", "median", "upper_quart", "upper")) %>%
    mutate(form=output$x$form, truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type)
  
  vl_posterior_plot <- ggplot(vl_posterior_df %>% filter(is.na(median)==F), aes(x=day)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_ribbon(aes(ymin=lower_quart, ymax=upper_quart), alpha=0.2) +
    labs(y=expression(log[10](vl)), x="days from peak VL") +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle("Viral load value by day of infection relative to peak")

  cum_inc_posterior <- data.frame(day=symp_days,
                                       x=rowQuantiles(cum_inc_mat, probs=c(0.025, 0.5, 0.975)),
                                       x2=rowQuantiles(inc_mat, probs=c(0.025, 0.5, 0.975))) %>%
    magrittr::set_colnames(c("day", "lower", "median", "upper", "lower2", "median2", "upper2")) %>%
    mutate(form=output$x$form, truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type)
  
  cum_inc_plot <- ggplot(cum_inc_posterior, aes(x=day, fill=truncation, linetype=kinetics)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower2/max(median2), ymax=upper2/max(median2)), alpha=0.2, fill="green") +
    geom_line(aes(y=median2/max(median2)), color="green") +
    theme_bw() +
    labs(y="proportion") +
    theme(legend.position="none") +
    ggtitle(ifelse(output$x$form=="incidence", "Model fit of incubation distribution with 95% CI", "Cumulative probability of symptom onset from peak"))

  tables <- output_comparison_table(tag, strain, vl_int_grad, n_samples)
  
  top <- ggarrange(plot_ct, plot_ntest, plot_npos, ncol=3, labels=c("A", "B", "C"))
  second <- ggarrange(plot_ct_day, cum_inc_plot, vl_posterior_plot, ncol=3, labels=c("D", "E", "F"))
  third <- ggarrange(plot_symp_ct, plot3, ncol=2, labels=c("G", "H"))
  fourth <- ggarrange(pairs_plot2, prior_plot, ncol=2, labels=c("I", "J"))
  fifth <- ggarrange(tables$table1, tables$table2, ncol=2, labels=c("K", "L"))
  
  output_plot1 <- annotate_figure(ggarrange(top, second, third, fourth, nrow=4, heights=c(1, 1, 1, 1)), 
                                  top = text_grob(paste0(strain, ": mean posterior = ",round(mean(output$MCMC_posteriors[warmup:(n_it-2)],0)), "; n_it = ", n_it, "; n_samples = ", n_samples),color="red", size=14, face="bold"), 
                                  bottom=text_grob(paste("file = ", tag),color = "black", size = 14, hjust=1))
  
  output_plot <- annotate_figure(ggarrange(top, second, third, fourth, fifth, nrow=5, heights=c(1, 1, 1, 1, 1)), 
                                 top = text_grob(paste0(strain, ": mean posterior = ",round(mean(output$MCMC_posteriors[warmup:(n_it-2)],0)), "; n_it = ", n_it, "; n_samples = ", n_samples),color="red", size=14, face="bold"), 
                                 bottom=text_grob(paste("file = ", tag),color = "black", size = 14, hjust=1))
  
  
  ggsave(filename=paste0(output$x$strain,"/",output$x$strain,"_",output$x$gene,"_",output$x$form,"_",output$x$vl_kinetics,"_",output$x$vl_truncation,"_ignored_",paste(output$x$ignored_pars,sep="",collapse=","),"_",output$x$ct_range[2],"_final.png"), plot = output_plot1, device = "png", width=15, height=24)
  ggsave(filename=paste0(output$x$strain,"/",output$x$strain,"_",output$x$gene,"_",output$x$form,"_",output$x$vl_kinetics,"_",output$x$vl_truncation,"_ignored_",paste(output$x$ignored_pars,sep="",collapse=","),"_",output$x$ct_range[2],"_final.pdf"), output_plot, width=16, height = 25)
  
  #return(time_vec)
}


output_list_generator <- function(tag_list, strain){
  output_list <- list()
  for(i in tag_list) output_list[[i]] <- output_extractor(strain, i)
  return(output_list)
}

traceplot_func <- function(tag_list, strain, warmup){
  output_list <- output_list_generator(tag_list, strain)
  
  MCMC_output_plot <- list()
  for(j in 1:length(output_list)){
    output <- output_list[[j]]
    colnames(output$MCMC_output)[c(5,6)] <- c("ct_min_bar", "ct_min_sigma")
    n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
    MCMC_output_plot[[j]] <- data.frame(output$MCMC_output[1:n_it,],posterior=output$MCMC_posteriors[1:n_it,]) %>% dplyr::mutate(iteration=row_number()) %>% tidyr::gather("param", "value", 1:(ncol(.)-1)) %>% mutate(param=factor(param, levels=c(names(output$x$prior_mean), "posterior"))) %>% filter(iteration %% 100 == 0) %>%
      mutate(form=output$x$form, kinetics=output$x$vl_kinetics_type, truncation=output$x$vl_truncation, 
             ignored=output$x$ignored_pars)
  }
  
  joint_traceplot <- do.call("rbind", MCMC_output_plot)
  
  traceplot <- ggplot(joint_traceplot %>% filter(iteration %% 100 == 0), aes(x=iteration, y=value)) +
    geom_line() +
    facet_wrap(~param, scales="free_y") +
    theme_bw() +
    theme(legend.position=ifelse(length(tag_list)==1, "none")) +
    ggtitle("MCMC traceplots from post-burnin")
  
  return(traceplot)
}


traceplot_func(tag_list=c("N_incidence_Cov_ignored_8_lognorm_non-differentiable_non-truncated_skew_normal_"),
               strain="all")

tag_list <- c("N_peak_Cov_ignored_NA_lognorm_non-differentiable_non-truncated")
              

prior_plot_func <- function(tag_list, strain){
  prior_df <- list()
  output_list <- output_list_generator(tag_list, strain)
  
  for(i in 1:length(tag_list)){
    output <- output_list[[i]]  
    colnames(output$MCMC_output)[c(5,6)] <- c("ct_min_bar", "ct_min_sigma")
    n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
    warmup <- burnin_calculator(output, nroll=10000, n_it)
    
    prior_df[[i]] <- as.data.frame(output$MCMC_output[seq(warmup,n_it-1),1:11]) %>% tidyr::gather("param", "value", 1:ncol(.)) %>%
      mutate(truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type, ignored=output$x$ignored_pars, form=output$x$form) %>%
      mutate(prior=rnorm(n=nrow(.), mean=output$x$prior_mean[param], sd=output$x$prior_sds[param])) 
  }
  
  prior_df_plotting <- do.call(rbind, prior_df)
  
  prior_plot <- ggplot(prior_df_plotting) + 
    geom_density(aes(x=value, color=truncation, linetype=kinetics)) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9")) +
    geom_density(aes(x=prior), color="black", linetype="dashed") +
    facet_wrap(~param, scales="free") +
    theme_bw() +
    ggtitle("Prior vs posterior") +
    theme(legend.position=ifelse(length(tag_list)==1, "none"))
  
  return(prior_plot)
}

inc_period_plot_func <- function(tag_list, n_samples){
  output_list <- output_list_generator(tag_list, strain)
  cum_inc_posterior <- list()
  
  for(i in 1:length(output_list)){
    output <- output_list[[i]]
    n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
    warmup <- burnin_calculator(output, nroll=10000, n_it)
    
    if(output$x$form=="peak") symp_days <- -output$x$symp_delay_lim:output$x$symp_delay_lim
    if(output$x$form=="incidence") symp_days <- 0:output$x$symp_delay_lim
    
    inc_mat <- matrix(nrow=length(symp_days), ncol=n_samples)
    cum_inc_mat <- matrix(nrow=length(symp_days), ncol=n_samples)
    
    for(k in 1:n_samples){
      iter <- floor(seq(warmup,n_it-1,length.out=n_samples)[k])
      parameters <- output$MCMC_output[iter,]
      
      if(output$x$form=="peak"){
        peak_to_symp <- setNames(psnorm(symp_days, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3']))-psnorm(symp_days-1, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3'])),symp_days)
        norm_peak_to_symp <- setNames(cumsum(peak_to_symp/sum(peak_to_symp)),symp_days)
        cum_norm_inf_to_symp <- setNames(cumsum(norm_inf_to_symp),symp_days)
        inc_mat[,k] <- norm_peak_to_symp
        cum_inc_mat[,k] <- cum_norm_inf_to_symp
      }
      if(output$x$form=="incidence"){
        inf_to_symp <- setNames(psnorm(symp_days, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3']))-psnorm(symp_days-1, mean=parameters['inc1'], sd=parameters['inc2'], xi=exp(parameters['inc3'])),symp_days)
        norm_inf_to_symp <- setNames(inf_to_symp/sum(inf_to_symp),symp_days)
        cum_norm_inf_to_symp <- setNames(cumsum(norm_inf_to_symp),symp_days)
        inc_mat[,k] <- norm_inf_to_symp
        cum_inc_mat[,k] <- cum_norm_inf_to_symp
      }
    }
    
    cum_inc_posterior[[i]] <- data.frame(day=symp_days,
                                    x=rowQuantiles(cum_inc_mat, probs=c(0.025, 0.5, 0.975)),
                                    x2=rowQuantiles(inc_mat, probs=c(0.025, 0.5, 0.975))) %>%
      magrittr::set_colnames(c("day", "lower", "median", "upper", "lower2", "median2", "upper2")) %>%
      mutate(form=output$x$form, truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type,
             ignored=output$x$ignored_pars)
    
  }
  
  cum_inc_posterior_df <- do.call(rbind, cum_inc_posterior)
  
  cum_inc_plot <- ggplot(cum_inc_posterior_df, aes(x=day, fill=truncation, linetype=kinetics)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower2/max(median2), ymax=upper2/max(median2)), alpha=0.4) +
    geom_line(aes(y=median2/max(median2))) +
    theme_bw() +
    facet_wrap(~form) +
    labs(y="proportion") +
    theme(legend.position=ifelse(length(tag_list)==1, "none")) +
    ggtitle(ifelse(output$x$form=="incidence", "Model fit of incubation distribution with 95% CI", "Cumulative probability of symptom onset from peak"))
  
  return(cum_inc_plot)
  
}

viral_kinetics_plot_func <- function(tag_list, strain, n_samples){
  output_list <- output_list_generator(tag_list, strain)
  vl_posterior_df <- list()
  
  for(m in 1:length(output_list)) {
    output <- output_list[[m]]
    n_it <- pmin(which(is.na(output$MCMC_posteriors)==T)[1]-2, output$x$n_iterations, na.rm=T)
    warmup <- burnin_calculator(output, nroll=10000, n_it)
    ct_range <- output$x$ct_range
    
    max_inf <- output$x$max_inf
    lims <- floor(output$x$max_inf/2)
    inf_lims <- floor(max_inf/2)
    
    vg <- ifelse(output$x$vl_growth_type == "lognorm", 1, 0)
    vk <- ifelse(output$x$vl_kinetics_type == "non-differentiable", 1, 0)
    vt <- ifelse(output$x$vl_truncation == "truncated", 1, 0)
    
    vl_posterior_array <- array(dim=c(length(0:(max_inf-1)), n_samples, 5))
    
    pseeds <- unlist(generateSeedVectors(ncores,1)) 
    
    for(k in 1:n_samples){
      print(k)
      iter <- floor(seq(warmup,n_it-1,length.out=n_samples)[k])
      parameters <- output$MCMC_output[iter,]
      
      if(output$x$form=="peak"){
        p_ct_tau_mat <- ct_func2_cpp_all(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(output$x$max_inf-1, output$x$ct_range[2]+1, output$x$ct_range[1], output$x$test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt)
        dimnames(p_ct_tau_mat) <- list(-inf_lims:inf_lims, c(ct_range[1]:(ct_range[2]+1)))
      } 
      if(output$x$form=="incidence"){
        p_ct_tau_mat <- ct_func2_cpp_inc(a=parameters[c(1,2)], b=parameters[c(3,4)], c=parameters[c(5,6)], l=parameters[c(7,8)], p=c(output$x$max_inf-1, output$x$ct_range[2]+1, output$x$ct_range[1], output$x$test_pop), nthreads=ncores, pseeds=pseeds, vg, vk, vt, stoch=0.0001)
        dimnames(p_ct_tau_mat) <- list(0:(max_inf-1), c(ct_range[1]:(ct_range[2]+1)))
      }
      
      for(i in 1:nrow(p_ct_tau_mat)){
        ct_list <- rep(as.numeric(names(p_ct_tau_mat[i,])),round(p_ct_tau_mat[i,]*100000))
        vl_posterior_array[i,k,] <- quantile(ct_list, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
      } 
    }
    
    if(output$x$form=="peak") days <- -inf_lims:inf_lims
    if(output$x$form=="incidence") days <- 0:(max_inf-1)
    
    vl_posterior_df[[m]] <- data.frame(day=days,
                                  apply(vl_posterior_array, 3, rowMedians)) %>%
      magrittr::set_colnames(c("day", "lower", "lower_quart", "median", "upper_quart", "upper")) %>%
      mutate(form=output$x$form, truncation=output$x$vl_truncation, kinetics=output$x$vl_kinetics_type,
             ignored=output$x$ignored_pars)
    
    
  }
  
  vl_posterior_plot_df <- do.call(rbind, vl_posterior_df)
  
  vl_posterior_plot <- ggplot(vl_posterior_plot_df %>% filter(is.na(median)==F), aes(x=day)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_ribbon(aes(ymin=lower_quart, ymax=upper_quart), alpha=0.2) +
    labs(x="day of infection") +
    theme_bw() +
    theme(legend.position=ifelse(length(tag_list)==1, "none")) +
    ggtitle("Ct value by day of infection")
  
  return(vl_posterior_plot)

}

data_fits_plot_func <- function(output, warmup, n_it, n_samples){
  ct_dist_p_mat <- matrix(nrow=dim(output$x$data_array)[1], ncol=n_samples)
  npos_p_mat <- matrix(nrow=dim(output$x$data_array)[2]-35, ncol=n_samples)
  ntest_p_mat <- matrix(nrow=dim(output$x$data_array)[3], ncol=n_samples)
  
  onset_delays <- 0:(dim(output$x$data_array)[3]-1)
  ct_from_onset_mat <- matrix(NA, nrow=length(onset_delays), ncol=n_samples)
  ct_from_onset_median_mat <- matrix(NA, nrow=length(onset_delays), ncol=n_samples)
  
  time_vec <- vector(length=n_samples)
  ct_range <- output$x$ct_range
  
  for(k in 1:n_samples){
    p <- proc.time()
    print(k)
    iter <- floor(seq(warmup,n_it-1,length.out=n_samples)[k])
    parameters <- output$MCMC_output[iter,]
    p_array_raw <- p_array_func(parameters, knots=output$x$knots, vg=ifelse(output$x$vl_growth_type=="lognorm", 1, 0), vk=ifelse(output$x$vl_kinetics_type=="non-differentiable", 1, 0), vt=ifelse(output$x$vl_truncation=="truncated", 1, 0), delay_dist=output$x$delay_dist, max_day=output$x$max_day, population=6.7e7, data_array = output$x$data_array, test_pop=1e7, ncores=ncores, form=output$x$form, test_form=output$x$test_form, symp_delay_lim=output$x$symp_delay_lim, stoch=0.5)
    p_array_raw[is.na(p_array_raw)] <- 0
    
    p_array <-  array(rnbinom(n=length(p_array_raw), size=1/output$MCMC_output[iter,'rdisp'], mu=p_array_raw), dim=dim(p_array_raw))[,-c(1:35),]
    dimnames(p_array) = dimnames(output$x$data_array[,-c(1:35),])
    
    pseeds <- unlist(generateSeedVectors(ncores,1)) 
    
    #vl_posterior[,k] <- p_ct_tau_mat %*% c(10:35)
    time_vec[k] <- (proc.time()-p)["elapsed"]
    
    # 1D fits
    ct_dist_p_mat[,k] <- apply(p_array,1, sum)
    npos_p_mat[,k] <- apply(p_array,2, sum)
    ntest_p_mat[,k] <- apply(p_array,3,sum)
    
    # 2d fit (ct dist by symp day)
    p_array_cdays <- apply(p_array,c(1,3),sum)
    
    if(k==1) df_p_raw <-                 p_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1))) %>% mutate(ct=10:(dim(output$x$data_array)[1]+9)) %>% tidyr::gather(symp_day, value_p, 1:(ncol(.)-1))%>% magrittr::set_colnames(c("ct", "symp_day", paste0("value",k)))
    else df_p_raw <- left_join(df_p_raw, p_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1))) %>% mutate(ct=10:(dim(output$x$data_array)[1]+9)) %>% tidyr::gather(symp_day, value_p, 1:(ncol(.)-1)) %>% magrittr::set_colnames(c("ct", "symp_day", paste0("value",k))), by=c("ct", "symp_day"))
    
    # 2d fit (mean ct by symp day)
    ct_symp_day_mat <- apply(p_array, c(3,1), sum)
    ct_from_onset_mat[,k] <- (ct_symp_day_mat %*% (ct_range[1]:(ct_range[2]))) / rowSums(ct_symp_day_mat) #mean
  }
  
  ct_symp_day_data_mat <- apply(output$x$data_array, c(3,1), sum)
  ct_from_onset_data <- (ct_symp_day_data_mat %*% (ct_range[1]:(ct_range[2]))) / apply(output$x$data_array, 3, sum) 
  
  ct_from_onset_posterior <- data.frame(day=onset_delays,
                                        x=rowQuantiles(ct_from_onset_mat, probs=c(0.025, 0.5, 0.975)),
                                        data=ct_from_onset_data) %>%
    magrittr::set_colnames(c("day", "lower", "median", "upper", "data"))
  
  plot_ct_day <- ggplot(ct_from_onset_posterior, aes(x=day)) +
    geom_point(aes(y=data)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
    theme_bw() +
    ggtitle("Mean ct value by days after symptom onset")
  
  print(plot_ct_day)
  
  ct_dist_data <- apply(output$x$data_array, 1, sum)
  npos_data <- apply(output$x$data_array, 2, sum)[36:dim(output$x$data_array)[2]]
  ntest_data <- apply(output$x$data_array, 3, sum)
  
  ct_dist_plot <- data.frame(ct=ct_range[1]:ct_range[2],
                             rowQuantiles(ct_dist_p_mat, probs=c(0.025, 0.5, 0.975)),
                             data=ct_dist_data,
                             multinomialCI(ct_dist_data, alpha=0.95)*sum(ct_dist_data)) %>%
    magrittr::set_colnames(c("ct", "lower", "count", "upper", "data", "data_lower", "data_upper"))
  
  npos_dist_plot <- data.frame(day=36:dim(output$x$data_array)[2],
                               rowQuantiles(npos_p_mat, probs=c(0.025, 0.5, 0.975)),
                               data=npos_data) %>%
    magrittr::set_colnames(c("day", "lower", "positives", "upper", "data"))
  
  ntest_dist_plot <- data.frame(test_day=0:(dim(output$x$data_array)[3]-1),
                                rowQuantiles(ntest_p_mat, probs=c(0.025, 0.5, 0.975)),
                                data=ntest_data) %>%
    magrittr::set_colnames(c("test_day", "lower", "tests", "upper", "data"))
  
  df <- tibble(x = max(ct_dist_plot$ct), y = max(ct_dist_plot$count), tb = list(DIC_table))
  
  # using defaults
  plot_ct <- ggplot(ct_dist_plot, aes(x=ct)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=count)) +
    geom_pointrange(aes(y=data, ymin=data_lower, ymax=data_upper)) +
    #geom_point(aes(y=data)) +
    ggtitle("Count of ct values in data and model fit") +
    geom_table(data=df,aes(x=x, y=y, label=tb)) +
    theme_bw() 
  
  plot_npos <- ggplot(npos_dist_plot, aes(x=day)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=positives)) +
    geom_point(aes(y=data)) +
    theme_bw() +
    ggtitle("Count of positives in data and model fit")
  
  plot_ntest <- ggplot(ntest_dist_plot, aes(x=test_day)) +
    #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_line(aes(y=tests)) +
    geom_point(aes(y=data)) +
    theme_bw() +
    ggtitle("Count of positives by day of test since onset")
  
  data_array2 <- output$x$data_array[,-c(1:35),]
  
  d_array_cdays <- matrix(nrow=dim(data_array2)[1], ncol=dim(data_array2)[3])
  
  for(i in 1:nrow(p_array_cdays)) for(j in 1:ncol(p_array_cdays)){
    d_array_cdays[i,j] <- sum(output$x$data_array[i,,j])
  }
  
  df_d <- (d = d_array_cdays %>% as.data.frame() %>% magrittr::set_colnames(c(0:(ncol(.)-1)))%>% mutate(ct=output$x$ct_range[1]:output$x$ct_range[2]) %>% tidyr::gather(symp_day, value_d, 1:(ncol(.)-1)))
  df_p <- df_p_raw[,c(1,2)] %>% mutate(lower=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.025)),
                                       median=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.5)),
                                       upper=rowQuantiles(as.matrix(df_p_raw[c(3:ncol(df_p_raw))]), probs=c(0.975)),)
  
  ct_plot <- left_join(df_p, df_d, by=c("symp_day", "ct")) %>% mutate(symp_day=factor(symp_day, levels=onset_delays)) %>% mutate(symp_day=paste0("symp_day=",symp_day))
  
  plot_symp_ct <- ggplot(ct_plot, aes(x=ct)) +
    geom_line(aes(y=median)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    geom_point(aes(y=value_d)) +
    facet_wrap(~symp_day, scales="free_y") +
    theme_bw() +
    labs(y="count") +
    ggtitle("Count of ct values by day of test since symptoms onset in data (points) and model fit (line)")
  
  plot3 <- ggplot(ct_plot %>% filter(ct<40) %>% group_by(symp_day) %>% mutate(symp_day=substr(symp_day, 10, nchar(symp_day))) %>% mutate(median=median/sum(median)), aes(x=ct)) +
    geom_line(aes(y=median, color=symp_day)) +
    #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    #geom_point(aes(y=value_d)) +
    #facet_wrap(~symp_day, ncol=1, scales="free_y") +
    theme_bw() +
    labs(y="count") +
    ggtitle("Count of ct values by day of test since symptoms onset in model fit")
  
  return(list(plot_npos=plot_npos, plot_ntest=plot_ntest, plot_ct=plot_ct, plot_symp_ct=plot_symp_ct, plot_ct_day=plot_ct_day, plot3=plot3))
}

DIC_table_func <- function(output, warmup, n_it, vg, vk, vt){
  prior_vec <- vector(length=n_it)
  for(i in 1:n_it){
    #print(i)
    prior_vec[i] <- prior(output$MCMC_output[i,], output$x$prior_mean, output$x$prior_sds, output$x$form, output$x$delay_dist, vt)
  } 
  likelihood_vec <- output$MCMC_posteriors[1:n_it]-prior_vec
  
  dev <- -2*likelihood_vec[warmup:n_it]
  mean_dev <- mean(dev)
  
  iter_maximised_post <- which(output$MCMC_posteriors==max(output$MCMC_posteriors, na.rm=T))[1]
  pars_maximised_post <- output$MCMC_output[iter_maximised_post,]
  L_maximised_post <- vector(length=10)
  
  for(i in 1:length(L_maximised_post)){
    print(i)
    L_maximised_post[i] <- likelihood_function3(parameters=pars_maximised_post, output$x$knots, vg, vk, vt, output$x$delay_dist,
                                                output$x$max_day, output$x$population, weekly_aggregator(output$x$data_array),
                                                output$x$test_pop, ncores=8, output$x$form, output$x$test_form, output$x$symp_delay_lim, output$x$day_agg)
  } 
  
  post_maximised_post2 <- sort(unique(c(output$MCMC_posteriors)), decreasing=T)[1:20]
  L_maximised_post_vec2 <- vector(length=20)
  for(i in 1:length(L_maximised_post_vec2)) L_maximised_post_vec2[i] <- which(c(output$MCMC_posteriors)==post_maximised_post2[i])[1]
  L_maximised_post2 <- likelihood_vec[L_maximised_post_vec2]
  
  AIC_first_term <- 2*(length(output$x$parameters)-length(output$x$ignored_pars)-1)
  AIC1 <- AIC_first_term-2*mean(L_maximised_post)
  AIC2 <- AIC_first_term-2*mean(L_maximised_post2)
  
  DIC1 <- mean_dev+2*mean(L_maximised_post)
  DIC2 <- mean_dev+2*mean(L_maximised_post2)
  
  DIC_table <- data.frame(metric=c("AIC1", "AIC2", "DIC", "DIC2"), value=round(c(AIC1, AIC2, DIC1, DIC2),1))
  
}



sum_acceptances <- function(tag){
  output <- output_extractor("all", tag)
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  
  plotting <- data.frame(posterior=output$MCMC_posteriors[1:n_it], iteration=1:n_it)
  plot1 <- ggplot(data=plotting, aes(x=iteration, y=posterior)) + geom_line() + theme_bw()
  plot2 <- ggplot(data=plotting %>% filter(iteration>n_it/2), aes(x=iteration, y=posterior)) + geom_line() + theme_bw()
  
  print(ggarrange(plot1, plot2))
  
  print(c(pre_cov=paste0(round(sum(output$Acceptances[1:min(n_it,150000)], na.rm=T)*100/(length(which(is.na(output$Acceptances[1:min(150000,n_it)])==F))),2),"%"),
          post_cov=paste0(round(sum(output$Acceptances[150000:n_it])*100/(n_it-150000),2),"%")), quote=F)
  

}

































x_all_1 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_fixes_all_ab_chain1", test_form="empirical", 
                           ignored_pars=c(1,2,3,4), spline_burnin=0, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=110000)
x_all_1$parameters[c(1:4)] <- x_all_1$prior_mean[c(1:4)]
x_all_1$parameters[c(5:61)] <- readRDS("all/2023-03-25_all_536000_peak_Cov_empirical_ignored_NA_more_knots_chain1.rds")$MCMC_output[536000,c(5:61)]

x_all_2 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_chain2", test_form="empirical", 
                           ignored_pars=NA, spline_burnin=50000, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=250000)

x_all_2.2 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                             cov_matrix=NA, population=1e7, form="peak_trunc", tag="more_knots_trunc_chain2", test_form="empirical", 
                             ignored_pars=NA, spline_burnin=0, 
                             parameters=NA, 
                             day_agg="weekly", cov_start=110000)
x_all_2.2$parameters <- readRDS("all/2023-03-25_all_186000_peak_trunc_MH_empirical_ignored_NA_more_knots_trunc_chain2.rds")$MCMC_output[186000,]


x_all_2.3 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=24, proposal_type="MH", 
                             cov_matrix=NA, population=1e7, form="peak_trunc", tag="more_knots_trunc_chain2", test_form="empirical", 
                             ignored_pars=NA, spline_burnin=50000, 
                             parameters=NA, 
                             day_agg="weekly", cov_start=250000)
x_all_2.3$parameters[c(2,4)] <- x_all_2.3$parameters[c(2,4)]/10


x_all_3 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_no_spline_burnin", test_form="empirical", 
                           ignored_pars=NA, spline_burnin=0, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=250000)


x_all_4 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_fixed_ab_sigma_chain1", test_form="empirical", 
                           ignored_pars=c(2,4), spline_burnin=50000, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=250000)
x_all_4$parameters[c(2,4)] <- x_all_4$prior_mean[c(2,4)]

x_all_5 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_fixed_ab_sigma_chain2", test_form="empirical", 
                           ignored_pars=c(2,4), spline_burnin=50000, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=250000)
x_all_5$parameters[c(2,4)] <- x_all_5$prior_mean[c(2,4)]

x_all_6 <- input_generator(n_iterations=1e6, strain="all", test_pop=1e7, ncores=16, proposal_type="MH", 
                           cov_matrix=NA, population=1e7, form="peak", tag="more_knots_fixed_ab_sigma_chain2", test_form="empirical", 
                           ignored_pars=c(2,4), spline_burnin=50000, 
                           parameters=NA, 
                           day_agg="weekly", cov_start=250000)
x_all_6$parameters[c(2,4)] <- x_all_6$prior_mean[c(2,4)]*0.1


x2 <- readRDS("all/2023-04-23_all_476000_incidence_Cov_empirical_ignored_NA_lognorm_non-differentiable_non-truncated_skew_normal_ctmax39_maxtestday10.rds")
x2_in <- x2$x
x2_nit <- which(is.na(x2$MCMC_posteriors))[1]-2
x2_in$parameters <- x2$MCMC_output[x2_nit,]
x2_in$cov_matrix <- cov(x2$MCMC_output[376000:x2_nit,])*(2.38^2)/length(x2_in$parameters)

x4 <- readRDS("all/2023-04-23_all_470000_peak_Cov_empirical_ignored_NA_lognorm_non-differentiable_non-truncated_skew_normal_ctmax39_maxtestday10.rds")
x4_in <- x4$x
x4_nit <- which(is.na(x4$MCMC_posteriors))[1]-2
x4_in$parameters <- x4$MCMC_output[x4_nit,]
x4_in$cov_matrix <- cov(x4$MCMC_output[420000:x4_nit,])*(2.38^2)/length(x4_in$parameters)

x5 <- readRDS("all/2023-04-23_all_471000_peak_Cov_empirical_ignored_NA_lognorm_non-differentiable_non-truncated_skew_normal_ctmax34_maxtestday6.rds")
x5_in <- x5$x
x5_nit <- which(is.na(x5$MCMC_posteriors))[1]-2
x5_in$parameters <- x5$MCMC_output[x5_nit,]
x5_in$parameters['inc1'] <- -0.6
x5_in$ignored_pars <- c(9)
x5_in$proposal_type <- "MH"
x5_in$tag <- "ctmax34_maxtestday6_fixed_inc"

x6 <- readRDS("all/2023-04-25_all_153000_peak_Cov_empirical_ignored_9_lognorm_non-differentiable_non-truncated_skew_normal_ctmax34_maxtestday6_fixed_inc2.rds")
x6_in <- x6$x
x6_nit <- which(is.na(x6$MCMC_posteriors))[1]-2
x6_in$parameters <- x6$MCMC_output[x6_nit,]
x6_in$parameters['inc1'] <- -1.5
x6_in$ignored_pars <- c(9)
x6_in$proposal_type <- "Cov"
x6_in$cov_matrix <- cov(x6$MCMC_output[20000:153000,])*(2.38^2)/length(x6_in$parameters)
x6_in$tag <- "ctmax34_maxtestday6_fixed_inc2"

# if(output$x$form == "incidence") t <- seq(0,output$x$max_inf,0.1)
# if(output$x$form == "peak") t <- seq(-14,20, 0.1)
# 
# vl_posterior <- matrix(nrow=length(t), ncol=1000)
# growth_rate <- vector(length=1000)
# gap <- vector(length=1000)
# 
# vl_func2 <- function(output, iteration, t){
#   cti <- output$x$ct_range[2]+1
#   ct_min <- rnorm(n=1, mean=output$MCMC_output[iteration,5], sd=output$MCMC_output[iteration,6])
#   a=-100
#   b=-100
#   
#   if(output$x$vl_truncation=="truncated"){
#     if(output$x$vl_growth_type=="lognorm"){
#       while(a<(cti-ct_min)/14) a <- exp(rnorm(n=1, mean=output$MCMC_output[iteration,1], sd=output$MCMC_output[iteration,2]))
#       while(b<(cti-ct_min)/21) b <- exp(rnorm(n=1, mean=output$MCMC_output[iteration,3], sd=output$MCMC_output[iteration,4]))
#       l <- rnorm(n=1,output$MCMC_output[iteration,7],output$MCMC_output[iteration,8])
#     }
#     if(output$x$vl_growth_type=="norm"){
#       while(a<(cti-ct_min)/14) a <- rnorm(n=1, mean=output$MCMC_output[iteration,1], sd=output$MCMC_output[iteration,2])
#       while(b<(cti-ct_min)/21) b <- rnorm(n=1, mean=output$MCMC_output[iteration,3], sd=output$MCMC_output[iteration,4])
#       l <- rnorm(n=1,output$MCMC_output[iteration,7],output$MCMC_output[iteration,8])
#     }
#     if(output$x$form == "peak"){
#       if(output$x$vl_kinetics_type=="non-differentiable") ct <- ifelse(t<=-l/2, ct_min-(t+l/2)*a, ifelse(t>=l/2,ct_min+(t-l/2)*b, ct_min))
#       if(output$x$vl_kinetics_type=="differentiable") ct <- ct_min-log((a+b)/(b*exp(-a*t)+a*exp(b*t)))
#     }
#     if(output$x$form == "incidence"){
#       tau <- (cti-ct_min)/a
#       if(output$x$vl_kinetics_type=="non-differentiable") ct <- ifelse(t<=tau, cti-t*a, ifelse(t>=tau+l,ct_min+(t-tau-l)*b, ct_min))
#       if(output$x$vl_kinetics_type=="differentiable") ct <- ct_min-log((a+b)/(b*exp(-a*(t-tau))+a*exp(b*(t-tau))))
#     }
#   }
#   if(output$x$vl_truncation!="truncated"){
#     if(output$x$vl_growth_type=="lognorm"){
#       a <- exp(rnorm(n=1, mean=output$MCMC_output[iteration,1], sd=output$MCMC_output[iteration,2]))
#       b <- exp(rnorm(n=1, mean=output$MCMC_output[iteration,3], sd=output$MCMC_output[iteration,4]))
#       l <- rnorm(n=1,output$MCMC_output[iteration,7],output$MCMC_output[iteration,8])
#     }
#     if(output$x$vl_growth_type=="norm"){
#       a <- rnorm(n=1, mean=output$MCMC_output[iteration,1], sd=output$MCMC_output[iteration,2])
#       b <- rnorm(n=1, mean=output$MCMC_output[iteration,3], sd=output$MCMC_output[iteration,4])
#       l <- rnorm(n=1,output$MCMC_output[iteration,7],output$MCMC_output[iteration,8])
#     }
#     if(output$x$form == "peak"){
#       if(output$x$vl_kinetics_type=="non-differentiable") ct <- ifelse(t<=-l/2, ct_min-(t+l/2)*a, ifelse(t>=l/2,ct_min+(t-l/2)*b, ct_min))
#       if(output$x$vl_kinetics_type=="differentiable") ct <- ct_min-log((a+b)/(b*exp(-a*t)+a*exp(b*t)+exp(l)))
#     }
#     if(output$x$form == "incidence"){
#       tau <- (cti-ct_min)/a
#       if(output$x$vl_kinetics_type=="non-differentiable") ct <- ifelse(t<=tau, cti-t*a, ifelse(t>=tau+l,ct_min+(t-tau-l)*b, ct_min))
#       if(output$x$vl_kinetics_type=="differentiable") ct <- ct_min-log((a+b)/(b*exp(-a*(t-tau))+a*exp(b*(t-tau))+exp(l)))
#     }
#   }
#   
#   ct[ct>cti] <- cti
#   
#   return(list(a=a,b=b,ct=ct, l=l))
# }
# 
# for(i in 1:ncol(vl_posterior)){
#   iteration <- floor(seq(warmup,n_it-1,length.out=1000)[i])
#   simul <- vl_func2(output, iteration, t)
#   vl_posterior[,i] <- simul$ct
#   growth_rate[i] <- simul$a
#   gap[i] <- simul$l
# } 
# 
# vl_posterior_df <- data.frame(day=t,
#                               rowQuantiles(vl_posterior, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
#   magrittr::set_colnames(c("day", "lower2", "lower1", "ct", "upper1", "upper2"))
# 
# vl_posterior_plot <- ggplot(vl_posterior_df, aes(x=day)) +
#   geom_line(aes(y=ct)) +
#   geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.2) +
#   geom_ribbon(aes(ymin=lower2, ymax=upper2), alpha=0.2) +
#   labs(x="day of infection") +
#   theme_bw() +
#   ggtitle("Ct value by dy of infection")
# 
# print(vl_posterior_plot)
# 

wt_MH <- readRDS("wt/2023-06-08_wt_249000_incidence_MH_empirical_ignored_8_lognorm_non-differentiable_truncated_skew_normal_trunc_diff.rds")
input_wt <- wt_MH$x
input_wt$parameters <- wt_MH$MCMC_output[249000,]
input_wt$proposal_type <- "Cov"
input_wt$cov_matrix <- cov(wt_MH$MCMC_output[200000:249000,])*(2.38^2)/length(input_wt$parameters)
input_wt$tag <- ""

input_wt2 <- input_wt
input_wt2$ignored_pars <- c(7,8)
input_wt2$parameters[7] <- 0

input_wt3 <- input_wt
input_wt3$vl_truncation <- "not_truncated"

input_wt4 <- input_wt
input_wt4$vl_truncation <- "non_truncated"
input_wt4$ignored_pars <- c(7,8)
input_wt4$parameters[7] <- 0


delta_MH <- readRDS("delta/2023-06-08_delta_249000_incidence_MH_empirical_ignored_8_lognorm_non-differentiable_truncated_skew_normal_trunc_diff.rds")
input_delta <- delta_MH$x
input_delta$parameters <- delta_MH$MCMC_output[249000,]
input_delta$proposal_type <- "Cov"
input_delta$cov_matrix <- cov(delta_MH$MCMC_output[200000:249000,])*(2.38^2)/length(input_delta$parameters)
input_delta$tag <- ""

input_delta2 <- input_delta
input_delta2$ignored_pars <- c(7,8)
input_delta2$parameters[7] <- 0

input_delta3 <- input_delta
input_delta3$vl_truncation <- "not_truncated"

input_delta4 <- input_delta
input_delta4$vl_truncation <- "non_truncated"
input_delta4$ignored_pars <- c(7,8)
input_delta4$parameters[7] <- 0


omicron_MH <- readRDS("omicron/2023-06-08_omicron_248000_incidence_MH_empirical_ignored_8_lognorm_non-differentiable_truncated_skew_normal_trunc_diff.rds")
input_omicron <- omicron_MH$x
input_omicron$parameters <- omicron_MH$MCMC_output[248000,]
input_omicron$proposal_type <- "Cov"
input_omicron$cov_matrix <- cov(omicron_MH$MCMC_output[200000:248000,])*(2.38^2)/length(input_omicron$parameters)
input_omicron$tag <- ""

input_omicron2 <- input_omicron
input_omicron2$ignored_pars <- c(7,8)
input_omicron2$parameters[7] <- 0

input_omicron3 <- input_omicron
input_omicron3$vl_truncation <- "not_truncated"

input_omicron4 <- input_omicron
input_omicron4$vl_truncation <- "non_truncated"
input_omicron4$ignored_pars <- c(7,8)
input_omicron4$parameters[7] <- 0

plotting_func2 <- function(strain, tag="", warmup_raw=NULL, n_samples=100, vl_int_grad=c(13.698,0.328)){
  output <- output_extractor(strain, tag)
  colnames(output$MCMC_output)[c(5,6)] <- c("ct_min_bar", "ct_min_sigma")
  output$x$ct_range <- c(10,39)
  
  vg <- ifelse(output$x$vl_growth_type == "lognorm", 1, 0)
  vk <- ifelse(output$x$vl_kinetics_type == "non-differentiable", 1, 0)
  vt <- ifelse(output$x$vl_truncation == "truncated", 1, 0)
  
  colnames(output$MCMC_output) <- names(output$x$prior_mean)
  
  sourceCpp("will3p_vlmax_all_tstoch.cpp")
  sourceCpp("will3p_inc.cpp")
  
  source("ct_symptoms_key_functions.R")
  
  ncores <- 8
  
  n_it <- which(is.na(output$MCMC_posteriors))[1]-2
  if(is.na(n_it)) n_it <- output$x$n_iterations
  
  #plot(output$MCMC_posteriors[1:n_it], type="l")
  
  if(is.null(warmup_raw)==F) warmup <- warmup_raw
  if(is.null(warmup_raw)==T) warmup <- burnin_calculator(output, nroll=10000, n_it)
  
  print(warmup)
  #plot(output$MCMC_posteriors[warmup:n_it], type="l")
  
  ct_range <- output$x$ct_range
  
  symp_delay_lim <- output$x$symp_delay_lim
  if(grepl("peak", output$x$form)) symp_days <- -symp_delay_lim:symp_delay_lim
  if(grepl("incidence", output$x$form)) symp_days <- 0:symp_delay_lim
  
  ## DIC / AIC
  ## Deviance infomration criterion
  #ggplot(data.frame(post=output$MCMC_posteriors[warmup:n_it]), aes(x=post)) + geom_density()
  
  DIC_table <- DIC_table_func(output, warmup, n_it, vg, vk, vt)
  
  ## traceplots
  traceplot <- traceplot_func(tag_list=tag, strain, warmup)
  
  ## PRIOR vs POSTERIOR
  prior_plot <- prior_plot_func(tag_list=tag, strain)
  
  data_fits_plot <- data_fits_plot_func(output, warmup, n_it, n_samples)
  
  ## Overall counts
  vl_posterior_plot <- viral_kinetics_plot_func(tag, strain, n_samples)
  
  plot_inc <- inc_period_plot_func(tag, n_samples)
  
  tables <- output_comparison_table(tag, strain, vl_int_grad=vl_int_grad, n_samples)
  
  top <- ggarrange(data_fits_plot$plot_ct, data_fits_plot$plot_ntest, data_fits_plot$plot_npos, ncol=3)
  second <- ggarrange(data_fits_plot$plot_ct_day, plot_inc, vl_posterior_plot, ncol=3)
  third <- ggarrange(data_fits_plot$plot_symp_ct, data_fits_plot$plot3, ncol=2)
  fourth <- ggarrange(traceplot, prior_plot, ncol=2)
  fifth <- ggarrange(tables$table1, tables$table2, ncol=2)
  
  output_plot1 <- annotate_figure(ggarrange(top, second, third, fourth, nrow=4, heights=c(1, 1, 1, 1)), 
                                  top = text_grob(paste0(strain, ": mean posterior = ",round(mean(output$MCMC_posteriors[warmup:(n_it-2)],0)), "; n_it = ", n_it),color="red", size=14, face="bold"), 
                                  bottom=text_grob(paste("file = ", tag),color = "black", size = 14, hjust=1))
  
  output_plot <- annotate_figure(ggarrange(top, second, third, fourth, fifth, nrow=5, heights=c(1, 1, 1, 1, 1)), 
                                 top = text_grob(paste0(strain, ": mean posterior = ",round(mean(output$MCMC_posteriors[warmup:(n_it-2)],0)), "; n_it = ", n_it),color="red", size=14, face="bold"), 
                                 bottom=text_grob(paste("file = ", tag),color = "black", size = 14, hjust=1))
  
  
  ggsave(filename=paste0(output$x$strain,"/",output$x$strain,"_",output$x$gene,"_",output$x$form,"_",output$x$vl_kinetics,"_",output$x$vl_truncation,"_ignored_",paste(output$x$ignored_pars,sep="",collapse=","),"_",output$x$ct_range[2],"_2.png"), plot = output_plot1, device = "png", width=15, height=24)
  ggsave(filename=paste0(output$x$strain,"/",output$x$strain,"_",output$x$gene,"_",output$x$form,"_",output$x$vl_kinetics,"_",output$x$vl_truncation,"_ignored_",paste(output$x$ignored_pars,sep="",collapse=","),"_",output$x$ct_range[2],"_2.pdf"), output_plot, width=16, height = 25)
  
  #return(time_vec)
}

