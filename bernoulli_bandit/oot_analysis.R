## ----------------------------------------------------------------- ##
## oot_analysis.R -------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: evaluate out of trial performance ---------------------- ##
## ----------------------------------------------------------------- ##

library(dplyr)

setwd("-------")
source("general_functions.R")


## oot_rate
## Purpose: evalue how often we assign optimal trt based off
##          of trial dataset, with different subject distributions
##          that may change
## param fit: glm fit from trial data
## param x_shape1: first parameter X ~ beta(shape1,shape2)
## param x_shape2: second parameter X ~ beta(shape1,shape2)
## param N: number of subjects to treat
## return info: info with success rate, opt trt rate, and regret
oot_rate <- function(fit,x_shape1,x_shape2,N=500){
  
  ## generate random vector of the x cov
  x1_vec <- rbeta(N,x_shape1,x_shape2)
  ## empty vector for treatments
  A_vec <- c()
  ## parameter estimates from fit
  beta_est <- fit$coefficients
  cov_est <- vcov(fit)
  
  ## looping through all individuals to find estimated opt trt
  for(i in 1:N){
    a_info <- matrix(NA,2,2)
    for(a in c(0,1)){
      x1 <- x1_vec[i]
      a_info[a+1,1] <- a
      a_info[a+1,2] <- c(1,a,x1,x1**2,a*x1,a*(x1**2)) %*% 
                       as.matrix(beta_est)
    }
    A_vec[i] <-  a_info[which.max(a_info[,2]),1]
  }
    
  ## creating a dataframe
  dat <- data.frame(sub=1:N,x1=x1_vec,A=A_vec)
  ## calculating log odds
  dat$log_odds <- mean_logit(A=dat$A,x=dat$x1,
                             beta0=0.00,beta1=-0.75,beta2=1.50,
                             beta3=-3.00,beta4=-3.25,beta5=6.67)
  ## probability of success
  dat$p <- get_prob(dat$log_odds)
  ## response
  dat$Y <- NA
  for(i in 1:N){
    dat$Y[i] <- rbinom(1,1,dat$p[i])
  }
  ## what the optimal treatment is
  dat$A_opt <- ifelse(dat$x1>0.66,1,0)
  ## how often we getit correct
  dat$correct <- ifelse(dat$A_opt==dat$A,1,0)
  ## calculating the regret for all individuals
  dat$regret <- mean_logit(A=dat$A_opt,x=dat$x1,
                beta0=0.00,beta1=-0.75,beta2=1.50,
                beta3=-3.00,beta4=-3.25,beta5=6.67) - dat$log_odds
  

  
  info <- c(mean(dat$Y),mean(dat$correct),
            mean(dat$regret),mean(dat[dat$correct==0,]$regret))

  
  return(info)
  
}  


## oot_analysis
## Purpose: take a dataset wtih multiple trials and extract
##          information run oot_rate on each one
## param trials: dataframe with multiple trials
## param method_str: string of different methods to evaluate
## param num_reps: number of replications per method
## param num_sub: number of subjects in each trial to consider
## param x_shape1: first parameter X ~ beta(shape1,shape2)
## param x_shape2: second parameter X ~ beta(shape1,shape2)
## param N_oot: number of subjects outside the trial to evaluate
## return big_info: info returned from oot_rate for each replication
oot_analysis <- function(trials,method_str,num_reps,num_sub,
                         x_shape1,x_shape2,N_oot=500){
  
  df <- trials %>% filter(method==method_str,subject<=num_sub)
  big_info <- matrix(NA,nrow=num_reps,ncol=4)
  for(i in 1:num_reps){
    temp <- df %>% filter(rep==i)
    fit <- glm(Y~A+x1+I(x1**2)+x1:A+I(x1**2):A,family=binomial(link="logit"),
               temp)
    big_info[i,] <- oot_rate(fit=fit,x_shape1=x_shape1,x_shape2=x_shape2,N=N_oot)
  }
    
  big_info <- as.data.frame(big_info)  
  colnames(big_info) <- c("success","correct","regret_0","regret_no0")
  
  big_info$rep <- 1:num_reps
  big_info$method <- method_str
  big_info$N <- num_sub
  
  return(big_info)
  
}


## power_calc
## Purpose: run hypothesis tests based on one glm fit
## param fit: GLM fit from a trial
## return big_info: p-value and rejection indicator (alpha=0.05) for each test
power_calc <- function(fit){
  
  beta_est <- fit$coefficients
  cov_est <- vcov(fit)
  
  ## Conducting three different hypothesis tests
  c1 <- c(0,1,0,0,0.6585,0.6585**2) ## h0 true
  t1 <- c1 %*% beta_est
  se1 <- sqrt(t(c1) %*% cov_est %*% c1)
  p1 <- 1-pnorm(t1/se1)
  r1 <- ifelse(p1<0.05,1,0)
  
  c2 <- c(0,1,0,0,0.70,0.77**2) ## small trt effect
  t2 <- c2 %*% beta_est
  se2 <- sqrt(t(c2) %*% cov_est %*% c2)
  p2 <- 1-pnorm(t2/se2)
  r2 <- ifelse(p2<0.05,1,0)
  
  c3 <- c(0,1,0,0,0.90,0.90**2) ## large trt effect
  t3 <- c3 %*% beta_est
  se3 <- sqrt(t(c3) %*% cov_est %*% c3)
  p3 <- 1-pnorm(t3/se3)
  r3 <- ifelse(p3<0.05,1,0)
  
  info <- c(p1,r1,p2,r2,p3,r3)
  return(info)
}

## power_analysis
## Purpose: run power analysis for a number of trials
## param trials: dataframe with multiple trials
## param method_str: string of different methods to evaluate
## param num_reps: number of replications per method
## param num_sub: number of subjects in each trial to consider
## return big_info: info from power_calc for each trial
power_analysis <- function(trials,method_str,num_reps,num_sub){
  
  df <- trials %>% filter(method==method_str,subject<=num_sub)
  big_info <- matrix(NA,nrow=num_reps,ncol=6)
  for(i in 1:num_reps){
    temp <- df %>% filter(rep==i)
    fit <- glm(Y~A+x1+I(x1**2)+x1:A+I(x1**2):A,family=binomial(link="logit"),
               temp)
    big_info[i,] <- power_calc(fit=fit)
  }
  
  big_info <- as.data.frame(big_info)  
  colnames(big_info) <- c("p1","r1","p2","r2","p3","r3")
  
  big_info$rep <- 1:num_reps
  big_info$method <- method_str
  big_info$N <- num_sub
  
  return(big_info)
  
}
