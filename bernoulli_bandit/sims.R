## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the logit reg example -------------- ##
## ----------------------------------------------------------------- ##

setwd("------")
source("general_functions.R")
source("simple.R")
source("kl_simple.R")
source("e_greedy.R")
source("kl_e_greedy.R")


library(parallel)
library(dplyr)
library(xtable)
library(ggplot2)

## wrapper function
sim <- function(N,burn_in,p1=0.5,
                beta0,beta1,beta2,
                beta3,beta4,beta5,
                x1_shape1=1,x1_shape2=2.5){
  
  
  
  ## generate training set
  train_set <- gen_train(N=N+5,p1=0.5,
                         beta0=beta0,beta1=beta1,beta2=beta2,
                         beta3=beta3,beta4=beta4,beta5=beta5,
                         x1_shape1=1,x1_shape2=2.5)
  
  simple_sim <- try( simple(train_set=train_set,burn_in=burn_in,N=N,
                            beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                            beta4=beta4,beta5=beta5) )
  
  kl_simple_sim <- try( kl_simple(train_set=train_set,burn_in=burn_in,N=N,
                               beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                               beta4=beta4,beta5=beta5) )
  
  e_greedy_sim_25 <- try( e_greedy(train_set=train_set,burn_in=burn_in,N=N,
                              beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                              beta4=beta4,beta5=beta5,eps=0.25) )
  
  e_greedy_sim_10 <- try( e_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                   beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                                   beta4=beta4,beta5=beta5,eps=0.10) )
  
  kl_e_greedy_sim_25 <- try( kl_e_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                   beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                                   beta4=beta4,beta5=beta5,eps=0.25) )
  
  kl_e_greedy_sim_10 <- try( kl_e_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                   beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                                   beta4=beta4,beta5=beta5,eps=0.10) )
  

  
  output <- try( list(simple=simple_sim,
                      kl_simple=kl_simple_sim,
                      e_greedy_25=e_greedy_sim_25,
                      e_greedy_10=e_greedy_sim_10,
                      kl_e_greedy_25=kl_e_greedy_sim_25,
                      kl_e_greedy_10=kl_e_greedy_sim_10))
  return(output)
  
}


## how many simulations to run
r <- 5000

## true parameter values
beta0=0.00
beta1=-0.75
beta2=1.50
beta3=-3.00
beta4=-3.25
beta5=6.67

burn_in=25
N=500
## Simulate the SMARTs
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(N=N,burn_in=burn_in,p1=0.5,
                                 beta0=beta0,beta1=beta1,beta2=beta2,
                                 beta3=beta3,beta4=beta4,beta5=beta5,
                                 x1_shape1=1,x1_shape2=2.5)},
                 mc.cores = 8
)
end <- Sys.time()

## save the information
trials <- data.frame()
param_ests <- data.frame()
std_errors <- data.frame()
rand_probs <- data.frame()

## toss out the simulations which crashed due to numerical
## instability early in the trial (burn in wasn't sufficient)
good_set <- c()
j=1
for(i in 1:r){
  info <- c(
  is.data.frame(sims[[i]]$simple[[1]]),
  is.data.frame(sims[[i]]$kl_simple[[1]]),
  is.data.frame(sims[[i]]$e_greedy_25[[1]]),
  is.data.frame(sims[[i]]$e_greedy_10[[1]]),
  is.data.frame(sims[[i]]$kl_e_greedy_25[[1]]),
  is.data.frame(sims[[i]]$kl_e_greedy_10[[1]]))
  if(mean(info)==1){
    good_set[j] = i
    j=j+1
  }
}

## loop through the "good" trials to save the information
j=1
for(i in good_set){
  
  ## trial information
  s_simple <- sims[[i]]$simple[[1]]
  s_simple$subject <- 1:nrow(s_simple)
  s_simple$method <- "simple"
  s_simple$rep <- j
  
  s_kl_simple <- sims[[i]]$kl_simple[[1]]
  s_kl_simple$subject <- 1:nrow(s_kl_simple)
  s_kl_simple$method <- "kl_simple"
  s_kl_simple$rep <- j
  
  s_e_greedy_25 <- sims[[i]]$e_greedy_25[[1]]
  s_e_greedy_25$subject <- 1:nrow(s_e_greedy_25)
  s_e_greedy_25$method <- "e_greedy_25"
  s_e_greedy_25$rep <- j
  
  s_e_greedy_10 <- sims[[i]]$e_greedy_10[[1]]
  s_e_greedy_10$subject <- 1:nrow(s_e_greedy_10)
  s_e_greedy_10$method <- "e_greedy_10"
  s_e_greedy_10$rep <- j
  
  s_kl_e_greedy_25 <- sims[[i]]$kl_e_greedy_25[[1]]
  s_kl_e_greedy_25$subject <- 1:nrow(s_kl_e_greedy_25)
  s_kl_e_greedy_25$method <- "kl_e_greedy_25"
  s_kl_e_greedy_25$rep <- j
  
  s_kl_e_greedy_10 <- sims[[i]]$kl_e_greedy_10[[1]]
  s_kl_e_greedy_10$subject <- 1:nrow(s_kl_e_greedy_10)
  s_kl_e_greedy_10$method <- "kl_e_greedy_10"
  s_kl_e_greedy_10$rep <- j
  
  
  trials <- rbind(trials,s_simple,s_kl_simple,
                  s_e_greedy_25,s_e_greedy_10,
                  s_kl_e_greedy_25,s_kl_e_greedy_10)
  
  ## parameter estimates
  est_simple <- sims[[i]]$simple[[2]] %>% as.data.frame()
  est_simple$subject <- 1:nrow(est_simple)
  est_simple$method <- "simple"
  est_simple$rep <- j
  
  est_kl_simple <- sims[[i]]$kl_simple[[2]] %>% as.data.frame()
  est_kl_simple$subject <- 1:nrow(est_kl_simple)
  est_kl_simple$method <- "kl_simple"
  est_kl_simple$rep <- j
  
  est_e_greedy_25 <- sims[[i]]$e_greedy_25[[2]] %>% as.data.frame()
  est_e_greedy_25$subject <- 1:nrow(est_e_greedy_25)
  est_e_greedy_25$method <- "e_greedy_25"
  est_e_greedy_25$rep <- j
  
  est_e_greedy_10 <- sims[[i]]$e_greedy_10[[2]] %>% as.data.frame()
  est_e_greedy_10$subject <- 1:nrow(est_e_greedy_10)
  est_e_greedy_10$method <- "e_greedy_10"
  est_e_greedy_10$rep <- j
  
  est_kl_e_greedy_25 <- sims[[i]]$kl_e_greedy_25[[2]] %>% as.data.frame()
  est_kl_e_greedy_25$subject <- 1:nrow(est_kl_e_greedy_25) 
  est_kl_e_greedy_25$method <- "kl_e_greedy_25"
  est_kl_e_greedy_25$rep <- j
  
  est_kl_e_greedy_10 <- sims[[i]]$kl_e_greedy_10[[2]] %>% as.data.frame()
  est_kl_e_greedy_10$subject <- 1:nrow(est_kl_e_greedy_10) 
  est_kl_e_greedy_10$method <- "kl_e_greedy_10"
  est_kl_e_greedy_10$rep <- j
  
  param_ests <- rbind(param_ests,est_simple,est_kl_simple,
                      est_e_greedy_25,est_e_greedy_10,
                      est_kl_e_greedy_25,est_kl_e_greedy_10)
  
  ## standard errors
  std_simple <- sims[[i]]$simple[[3]] %>% as.data.frame()
  std_simple$subject <- 1:nrow(std_simple)
  std_simple$method <- "simple"
  std_simple$rep <- j
  
  std_kl_simple <- sims[[i]]$kl_simple[[3]] %>% as.data.frame()
  std_kl_simple$subject <- 1:nrow(std_kl_simple)
  std_kl_simple$method <- "kl_simple"
  std_kl_simple$rep <- j
  
  std_e_greedy_25 <- sims[[i]]$e_greedy_25[[3]] %>% as.data.frame()
  std_e_greedy_25$subject <- 1:nrow(std_e_greedy_25)
  std_e_greedy_25$method <- "e_greedy_25"
  std_e_greedy_25$rep <- i
  
  std_e_greedy_10 <- sims[[i]]$e_greedy_10[[3]] %>% as.data.frame()
  std_e_greedy_10$subject <- 1:nrow(std_e_greedy_10)
  std_e_greedy_10$method <- "e_greedy_10"
  std_e_greedy_10$rep <- j
  
  std_kl_e_greedy_25 <- sims[[i]]$kl_e_greedy_25[[3]] %>% as.data.frame()
  std_kl_e_greedy_25$subject <- 1:nrow(std_kl_e_greedy_25) 
  std_kl_e_greedy_25$method <- "kl_e_greedy_25"
  std_kl_e_greedy_25$rep <- j
  
  std_kl_e_greedy_10 <- sims[[i]]$kl_e_greedy_10[[3]] %>% as.data.frame()
  std_kl_e_greedy_10$subject <- 1:nrow(std_kl_e_greedy_10) 
  std_kl_e_greedy_10$method <- "kl_e_greedy_10"
  std_kl_e_greedy_10$rep <- j
  
  j=j+1
  
  std_errors <- rbind(std_errors,std_simple,std_kl_simple,
                      std_e_greedy_25,std_e_greedy_10,
                      std_kl_e_greedy_25,std_kl_e_greedy_10)
  
}

## save the information
save(trials,file="------.RData")
save(param_ests,file="------.RData")
save(std_errors,file="------.RData")

