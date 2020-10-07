## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the anova/mab scenario ------------- ##
## ----------------------------------------------------------------- ##

setwd("-------")
source("general_functions.R")
source("simple.R")
source("det_simple.R")
source("thompson.R")
source("det_thompson.R")

library(parallel)

## wrapper function
sim <- function(burn_in,N,
                p1=0.2,p2=0.2,p3=0.2,p4=0.2,p5=0.2,
                mu1,mu2,mu3,mu4,mu5,
                sigma,
                B=500,c=0.5,br=100){

  
  
  ## generate training set
  train_set <- gen_train(N=burn_in,
                      p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,
                      mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,
                      sigma=sigma)
  
  simple_sim <- try( simple(train_set=train_set,burn_in=burn_in,N=N,
                        p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,
                        mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                        sigma=sigma) )
  
  det_simple_sim <- try( det_simple(train_set=train_set,burn_in=burn_in,N=N,
                                         mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                                         sigma=sigma,br=br) )
  
  thompson_c_sim <- try( thompson(train_set=train_set,burn_in=burn_in,N=N,
                                  mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                                  sigma=sigma,B=B,c=c) )
  
  thompson_1_sim <- try( thompson(train_set=train_set,burn_in=burn_in,N=N,
                                  mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                                  sigma=sigma,B=B,c=1) )
  
  det_thompson_c_sim <- try( det_thompson(train_set=train_set,burn_in=burn_in,N=N,
                                mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                                sigma=sigma,B=B,c=c,br=br) )
  
  det_thompson_1_sim <- try( det_thompson(train_set=train_set,burn_in=burn_in,N=N,
                                          mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5,
                                          sigma=sigma,B=B,c=1,br=br) )
  
  
  
  output <- try( list(simple=simple_sim,
                      det_simple=det_simple_sim,
                      thompson_c=thompson_c_sim,
                      thompson_1=thompson_1_sim,
                      det_thompson_c=det_thompson_c_sim,
                      det_thompson_1=det_thompson_1_sim))
  return(output)
  
}


## how many simulations to run
r <- 5000

burn_in=7
N=500
## Simulate the trials
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(burn_in=7,N=N,
                                 p1=0.2,p2=0.2,p3=0.2,p4=0.2,p5=0.2,
                                 mu1=1,mu2=1,mu3=1,mu4=1.8,mu5=2.0,
                                 sigma=1,
                                 B=500,c=0.5,br=100)},
                 mc.cores = 8
)
end <- Sys.time()


## compile the data into nice dataframes
trials <- data.frame()
param_ests <- data.frame()
std_errors <- data.frame()
rand_probs <- data.frame()
for(i in 1:(r)){
  
  ## trial information
  s_simple <- sims[[i]]$simple[[1]]
  s_simple$subject <- 1:nrow(s_simple)
  s_simple$method <- "simple"
  s_simple$rep <- i
  
  s_det_simple <- sims[[i]]$det_simple[[1]]
  s_det_simple$subject <- 1:nrow(s_det_simple)
  s_det_simple$method <- "det_simple"
  s_det_simple$rep <- i
  
  s_thompson_c <- sims[[i]]$thompson_c[[1]]
  s_thompson_c$subject <- 1:nrow(s_thompson_c)
  s_thompson_c$method <- "thompson_c"
  s_thompson_c$rep <- i
  
  s_thompson_1 <- sims[[i]]$thompson_1[[1]]
  s_thompson_1$subject <- 1:nrow(s_thompson_1)
  s_thompson_1$method <- "thompson_1"
  s_thompson_1$rep <- i
  
  s_det_thompson_c <- sims[[i]]$det_thompson_c[[1]]
  s_det_thompson_c$subject <- 1:nrow(s_det_thompson_c)
  s_det_thompson_c$method <- "det_thompson_c"
  s_det_thompson_c$rep <- i
  
  s_det_thompson_1 <- sims[[i]]$det_thompson_1[[1]]
  s_det_thompson_1$subject <- 1:nrow(s_det_thompson_1)
  s_det_thompson_1$method <- "det_thompson_1"
  s_det_thompson_1$rep <- i
  
  
  trials <- rbind(trials,s_simple,s_det_simple,
                  s_thompson_c,s_thompson_1,
                  s_det_thompson_c,s_det_thompson_1)
  
  ## parameter estimates
  est_simple <- sims[[i]]$simple[[2]] %>% as.data.frame()
  est_simple$subject <- 1:nrow(est_simple)
  est_simple$method <- "simple"
  est_simple$rep <- i
  
  est_det_simple <- sims[[i]]$det_simple[[2]] %>% as.data.frame()
  est_det_simple$subject <- 1:nrow(est_det_simple)
  est_det_simple$method <- "det_simple"
  est_det_simple$rep <- i
  
  est_thompson_c <- sims[[i]]$thompson_c[[2]] %>% as.data.frame()
  est_thompson_c$subject <- 1:nrow(est_thompson_c)
  est_thompson_c$method <- "thompson_c"
  est_thompson_c$rep <- i
  
  est_thompson_1 <- sims[[i]]$thompson_1[[2]] %>% as.data.frame()
  est_thompson_1$subject <- 1:nrow(est_thompson_1)
  est_thompson_1$method <- "thompson_1"
  est_thompson_1$rep <- i
  
  est_det_thompson_c <- sims[[i]]$det_thompson_c[[2]] %>% as.data.frame()
  est_det_thompson_c$subject <- 1:nrow(est_det_thompson_c) 
  est_det_thompson_c$method <- "det_thompson_c"
  est_det_thompson_c$rep <- i
  
  est_det_thompson_1 <- sims[[i]]$det_thompson_1[[2]] %>% as.data.frame()
  est_det_thompson_1$subject <- 1:nrow(est_det_thompson_1) 
  est_det_thompson_1$method <- "det_thompson_1"
  est_det_thompson_1$rep <- i
  
  param_ests <- rbind(param_ests,est_simple,est_det_simple,
                      est_thompson_c,est_thompson_1,
                      est_det_thompson_c,est_det_thompson_1)
  
  ## standard errors
  std_simple <- sims[[i]]$simple[[3]] %>% as.data.frame()
  std_simple$subject <- 1:nrow(std_simple)
  std_simple$method <- "simple"
  std_simple$rep <- i
  
  std_det_simple <- sims[[i]]$det_simple[[3]] %>% as.data.frame()
  std_det_simple$subject <- 1:nrow(std_det_simple)
  std_det_simple$method <- "det_simple"
  std_det_simple$rep <- i
  
  std_thompson_c <- sims[[i]]$thompson_c[[3]] %>% as.data.frame()
  std_thompson_c$subject <- 1:nrow(std_thompson_c)
  std_thompson_c$method <- "thompson_c"
  std_thompson_c$rep <- i
  
  std_thompson_1 <- sims[[i]]$thompson_1[[3]] %>% as.data.frame()
  std_thompson_1$subject <- 1:nrow(std_thompson_1)
  std_thompson_1$method <- "thompson_1"
  std_thompson_1$rep <- i
  
  std_det_thompson_c <- sims[[i]]$det_thompson_c[[3]] %>% as.data.frame()
  std_det_thompson_c$subject <- 1:nrow(std_det_thompson_c) 
  std_det_thompson_c$method <- "det_thompson_c"
  std_det_thompson_c$rep <- i
  
  std_det_thompson_1 <- sims[[i]]$det_thompson_1[[3]] %>% as.data.frame()
  std_det_thompson_1$subject <- 1:nrow(std_det_thompson_1) 
  std_det_thompson_1$method <- "det_thompson_1"
  std_det_thompson_1$rep <- i
  
  std_errors <- rbind(std_errors,std_simple,std_det_simple,
                      std_thompson_c,std_thompson_1,
                      std_det_thompson_c,std_det_thompson_1)
  
  ## randomization probabilities
  rand_simple <- sims[[i]]$simple[[4]] %>% as.data.frame()
  rand_simple$subject <- 1:nrow(rand_simple)
  rand_simple$method <- "simple"
  rand_simple$rep <- i
  
  rand_det_simple <- sims[[i]]$det_simple[[4]] %>% as.data.frame()
  rand_det_simple$subject <- 1:nrow(rand_det_simple)
  rand_det_simple$method <- "det_simple"
  rand_det_simple$rep <- i
  
  rand_thompson_c <- sims[[i]]$thompson_c[[4]] %>% as.data.frame()
  rand_thompson_c$subject <- 1:nrow(rand_thompson_c)
  rand_thompson_c$method <- "thompson_c"
  rand_thompson_c$rep <- i
  
  rand_thompson_1 <- sims[[i]]$thompson_1[[4]] %>% as.data.frame()
  rand_thompson_1$subject <- 1:nrow(rand_thompson_1)
  rand_thompson_1$method <- "thompson_1"
  rand_thompson_1$rep <- i
  
  rand_det_thompson_c <- sims[[i]]$det_thompson_c[[4]] %>% as.data.frame()
  rand_det_thompson_c$subject <- 1:nrow(rand_det_thompson_c) 
  rand_det_thompson_c$method <- "det_thompson_c"
  rand_det_thompson_c$rep <- i
  
  rand_det_thompson_1 <- sims[[i]]$det_thompson_1[[4]] %>% as.data.frame()
  rand_det_thompson_1$subject <- 1:nrow(rand_det_thompson_1) 
  rand_det_thompson_1$method <- "det_thompson_1"
  rand_det_thompson_1$rep <- i
  
  rand_probs <- rbind(rand_probs,rand_simple,rand_det_simple,
                      rand_thompson_c,rand_thompson_1,
                      rand_det_thompson_c,rand_det_thompson_1)
  
}

## save the information
save(trials,file="--------.RData")
save(rand_probs,file="-----------.RData")
save(param_ests,file="-----------.RData")
save(std_errors,file="-----------.RData")

