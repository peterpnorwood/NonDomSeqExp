## ----------------------------------------------------------------- ##
## power_one.R ----------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the two arm bandit ----------------- ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("sims.R")

library(parallel)
library(tidyverse)

## parameters 
burn_in=5
theta0=0
sigma=1
batch_size=50
batches=100
eps=0.1
clip=0.1

## reps
r=5250

## wrapper function for different theta1
wrapper <- function(theta1){
  res <- run_sim(burn_in=burn_in,theta0=theta0,
                   theta1=theta1,sigma=sigma,
                   batch_size=batch_size,batches=batches, 
                   eps=eps, clip=clip,r=r)
  
  return(res)
}

## run sims
length_out=10
theta1_list <- seq(0,0.25,length.out=length_out)
results <- mclapply(theta1_list,wrapper,mc.cores=1)

## gather results
power <- data.frame()
trial <- data.frame()
for(i in 1:length_out){
  temp <- results[[i]]$pvals
  temp <- temp %>% filter(method=="bols")
  power <- rbind(power,temp)
  
  temp_trial <- results[[i]]$trials
  trial <- rbind(trial,temp_trial)
}

## save data
power_one <- power %>% 
  filter(rep<=5000) %>%
  mutate(reject=ifelse(pval<0.05,1,0)) %>%
  group_by(rand_method,theta1) %>%
  summarise(power=mean(reject))
  
save(power_one,file="power_one.RData")

in_trial_one <- trial %>%
  filter(rep<=5000) %>% 
  group_by(method,theta1) %>%
  summarise(mean_correct=mean(A),
            mean_reward=mean(Y))

save(in_trial_one,file="in_trial_one.RData")


