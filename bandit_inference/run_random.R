## ----------------------------------------------------------------- ##
## run_random.R ---------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the two arm bandit w/ random ------- ##
## batch sizes ----------------------------------------------------- ##
## ----------------------------------------------------------------- ##

setwd("-------------------------")
source("sims.R")

library(parallel)
library(tidyverse)

## parameters 
burn_in=5
theta0=0
sigma=1
batch_size=NULL
batches=NULL
batch_type="random"
batch_range=3:20
clip=0.1
N_total=1000
## reps
r=10

## wrapper function for different theta1
wrapper <- function(theta1){
  res <- run_sim(r=r,burn_in=burn_in,
                 theta0=theta0,theta1=theta1,sigma=sigma,clip=clip,
                 batch_type=batch_type,
                 batch_size=batch_size,batches=batches,
                 batch_range=batch_range,N_total=N_total)
  
  return(res)
}

## run sims
length_out=5
theta1_list <- seq(0,0.5,length.out=length_out)
results <- mclapply(theta1_list,wrapper,mc.cores=1)
save(results,file="random_results.RData")

## gather results
power <- data.frame()
trial <- data.frame()
for(i in 1:length_out){
  temp <- results[[i]]$pvals
  temp <- temp
  power <- rbind(power,temp)
  
  temp_trial <- results[[i]]$trials
  trial <- rbind(trial,temp_trial)
}

## save power information
random_power <- power %>% 
  filter(rep<=5000) %>%
  mutate(reject=ifelse(pval<0.05,1,0)) %>%
  group_by(method,rand_method,theta1) %>%
  summarise(power=mean(reject))
  
save(random_power,file="random_power.RData")

## save randomization probabilites at the last batch
random_pi <- trial %>%
  filter(rep<=5000 & subject==N_total+1)

save(random_pi,file="random_pi.RData")

## save reward and correct assignment
random_reward <- trial %>%
                 filter(rep<=5000) %>%
                 group_by(method,rep,theta1) %>%
                 summarise(reward=mean(Y),
                           correct=mean(A)) %>%
                 group_by(method,theta1) %>%
                 summarise(mean_reward=mean(reward),
                           sd_reward=sd(reward),
                           mean_correct=mean(correct),
                           sd_correct=sd(correct))

save(random_reward,file="random_reward.RData")

