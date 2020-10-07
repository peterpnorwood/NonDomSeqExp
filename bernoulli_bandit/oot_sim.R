## ----------------------------------------------------------------- ##
## oot_sim.R ------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: out of trial analysis simulation ----------------------- ##
## ----------------------------------------------------------------- ##


library(dplyr)

setwd("-------")
source("oot_analysis.R")

## load the trials
load("-------.RData")

## fix the parameters as we would like
x_shape1=1
x_shape2=1

## methods and sample sizes we loop through
methods <- c("simple","kl_simple","e_greedy_25","e_greedy_10","kl_e_greedy_25","kl_e_greedy_10")
Ns <- c(50,100,250,500)

## loop through
big_df <- data.frame()
for(n in Ns){
  lapp <- lapply(methods,oot_analysis,
                 trials=trials,num_reps=max(trials$rep),
                 num_sub=n,x_shape1=x_shape1,x_shape2=x_shape2,N_oot=500)
  df <- data.frame()
  for(j in 1:length(lapp)){
    temp <- lapp[[j]]
    df <- rbind(df,temp)
  }
  big_df <- rbind(big_df,df)
}

big_df$dist <- "----"
save(big_df,file="-----.RData")


### power analysis
power_df <- data.frame()
for(n in Ns){
  lapp <- lapply(methods,power_analysis,
                 trials=trials,num_reps=max(trials$rep),
                 num_sub=n)
  df <- data.frame()
  for(j in 1:length(lapp)){
    temp <- lapp[[j]]
    df <- rbind(df,temp)
  }
  power_df <- rbind(power_df,df)
}

save(power_df,file="------.RData")