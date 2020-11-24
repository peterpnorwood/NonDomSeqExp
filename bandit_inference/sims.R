## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: functions to rum sims ---------------------------------- ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("general_functions.R")
source("SR.R")
source("RSR.R")
source("TS.R")
source("RTS.R")

library(parallel)

## wrapper function
sim <- function(burn_in,theta0,theta1,sigma,
                batch_size=50,batches=10,
                eps=0.1,clip=0.1){
  
  
  ## generate training set
  train_set <- gen_data(burn_in,theta0,theta1,sigma)
  
  SR_sim <- try(SR(train_set=train_set,
                    theta0=theta0,theta1=theta1,sigma=sigma,
                    batches=batches,batch_size=batch_size))
  
  RSR_sim <- try(RSR(train_set=train_set,
                     theta0=theta0,theta1=theta1,sigma=sigma,
                     batches=batches,batch_size=batch_size))
  
  TS_sim <- try(TS(train_set=train_set,
                   theta0=theta0,theta1=theta1,sigma=sigma,
                   batches=batches,batch_size=batch_size,clip=clip))
  
  RTS_sim <- try(RTS(train_set=train_set,
                     theta0=theta0,theta1=theta1,sigma=sigma,
                     batches=batches,batch_size=batch_size,clip=clip))
  
  
  output <- try(list(SR=SR_sim,
                     RSR=RSR_sim,
                     TS=TS_sim,
                     RTS=RTS_sim))
  return(output)
  
}

## another wrapper function
run_sim <- function(burn_in,theta0,
                    theta1,sigma,
                    batch_size,batches,
                    eps,clip,r){
  
  ## Simulate the trials
  sims <- mclapply(X=1:r, 
                   function(X){sim(burn_in=burn_in,
                                   theta0=theta0,theta1=theta1,sigma=sigma,
                                   batch_size=batch_size,batches=batches,
                                   eps=eps,clip=clip)},
                   mc.cores = 1)
  
  
  ## throw out bad sets
  bad_sets <- c()
  tick=1
  for(i in 1:r){
    for(j in 1:4){
      if(!is.list(sims[[i]][[j]])){
        bad_sets[tick]=i
        tick=tick+1
      }else{
        next
      }
    }
  }
  good_set <- c(1:r)
  good_set[bad_sets] <- NA
  good_set <- na.omit(good_set)
  
  ## compile the data into nice dataframes
  trials <- data.frame()
  pvals <- data.frame()
  for(i in good_set){
    
    ## trial information
    s_SR <- sims[[i]]$SR$trial
    s_SR$method <- "SR"
    s_SR$rep <- i
    
    s_RSR <- sims[[i]]$RSR$trial
    s_RSR$method <- "RSR"
    s_RSR$rep <- i
    
    s_TS <- sims[[i]]$TS$trial
    s_TS$method <- "TS"
    s_TS$rep <- i
    
    s_RTS <- sims[[i]]$RTS$trial
    s_RTS$method <- "RTS"
    s_RTS$rep <- i
    
    trials <- rbind(trials,s_SR,s_RSR,s_TS,s_RTS)
    
    ## pvalues
    
    ## trial information
    p_SR <- sims[[i]]$SR$pval_df
    p_SR$rand_method <- "SR"
    p_SR$rep <- i
    
    p_RSR <- sims[[i]]$RSR$pval_df
    p_RSR$rand_method <- "RSR"
    p_RSR$rep <- i
    
    
    p_TS <- sims[[i]]$TS$pval_df
    p_TS$rand_method <- "TS"
    p_TS$rep <- i
    
    p_RTS <- sims[[i]]$RTS$pval_df
    p_RTS$rand_method <- "RTS"
    p_RTS$rep <- i
    
    pvals <- rbind(pvals,p_SR,p_RSR,p_TS,p_RTS)
    
  }
  
  pvals$theta1=theta1
  trials$theta1=theta1
  
  return(list(pvals=pvals,
              trials=trials))
  
}
