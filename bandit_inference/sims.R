## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the anova/mab scenario ------------- ##
## ----------------------------------------------------------------- ##

setwd("-------------------------")
source("general_functions.R")
source("SR.R")
source("RSR.R")
source("TS.R")
source("RTS.R")

library(parallel)

## wrapper function
sim <- function(burn_in,
                theta0,theta1,sigma,clip=0.1,
                batch_type,
                batch_size=NULL,batches=NULL,
                batch_range=NULL,N_total=NULL){
  
  
  ## generate training set
  train_set <- gen_data(burn_in,theta0,theta1,sigma)
  
  SR_sim <- try(SR(train_set=train_set,
                   theta0=theta0,theta1=theta1,sigma=sigma,
                   batch_type=batch_type,
                   batch_size=batch_size,batches=batches,
                   batch_range=batch_range,N_total=N_total))
  
  RSR_sim <- try(RSR(train_set=train_set,
                     theta0=theta0,theta1=theta1,sigma=sigma,
                     batch_type=batch_type,
                     batch_size=batch_size,batches=batches,
                     batch_range=batch_range,N_total=N_total))
  
  TS_sim <- try(TS(train_set=train_set,
                   theta0=theta0,theta1=theta1,sigma=sigma,clip=clip,
                   batch_type=batch_type,
                   batch_size=batch_size,batches=batches,
                   batch_range=batch_range,N_total=N_total))
  
  RTS_sim <- try(RTS(train_set=train_set,
                     theta0=theta0,theta1=theta1,sigma=sigma,clip=clip,
                     batch_type=batch_type,
                     batch_size=batch_size,batches=batches,
                     batch_range=batch_range,N_total=N_total))
  
  
  output <- try(list(SR=SR_sim,
                     RSR=RSR_sim,
                     TS=TS_sim,
                     RTS=RTS_sim))
  return(output)
  
}

## another wrapper function
run_sim <- function(r,burn_in,
                    theta0,theta1,sigma,clip=0.1,
                    batch_type,
                    batch_size=NULL,batches=NULL,
                    batch_range=NULL,N_total=NULL){
  
  ## Simulate the trials
  sims <- mclapply(X=1:r, 
                   function(X){sim(burn_in=burn_in,
                                   theta0=theta0,theta1=theta1,sigma=sigma,clip=clip,
                                   batch_type=batch_type,
                                   batch_size=batch_size,batches=batches,
                                   batch_range=batch_range,N_total=N_total)},
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
  
  ## add new columns
  pvals$theta1 <- theta1
  pvals$batch_type <- batch_type
  trials$theta1 <- theta1
  trials$batch_type <- batch_type
  
  return(list(pvals=pvals,
              trials=trials))
  
}
