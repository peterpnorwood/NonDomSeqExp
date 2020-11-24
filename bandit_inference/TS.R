## ----------------------------------------------------------------- ##
## TS.R ------------------------------------------------------------ ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run a batched bandit with thompson sampling ------------ ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("general_functions.R")

## TS
## Purpose: run a batched bandit with thompson sampling
## param train_set: data to start randomizing on
## param theta0: mean reward for arm 0
## param theta1: mean reward for arm 1
## param sigma: std deviation fo reward distribution
## param batch_size: size of batch before updating model
## param batches: number of batches
## clip: lower bound for randomization probs
## return lst: list with trial df and p-value df
TS <- function(train_set,theta0,theta1,batch_size,batches,sigma,clip){
  
  ## save params
  theta <- c(theta0,theta1)
  
  ## loop through the batches
  trial <- data.frame()
  for(b in 1:batches){
    ## get the temporary fit
    if(b==1){ 
      temp <- train_set
    }else{
      temp <- trial
    }
    fit <- lm(Y~-1+as.factor(A),data=temp)
    
    ## randomization probability
    ts_probs <- thompson_probs(fit=fit,B=500)
    rand_prob <- ts_probs$pi[2]
    
    ## clip the randomization probability 
    if(rand_prob>1-clip){
      rand_prob <- 1-clip
    }else if(rand_prob<clip){
      rand_prob <- clip
    }
    
    ## assign treatment
    A <- rbinom(batch_size,1,rand_prob)
    
    ## make sure they all aren't the same
    if(mean(A)==1){
      A[1] <- 0
    }else if(mean(A)==0){
      A[1] <-1
    }
    
    ## mean rewards
    mean_rewards <- mean_reward(theta=theta,A=A)
    ## observed rewards
    y <- rnorm(batch_size,mean=mean_rewards,sd=sigma)
    ## batch dataset
    temp <- data.frame(A=A,pi=rep(rand_prob,batch_size),mu=mean_rewards,Y=y,batch=b)
    trial <- rbind(trial,temp)
  }
  
  ## assign subject 
  trial$subject <- 1:(batches*batch_size)
  
  ## OLS test
  ols_p <- ols_pvalue(trial=trial)
  ## BOLS test
  bols_p <- bols_pvalue(trial=trial)
  ## pvalue dataframe
  pval_df <- data.frame(method=c("ols","bols"),
                        pval=c(ols_p,bols_p))
  
  ## return
  lst <- list(trial=trial,
              pval_df=pval_df)
  return(lst)
}
