## ----------------------------------------------------------------- ##
## SR.R ------------------------------------------------------------ ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run a batched "bandit" with simple randomization  ------ ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("general_functions.R")

## SR
## Purpose: run a batched "bandit" with simple randomization
## param train_set: data to start randomizing on
## param theta0: mean reward for arm 0
## param theta1: mean reward for arm 1
## param sigma: std deviation fo reward distribution
## param batch_size: size of batch before updating model
## param batches: number of batches
## return lst: list with trial df and p-value df
SR <- function(train_set,theta0,theta1,batch_size,batches,sigma){
  
  ## save params
  theta <- c(theta0,theta1)
  
  ## first fit
  fit <- lm(Y~-1+as.factor(A),data=train_set)
  trial <- data.frame()
  for(b in 1:batches){
    ## assign treatment
    A <- rbinom(batch_size,1,0.5)
    ## randomization probabilites 
    pi <- rep(0.5,batch_size)
    ## mean rewards
    mean_rewards <- mean_reward(theta=theta,A=A)
    ## observed rewards
    y <- rnorm(batch_size,mean=mean_rewards,sd=sigma)
    ## batch dataset
    temp <- data.frame(A=A,pi=pi,mu=mean_rewards,Y=y,batch=b)
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

