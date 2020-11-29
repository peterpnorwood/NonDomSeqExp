## ----------------------------------------------------------------- ##
## SR.R ------------------------------------------------------------ ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run a batched "bandit" with simple randomization  ------ ##
## ----------------------------------------------------------------- ##

setwd("-------------------------")
source("general_functions.R")

## SR
## Purpose: run a batched "bandit" with simple randomization
## param train_set: data to start randomizing on
## param theta0: mean reward for arm 0
## param theta1: mean reward for arm 1
## param sigma: std deviation fo reward distribution
## param batch_type: "random" or "fixed" to denote type of batches
## param batch_size: size of batch before updating model
## param batches: number of batches
## param batch_range: range of random batches
## Parma N_total: number of total subjects if using random batches
## return lst: list with trial df and p-value df
SR <- function(train_set,theta0,theta1,sigma,
               batch_type,
               batch_size=NULL,batches=NULL,
               batch_range=NULL,N_total=NULL){
  
  ## save params
  theta <- c(theta0,theta1)
  
  ## first fit
  fit <- lm(Y~-1+as.factor(A),data=train_set)
  
  ## run different types of trials based on fixed or random
  if(batch_type=="fixed"){
    
    trial <- data.frame()
    for(b in 1:batches){
      ## assign treatment
      A <- rbinom(batch_size,1,0.5)
      ## make sure they're not all the same for BOLS
      if(mean(A)==1){
        A[1] <- 0
      }else if(mean(A)==0){
        A[1] <- 1
      }
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
    
  }else if(batch_type=="random"){
    
    ## initialize
    b <- 1; sub <- 0;
    ## first random batches
    trial <- data.frame()
    while(sub < N_total){
      ## get random batch_size
      batch_size2 <- sample(batch_range,1)
      ## assign treatment
      A <- rbinom(batch_size2,1,0.5)
      ## make sure they're not all the same for BOLS
      if(mean(A)==1){
        A[1] <- 0
      }else if(mean(A)==0){
        A[1] <- 1
      }
      ## randomization probabilites 
      pi <- rep(0.5,batch_size2)
      ## mean rewards
      mean_rewards <- mean_reward(theta=theta,A=A)
      ## observed rewards
      y <- rnorm(batch_size2,mean=mean_rewards,sd=sigma)
      ## batch dataset
      temp <- data.frame(A=A,pi=pi,mu=mean_rewards,Y=y,batch=b)
      trial <- rbind(trial,temp)
      ## update
      b <- b+1
      sub <- sub+batch_size2
    }
    
  }
  
  ## assign subject 
  trial$subject <- 1:nrow(trial)

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

# test_SR <- SR(train_set=train_set,theta0=0,theta1=0,
#                sigma=1,batch_type="random",batch_=20,
#                batches=50)
