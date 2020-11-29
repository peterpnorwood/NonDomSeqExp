## ----------------------------------------------------------------- ##
## RSR.R ----------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run a batched bandit with restricted ------------------- ##
## simple randomization  ------------------------------------------- ##
## ----------------------------------------------------------------- ##

setwd("-------------------------")
source("general_functions.R")

## RSR
## Purpose: run a batched bandit with RSR
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
RSR <- function(train_set,theta0,theta1,sigma,
                batch_type,
                batch_size=NULL,batches=NULL,
                batch_range=NULL,N_total=NULL){
  
  ## save params
  theta <- c(theta0,theta1)
  
  ## loop through the batches
  if(batch_type=="fixed"){
    trial <- data.frame()
    for(b in 1:batches){
      ## get the temporary fit
      if(b==1){ 
        temp <- train_set
      }else{
        temp <- trial
      }
      fit <- lm(Y~-1+as.factor(A),data=temp)
    
      ## reward and info gain matrix
      info <- matrix(NA,2,3)
      for(a in 0:1){

        ## det info
        temp_a <- matrix(0,nrow=batch_size,ncol=2)
        temp_a[,a+1] <- rep(1,batch_size)
        temp_X <- rbind(model.matrix(fit),temp_a)
        det <- det(t(temp_X) %*% temp_X)
      
        ## reward info
        reward <- coef(fit)[a+1]
      
        ## matrix of info
        info[a+1,1] <- a
        info[a+1,2] <- det
        info[a+1,3] <- reward
      }
    
      ## find convex combinations
      conx <- comb(x=info[,2],y=info[,3])-1
    
      ## assign rand probs
      if(length(conx)==2){
        rand_prob=0.5
      }else if(conx==0){
        rand_prob=0.1
      }else{
        rand_prob=0.9
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
  }else if(batch_type=="random"){
    
    ## initialize
    b <- 1; sub <- 0;
    trial <- data.frame()
    while(sub < N_total){
      ## get the temporary fit
      if(b==1){ 
        temp <- train_set
      }else{
        temp <- trial
      }
      fit <- lm(Y~-1+as.factor(A),data=temp)
      
      ## get random batch size
      batch_size2 <- sample(batch_range,1)
      
      ## reward and info gain matrix
      info <- matrix(NA,2,3)
      for(a in 0:1){
        ## det info
        temp_a <- matrix(0,nrow=batch_size2,ncol=2)
        temp_a[,a+1] <- rep(1,batch_size2)
        temp_X <- rbind(model.matrix(fit),temp_a)
        det <- det(t(temp_X) %*% temp_X)
        
        ## reward info
        reward <- coef(fit)[a+1]
        
        ## matrix of info
        info[a+1,1] <- a
        info[a+1,2] <- det
        info[a+1,3] <- reward
      }
      
      ## find convex combinations
      conx <- comb(x=info[,2],y=info[,3])-1
      
      ## assign rand probs
      if(length(conx)==2){
        rand_prob=0.5
      }else if(conx==0){
        rand_prob=0.1
      }else{
        rand_prob=0.9
      }
      ## assign treatment
      A <- rbinom(batch_size2,1,rand_prob)
      ## make sure they all aren't the same
      if(mean(A)==1){
        A[1] <- 0
      }else if(mean(A)==0){
        A[1] <-1
      }
      ## mean rewards
      mean_rewards <- mean_reward(theta=theta,A=A)
      ## observed rewards
      y <- rnorm(batch_size2,mean=mean_rewards,sd=sigma)
      ## batch dataset
      temp <- data.frame(A=A,pi=rep(rand_prob,batch_size2),mu=mean_rewards,Y=y,batch=b)
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

# test_RSR <- RSR(train_set=train_set,theta0=0,theta1=0,
#                 sigma=1,batch_type="random",batch_range=3:20,N_total=1000)
# 
# test_RSR <- RSR(train_set=train_set,theta0=0,theta1=0,
#                sigma=1,batch_type="fixed",batch_size=20,batches=50)
