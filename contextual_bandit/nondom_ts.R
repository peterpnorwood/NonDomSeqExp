## ----------------------------------------------------------------- ##
## nondom_ts.R ----------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run an experiment with nondom thompson sampling -------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/NonDomSeqExp/NonDomSeqExp/contextual_bandit")
source("funcs.R")

library(MASS)

## nondom_ts
## Purpose: run an experiment with thompson sampling
## param train_set: dataset with context for N individuals
## param burn_in: sample size of simple randomization
## param A: vector of possible treatments
## param theta: true mean outcome parameter vector
## param sd_Y: standard deviation for response
## return dat: dataframe with X,A,mu,Y,regret,norm
nondom_ts <- function(train_set,burn_in,A,theta,sd_Y){
  
  ## number of subjects
  N <- nrow(train_set)
  ## dimension of context
  p <- ncol(train_set)-3
  ## number of arms
  K <- length(A)
  
  ## trial dataset
  dat <- matrix(NA,nrow=N,ncol=p+6)
  ## context
  dat[1:N,1:p] <- as.matrix(train_set)[1:N,1:p]
  ## first burn_in interventions
  dat[1:burn_in,p+1] <- train_set$A[1:burn_in]
  ## first burn_in means
  dat[1:burn_in,p+2] <- train_set$mu[1:burn_in]
  ## first burn_in outcomes
  dat[1:burn_in,p+3] <- train_set$Y[1:burn_in]
  ## name the same colnames
  colnames(dat) <- c(colnames(train_set),"regret","norm","non_dom")
  
  ## loop through the new patients
  for(i in (burn_in+1):N){
    
    ## fit the outcome model
    X_temp <- dat[1:(i-1),1:p]
    A_temp <- dat[1:(i-1),p+1]
    Y <- dat[1:(i-1),p+3]
    
    temp <- data.frame(X_temp,A=A_temp,Y)
    
    fit <- lm(Y~-1+as.factor(A)+as.factor(A):.-
                as.factor(A):A,
              data=temp)
    
    ## gather parameter convergence information
    coef_fit <- coef(fit)
    #Sigma <- vcov(fit)
    theta_hat <- c()
    ## put them in the same format as the theta vector
    tik <- 1
    for(ii in 1:K){
      for(jj in 0:p){
        theta_hat[tik] <- coef_fit[ii+(K)*jj]
        tik=tik+1
      }
    }
    
    ## measure the euclidean norm between theta and theta_hat
    dat[i,p+5] <- norm(matrix(theta-theta_hat),type="F")
    
    ## loop through interventions to find greedy intevention
    info <- matrix(NA,nrow=length(A),ncol=4)
    tick=1
    for(a in A){
      ## gather ests if a is assigned
      temp_dat <- data.frame(t(dat[i,1:p]),A=a,Y=0)
      ## estiamted mean outcome given a
      mu_hat <- predict(fit,temp_dat)
      ## true mean outcome given a
      mu <- mean_outcome(X=dat[i,1:p],A=A,a=a,theta=theta)
      ## new design
      temp_df <- rbind(temp,temp_dat)
      
      temp_X <- model.matrix(Y~-1+as.factor(A)+as.factor(A):.-
                               as.factor(A):A,temp_df)
      
      XtX <- t(temp_X) %*% temp_X
      XtXi <- solve(XtX)
      info_gain <- 1/sum(diag(XtXi))
      
      ## save info
      info[tick,] <- c(a,mu_hat,mu,info_gain)
      tick=tick+1
    }
    ## save info as dataframe
    info <- data.frame(info)
    colnames(info) <- c("A","mu_hat","mu","info_gain")
    
    ## true non-dominated
    true_nondom <- comb(info$mu,info$info_gain)
    est_nondom <- comb(info$mu_hat,info$info_gain)
    
    ## randomize via thompson sampling
    ts_probs <- thompson_probs(fit=fit,txt=est_nondom,
                               new_sub=data.frame(t(dat[i,1:p]),A=1,Y=0))
    
    ts_probs$A <- as.numeric(as.character(ts_probs$A))
    
    
    ## assign intervention (e-greedy)
    if(nrow(ts_probs)==1){
      dat[i,p+1] <- ts_probs$A[1]
    }else{
      dat[i,p+1] <- sample(ts_probs$A,1,prob=ts_probs$probs)
    }
    ## find mean outcome
    dat[i,p+2] <- info$mu[dat[i,p+1]]
    ## find outcome
    dat[i,p+3] <- rnorm(1,dat[i,p+2],sd_Y)
    ## find regret
    dat[i,p+4] <- max(info$mu) - dat[i,p+2]
    ## determine if it was non-dominated
    dat[i,p+6] <- ifelse(dat[i,p+1] %in% true_nondom,1,0)
    
  }
  
  dat <- data.frame(dat)
  dat$sub <- 1:nrow(dat)
  return(dat)
  
}
