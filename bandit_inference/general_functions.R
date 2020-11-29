## ----------------------------------------------------------------- ##
## general_functions.R --------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: store a list of common functions that are used in ------ ##
## multiple other scripts ------------------------------------------ ##
## ----------------------------------------------------------------- ##

library(MASS)
library(dplyr)

## comb
## Purpose: report which rows maximize a convex combo
## of x and y
## param x: the first variable in the combination
## param y: the second variable in the combination
## param br: number of breaks acrosse [0,1] in a grid search
## return max_combos: list of rows which return a maximum combo
comb <- function(x,y,br=100){
  ## breaks
  alphas <- seq(0,1,by=1/br)
  maxs <- matrix(NA,ncol=1,nrow=br)
  i=1
  ## grid search through alpha
  for(a in alphas){
    z <- a*x + (1-a)*y
    maxs[i] <- which.max(z)
    i=i+1
  }
  max_combos <- unique(maxs)
  return(max_combos)
}

## mean_reward
## Purpose: generate mean reward based on A
## param theta: mean vector
## param A: arm
## return mu: mean reward
mean_reward <- function(theta,A){
  mu <- ifelse(A==0,theta[1],theta[2])
  return(mu)
}

mean_reward <- Vectorize(mean_reward,vectorize.args = "A")

## gen_data
## Purpose: develop training set
## param burn_in: number of observations
## param theta0: mean reward for arm 0
## param theta1: mean reward for arm 1
## param sigma: std deviation fo reward distribution
## return trial: dataframe with A, Pr(A=1), mean & observed reward
gen_data <- function(burn_in,theta0,theta1,sigma){
  ## arms
  A <- c(0,1,rbinom(burn_in-2,1,0.5))
  ## probs
  pi <- rep(0.5,burn_in)
  ## mean rewards
  theta <- c(theta0,theta1)
  mean_rewards <- mean_reward(theta=theta,A=A)
  ## observed rewards
  y <- rnorm(burn_in,mean=mean_rewards,sd=sigma)
  
  ## trial
  trial <- data.frame(A=A,pi=pi,mu=mean_rewards,Y=y)
  return(trial)
}

## thompson_probs
## Purpose: generate posterior probabilites based on thompson sampling
## param fit: fitted linear model object
## param B: samples to take from MVN normal
## return probs: posterior probabilites
thompson_probs <- function(fit,B){
  
  ## sample from "posterior"
  mu_vec <- fit$coefficients
  Sigma <- vcov(fit)
  posterior <- mvrnorm(n=B,mu=mu_vec,Sigma=Sigma)
  
  ## check if the arm one mean is greater than arm 0
  gt <- ifelse(posterior[,2]-posterior[,1]>0,1,0)
  
  ## P(A1>A0)
  probs <- data.frame(A=c(0,1),pi=c(1-mean(gt),mean(gt)))
  return(probs)
  
}


## ols_pvalue
## Purpose: get p-value from OLS testing mu1=mu0 vs mu1>mu0
## param trial: trial data
## return p_val: pvalue from one sided test
ols_pvalue <- function(trial){
  
  ## OLS test
  ols_fit <- lm(Y~-1+as.factor(A),data=trial)
  ols_diff <- c(-1,1) %*% coef(ols_fit)
  ols_var <- c(-1,1) %*% vcov(ols_fit) %*% c(-1,1)
  ols_se <- sqrt(ols_var)
  T_stat <- ols_diff/ols_se
  p_val <- 1-pnorm(T_stat)
  
  ## return
  return(p_val)
}


## bols_pvalue
## Purpose: get p-value from Batched OLS testing mu1=mu0 vs mu1>mu0
## param trial: trial data
## return p_val: pvalue from one sided test
bols_pvalue <- function(trial){
  
  ## get batch size info
  batch_size <- nrow(trial %>% filter(batch==1))
  batches <- max(trial$batch)
  
  ## sample sizes
  N0 <- c()
  N1 <- c()
  
  ## batch level differences
  deltas <- c()
  
  ## loop through the batches
  for(b in 1:batches){
    ## batch data
    temp <- trial %>% filter(batch==b)
  
    ## save counts
    N0[b] <- sum(1-temp$A)
    N1[b] <- sum(temp$A)
    
    ## save delta
    deltas[b] <- sum(temp$A*temp$Y)/N1[b]-sum((1-temp$A)*temp$Y)/N0[b]
  }
  
  ## variance estimate
  sigma_fit <- summary(lm(Y~-1+as.factor(A),data=trial))$sigma

  ## save test stat
  inner_sum <- (N0*N1)/(batch_size*(sigma_fit**2))
  T_stat <- (1/sqrt(batches))*sum(sqrt(inner_sum)*deltas) 
  p_val <- 1-pnorm(T_stat)
  
  ## return_pval
  return(p_val)
}

## W_d_pvalue
## Purpose from W-d-correlated estimator
## param trial: trial data
## return lst: list with diff, se, p-val
# W_d_pvalue <- function(trial){
#   
#   ## design matrix
#   X <- cbind(c(1-trial$A),trial$A)
#   
#   ## initial lambda
#   lambda_0 <- 1
#   
#   ## little w
#   w <- matrix(NA,nrow=2,ncol=nrow(trial))
#   ## big W
#   W <- matrix(NA,nrow=2,ncol=nrow(trial))
#   for(i in 1:nrow(trial)){
#     if(i==1){
#       ## little
#       w[,i] <- X[i,]/(lambda_0 + 1)
#       ## big
#       W[,i] <- w[,i]
#     }else{
#       ## little
#       w[,i] <- (diag(c(1,1))- W[,(i-1)] %*% t(X[i,]))%*%X[i,]/(lambda_0+1)
#       W[,i] <- W[,1:(i-1)] %*% w[,i] 
#     }
#   }
#     
#     
# }
