## ----------------------------------------------------------------- ##
## funcs.R --------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: general functions used for the study ------------------- ##
## ----------------------------------------------------------------- ##

library(MASS)
library(dplyr)

## make_design
## Purpose: reformat design matrix to align with lm
## param K: number of arms
## param p: dim(X)
## param a: treatment
## param X: context 
## return design: vector alining with lm
make_design <- function(K,p,a,X){
  
  ## create design
  design <- c()
  for(j in 1:K){
    if(a==j){
      x <- c(1,X)
    }else{
      x <- rep(0,p+1)
    }
    design <- c(design,x)
  }
  
  return(design)
  
}

## mean_outcome1
## Purpose: generate mean outcome
## param X: p-dimensional context
## param a: intervention
## param A: possible interventions
## param theta: mean parameter vector
## return mu: mean response
mean_outcome1 <- function(X,a,A,theta){
  
  ## dimensions
  K <- length(A)
  p <- length(X)
  
  design <- make_design(K=K,p=p,a=a,X=X)
  
  mu <- design %*% theta
  
  ## return the mean outcome
  return(mu)
  
}

## "vectorize" mean outcome for X and A
mean_outcome <- function(X,a,A,theta){
  N <- max(1,nrow(X))
  
  if(N==1){
    mu <- mean_outcome1(X=X,a=a,A=A,theta=theta)
  } else{
    mu <- c()
    ## loop through rows of x
    for(n in 1:N){
      mu[n] <- mean_outcome1(X=X[n,],a=a[n],A=A,theta=theta)
    }
  }
  return(mu)
}

## gen_data
## Purpose: generate datasets
## param N: sample size
## param p: dimension of the context
## param A: vector of possible interventions
## param theta: mean parameter vector
## param sd_X: standard deviation of random error
## return dat: dataset with columns: X, A, mu, Y
gen_data <- function(N,p,sd_X,A,sd_Y,theta){
  
  ## generate context
  X <- 0.5*mvrnorm(N,rep(0,p),sd_X*diag(p))

  ## create randomly assigned A
  A_vec <- sample(A,N,replace=TRUE)
  
  ## mean outcome
  mu <- mean_outcome(X=X,a=A_vec,A=A,theta=theta)
  ## true outcome
  Y <- rnorm(N,mu,sd_Y)
  
  ## dataset to return
  dat <- data.frame(X,A=A_vec,mu,Y)
  
  return(dat)
  
}

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

## thompson_probs
## Purpose: generate randomization probabilities via thompson sampling
## param fit: regression object
## param txt: vector of treatments
## param new_sub: new subject to evaluate
## param B: MC reps to estimate posterior probabilites
## param c: dampening probability
## return probs: randomization probabilites for each trt
thompson_probs <- function(fit,txt,new_sub,B=100,c=1){
  
  ## generate B sample from a MVN distribution
  cov_mat <- vcov(fit)
  mu_vec <- fit$coefficients
  samples <- mvrnorm(n=B,mu=mu_vec,Sigma=cov_mat)
  samples_txt <- samples[,c(txt)]
  
  ## loop through the samples
  maxs <- c()
  if(length(txt)==1){
    maxs <- txt
  }else{
    for(s in 1:B){
      a_info <- matrix(NA,nrow=length(txt),ncol=2)
      for(t in 1:length(txt)){
        new_sub$A <- txt[t]
        new_design <- matrix(model.matrix(fit,data=new_sub))
        p <- t(new_design) %*% samples[s,]
        a_info[t,1] <- txt[t]
        a_info[t,2] <- p
      }
      maxs[s] <- a_info[which.max(a_info[,2]),1]
    }
  }
  
  ## calculate probabilities of being superior
  #info <- data.frame(x=rep(1,B),maxs)
  info <- data.frame(table(factor(maxs, levels = txt)))
  colnames(info) <- c("A","n")
  info$probs <- info$n / B 
  probs <- info[,c(1,3)]
  
  return(probs)
  
}

#test <- gen_data(N=1000,p=5,sd_X=0.5,A=1:5,sd_Y=1,theta=mvrnorm(n=1,rep(0,25),diag(25)))

