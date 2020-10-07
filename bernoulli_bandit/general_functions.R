## ----------------------------------------------------------------- ##
## general_functions.R --------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: store a list of common functions that are used in ------ ##
## multiple other scripts ------------------------------------------ ##
## ----------------------------------------------------------------- ##

library(MASS)
library(dplyr)

## possible treatment choices
txt_choices <- list()
txt_choices[[1]] <- c(1,0)
txt_choices[[2]] <- cbind(c(1,1),c(1,0),c(0,1),c(0,0))
txt_choices[[3]] <- cbind(c(1,1,1),
                          c(1,1,0),c(1,0,1),c(0,1,1),
                          c(1,0,0),c(0,1,0),c(0,0,1),
                          c(0,0,0))

txt_choices[[4]] <- cbind(c(1,1,1,1),
                          c(1,1,1,0),c(1,1,0,1),c(1,0,1,1),c(0,1,1,1),
                          c(1,1,0,0),c(1,0,1,0),c(1,0,0,1),c(0,1,0,1),c(0,0,1,1),c(0,1,1,0),
                          c(0,0,0,1),c(0,0,1,0),c(0,1,0,0),c(1,0,0,0),
                          c(0,0,0,0))

## tr
## Purpose: calculate the trace of a matrix
## param X: matrix
## return tr: trace of matrix X
tr <- function(X){
  tr <- sum(diag(X))
  return(tr)
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


## mean_logit
## Purpose: generate the log-odds for an individual
## param A: treatment variable in {0,1}
## param x: the covariate
## param beta0-5: coefficients that change the log-odds
## return logit: mean response for the individual
mean_logit <- function(A,x,beta0,beta1,beta2,beta3,beta4,beta5) {
  logit <- beta0 + beta1*A + 
    beta2*(x) + beta3*(x**2) + 
    beta4*A*x + beta5*A*(x**2)
  return(logit)
}

## get_prob
## Purpose: takes the log odds and turns it into probability
## param logit: the log odds
## return prob: the probability of success
get_prob <- function(logit){
  prob <- exp(logit)/(1+exp(logit))
  return(prob)
}


## gen_train
## Purpose: generate a training set to use for the burn-in period
## and future covariates
## param N: number of patients to include
## param p1: randomization probability to A=1
## param beta0-6: coefficients that change the log-odds
## x1_shape1: first shape parameter of the x1 ~ beta(shape1,shape2) cov
## x1_shape2: first shape parameter of the x1 ~ beta(shape1,shape2) cov
## return train_set: dataset with N observations includes 
##                   subject,x1,x2,A,logit,p,Y
gen_train <- function(N,p1=0.5,beta0,beta1,beta2,beta3,beta4,beta5,
                      x1_shape1=1,x1_shape2=2.5){
  
  ## generate random A treatments
  A_vec <- rbinom(N,1,p1)
  x1_vec <- rbeta(N,shape1=x1_shape1,shape2=x1_shape2)
  
  ## empty training set, filling in predetermined entires
  train_set <- matrix(NA,nrow=N,ncol=6)
  train_set[,1:3] <- cbind(1:N,x1_vec,A_vec)
  
  for(i in 1:N){
    ## gaining the mean log odds for each individual
    train_set[i,4] <- mean_logit(x=train_set[i,2],
                                 A=train_set[i,3],
                                 beta0=beta0,beta1=beta1,beta2=beta2,
                                 beta3=beta3,beta4=beta4,beta5=beta5)
    
    ## probability of success
    train_set[i,5] <- get_prob(train_set[i,4])
    
    ## success or failure
    train_set[i,6] <- rbinom(1,1,p=train_set[i,5])
  
  }
  
  train_set <- data.frame(train_set)
  colnames(train_set) <- c("subject","x1","A","log_odds","p","Y")
  return(train_set)
}

## kl_div
## Purpose: calculates the kl-divergence between two normal distributions
## param mu0: mean vector of baseline distribution
## param mu1: mean vector of new distribution
## param Sigma0: covariance matrix of baseline distribution
## param Sigma1: covariance matrix of new distribution
## param k: number of parameters
## return kl_div: kl-divergence between two dists
kl_div <- function(mu0,mu1,Sigma0,Sigma1,k){
  Sigma1Inv <- solve(Sigma1)
  kl_div <- (0.5)*( tr(Sigma1Inv%*%Sigma0) + 
                  t(mu1-mu0)%*%(Sigma1Inv%*%(mu1-mu0)) -
                  k + log(det(Sigma1)/det(Sigma0)) )
  return(kl_div)
}

## p_kl
## Purpose: calculate predicited log-odds and kl-divergence for 
## each treatment combination
## param trts: matrix of different treatment combinations
## param beta: estimated coefficient vector  
## param X0: current design matrix before new treatments
## param Sigma0: current estimated covariance matrix (X0^T W x)^{-1}
## param fit: current glm object
## return: list of trts with estimated kl-div and reward
p_kl <- function(trts,beta,new_num,X1,Sigma0,fit){
  
  trt_info <- matrix(NA,nrow=2**new_num,ncol=3)
  
  ## previous predicted values
  prev_phat <- predict(fit,type="response")
  
  for(t in 1:2**new_num){
    design <- trts[,c(1,2,3,3+t)]
    colnames(design) <- c("int","cov1","cov2","A1")
    design <- design %>%
      mutate(A1 = as.numeric(as.character(A1)),
             A1 = ifelse(A1<0,0,1),
             cov1cov2=cov1*cov2,
             cov1A1 = cov1*A1,
             cov2A1=cov2*A1)
    
    est_log_odds <- as.matrix(design) %*% as.matrix(beta)
    est_ps <- get_prob(est_log_odds)
    
    ## new design and covariance matricies
    colnames(X1) <- colnames(design)
    X1 <- na.omit(as.matrix(rbind(X1,design)))
    W1 <- diag( c(prev_phat*(1-prev_phat),est_ps) )
    Sigma1 <- solve( t(X1) %*% W1 %*% X1  )
    
    # calculate kl divergence
    kl_div <- kl(mu0=beta,mu1=beta,Sigma0=covMat,Sigma1=covMatnew,k=nrow(covMatnew))
    trt_info[t,] <- c(t,
                      ifelse(is.na(mean(p)),0,mean(p)),
                      ifelse(is.na(kl_div),0,kl_div))
  }
  
  return(trt_info)
  
}




opt_grid <- function(trts,beta,new_num){
  
  trt_info <- matrix(NA,nrow=2**new_num,ncol=2)
  
  for(t in 1:2**new_num){
    design <- trts[,c(1,2,3,3+t)]
    colnames(design) <- c("int","cov1","cov2","A1")
    design <- design %>%
      mutate(A1 = as.numeric(as.character(A1)),
             A1 = ifelse(A1<0,0,1),
             cov1cov2=cov1*cov2,
             cov1A1 = cov1*A1,
             cov2A1=cov2*A1)
    
    #colnames(X1) <- colnames(design)
    p <- as.matrix(design) %*% as.matrix(beta)
    
    trt_info[t,] <- c(t,
                      ifelse(is.na(mean(p)),0,mean(p)))
  }
  return(trt_info)
  
}

