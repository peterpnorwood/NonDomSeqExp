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


## mean_response
## Purpose: generate the mean response for an individual
## param X1: covariate one
## param X2: covariate two
## param A: treatment variable in {1,2,3,4,5}
## param beta00-beta35: coefficients for linear model
## return mu: mean response for the individual
mean_response <- function(X1,X2,A,
                          beta00,beta10,beta20,beta30,
                          beta01,beta11,beta21,beta31,
                          beta02,beta12,beta22,beta32,
                          beta03,beta13,beta23,beta33,
                          beta04,beta14,beta24,beta34,
                          beta05,beta15,beta25,beta35) {
  
  A0=ifelse(A==0,1,0)
  A1=ifelse(A==1,1,0)
  A2=ifelse(A==2,1,0)
  A3=ifelse(A==3,1,0)
  
  mu <- beta00*A0 + beta01*A0*X1 + beta02*A0*X2 +
        beta03*A0*(X1**2) + beta04*A0*(X1*X2) +
        beta05*A0*X2*(X1**2) +
    
    
        beta10*A1 + beta11*A1*X1 + beta12*A1*X2 +
        beta13*A1*(X1**2) + beta14*A1*(X1*X2) +
        beta15*A1*X2*(X1**2) +
        beta20*A2 + beta21*A2*X1 + beta22*A2*X2 +
        beta23*A2*(X1**2) + beta24*A2*(X1*X2) +
        beta25*A2*X2*(X1**2) +
        beta30*A3 + beta31*A3*X1 + beta32*A3*X2 +
        beta33*A3*(X1**2) + beta34*A3*(X1*X2) +
        beta35*A3*X2*(X1**2)
    
    
    
  return(mu)
}

## gen_train
## Purpose: generate a training set to use for the burn-in period
## and future covariates
## param N: number of patients to include
## param shape1: first parameter for distribution of x2
## param shape2: second parameter for distribution of x2
## param pX1: P(X1=1)
## param pA: P(A=1) before adaptive randomization
## param beta00-beta35: coefficients for linear model
## param sigma: common standard deviation for response
## return train_set: dataset with N observations includes subject, A, mu, Y
gen_train <- function(N,
                      shape1,shape2,
                      pX1=0.5,pA=0.5,
                      beta00,beta10,beta20,beta30,
                      beta01,beta11,beta21,beta31,
                      beta02,beta12,beta22,beta32,
                      beta03,beta13,beta23,beta33,
                      beta04,beta14,beta24,beta34,
                      beta05,beta15,beta25,beta35,
                      sigma=0.1){
  
  ## generate random A treatments
  txt <- c(0,1)
  X1 <- rbeta(N,shape1,shape2)
  X2 <- rbinom(N,1,pX1)
  A <- sample(0:3,N,replace=TRUE)
  subs <- cbind(1:N)
  
  ## empty training set, filling in predetermined entires
  train_set <- matrix(NA,nrow=N,ncol=6)
  train_set[,1:4] <- c(subs,X1,X2,A)
  
  for(i in 1:N){
    ## gaining the mean response at stage one
    train_set[i,5] <- mean_response(X1=train_set[i,2],
                                    X2=train_set[i,3],
                                    A=train_set[i,4],
                                    beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                    beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                    beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                    beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                    beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                    beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35)
    
    ## determining the success at stage one based on the probability
    train_set[i,6] <- rnorm(1,mean=train_set[i,5],sd=sigma)
  }
  
  train_set <- data.frame(train_set)
  colnames(train_set) <- c("subject","X1","X2","A","mu","Y")
  return(train_set)
  
}

tr <- function(matrix){
  tr <- sum(diag(matrix))
  return(tr)
}



## thompson_probs
## Purpose: generate randomization probabilities via thompson sampling
## param fit: regression object
## param txt: vector of treatments
## param new_sub: new subject to evaluate
## param B: MC reps to estimate posterior probabilites
## param c: dampening probability
## return probs: randomization probabilites for each trt
thompson_probs <- function(fit,txt,new_sub,B=500,c=0.5){
  
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
  info <- info %>% 
    #mutate(n=table(factor(maxs,levels=1:5))/B ) %>%
    #group_by(maxs) %>% 
    #count() %>% 
    mutate(probs_raw = n/B,
           ## dampen the probabilites based on the constant
           probs_c = probs_raw**c)
  
  ## final probability dataframe to return
  probs <- data.frame(A=txt,
                      probs=NA)
  info$probs <- info$probs_c / sum(info$probs_c)
  probs <- info[,c(1,5)]
  
  return(probs)
  
}
