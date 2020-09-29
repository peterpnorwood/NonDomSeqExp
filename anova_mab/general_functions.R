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
## param A: treatment variable in {1,2,3,4,5}
## param mu1-mu5: mean response for each group 1-5
## return mu: mean response for the individual
mean_response <- function(A,mu1,mu2,mu3,mu4,mu5) {
  mu <- mu1*ifelse(A==1,1,0) + mu2*ifelse(A==2,1,0) + mu3*ifelse(A==3,1,0) + 
        mu4*ifelse(A==4,1,0) + mu5*ifelse(A==5,1,0)
  return(mu)
}


## gen_train
## Purpose: generate a training set to use for the burn-in period
## and future covariates
## param N: number of patients to include
## param p1-p5: randomization probabilities for each group
## param mu1-mu5: mean response for each group 1-5
## param sigma: common standard deviation for each group
## return train_set: dataset with N observations includes subject, A, mu, Y
gen_train <- function(N,
                      p1=0.2,p2=0.2,p3=0.2,p4=0.2,p5=0.2,
                      mu1=1,mu2=1,mu3=1,mu4=1,mu5=1,
                      sigma=1){
  
  ## generate random A treatments
  txt <- 1:5
  A <- c(1:5,sample(txt,N-5,prob=c(p1,p2,p3,p4,p5),replace=TRUE))
  subs <- cbind(1:N)
  
  ## empty training set, filling in predetermined entires
  train_set <- matrix(NA,nrow=N,ncol=4)
  train_set[,1:2] <- c(subs,A)
  
  for(i in 1:N){
    ## gaining the mean response at stage one
    train_set[i,3] <- mean_response(A=train_set[i,2],
                                    mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5)
    
    ## determining the success at stage one based on the probability
    train_set[i,4] <- rnorm(1,mean=train_set[i,3],sd=sigma)
  }
  
  train_set <- data.frame(train_set)
  colnames(train_set) <- c("subject","A","mu","Y")
  return(train_set)
  
}


## thompson_probs
## ## Purpose: generate randomization probabilities via thompson sampling
## param fit: lm fit based on the current data
## param txt: treatments to consider
## param B: samples to take from mvn distribution
## c: dampening constant for the randomization probabilities
## return train_set: dataset with N observations includes subject, A, mu, Y
thompson_probs <- function(fit,txt,B=500,c=0.5){
  
  ## generate B sample from a MVN distribution
  cov_mat <- vcov(fit)
  mu_vec <- fit$coefficients
  samples <- mvrnorm(n=B,mu=mu_vec,Sigma=cov_mat)
  samples_txt <- samples[,c(txt)]
  ## loop through the samples
  maxs <- c()
  if(length(txt)==1){
    for(s in 1:B){
      maxs[s] <- txt[which.max(samples_txt[s])]
    }
  }else{
    for(s in 1:B){
      maxs[s] <- txt[which.max(samples_txt[s,])]
    }
  }
  
  ## calculate probabilities of being superior
  #info <- data.frame(x=rep(1,B),maxs)
  info <- data.frame(table(factor(maxs, levels = 1:5)))
  colnames(info) <- c("A","n")
  info <- info %>% 
          #mutate(n=table(factor(maxs,levels=1:5))/B ) %>%
          #group_by(maxs) %>% 
          #count() %>% 
          mutate(probs_raw = n/B,
                 ## dampen the probabilites based on the constant
                 probs_c = probs_raw**c)
  
  ## final probability dataframe to return
  probs <- data.frame(A=1:5,
                      probs=NA)
  info$probs <- info$probs_c / sum(info$probs_c)
  probs <- info[,c(1,5)]
  
  return(probs)
  
}



