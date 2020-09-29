## ----------------------------------------------------------------- ##
## simple.R -------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: simulate an experiment with simple randomization ------- ##
## ----------------------------------------------------------------- ##

#setwd("~/Research/Non_Dominated/anova_mab")
source("general_functions.R")


## simple
## Purpose: run a trial with simple randomization
## param train_set: training set to base first patients off of
## param burn_in: number of subjects without adaptive randomization
## param N: total number of subjects
## param p1-p5: randomization probabilities
## param mu1-mu5: mean resposnes for each group
## sigma: standard deviation of the response
## return lst: list of trial information
simple <- function(train_set,burn_in,N,
                   p1=0.2,p2=0.2,p3=0.2,p4=0.2,p5=0.2,
                   mu1,mu2,mu3,mu4,mu5,
                   sigma){
  
  
  ## data we are interested in storing
  trial <- matrix(NA,nrow=N,ncol=5)
  true_param <- c(mu1,mu2,mu3,mu4,mu5)
  param_ests <- matrix(NA,nrow=N-burn_in,ncol=6)
  std_errors <- matrix(NA,nrow=N-burn_in,ncol=6)
  rand_probs <- matrix(NA,nrow=N-burn_in,ncol=6)
  
  trial <- data.frame(trial)
  colnames(trial) <- c("subject","A","mu","Y","non_dom")
  trial[1:burn_in,1:4] <- train_set[1:burn_in,]
  
  for(i in (burn_in+1):N){
    
    ## fit the anova model
    fit <- lm(Y~as.factor(A)-1,data=na.omit(trial[,1:4]))
    
    ## save estimates & standard errors
    ests <- fit$coefficients
    param_ests[i-burn_in,] <- c(i-burn_in,ests)
    ses <- diag(vcov(fit))
    std_errors[i-burn_in,] <- c(i-burn_in,ses)
    
    ## save probabilites
    rand_probs[i-burn_in,] <- c(i-burn_in,p1,p2,p3,p4,p5) 
    
    ## gather information on information gain
    det_info <- matrix(NA,nrow=5,ncol=3)
    #true_det_info <- matrix(NA,nrow=5,ncol=3)
    model_matrix <- model.matrix(fit)
    for(t in 1:5){
      
      ## adding extra row to design matrix
      design_row <- c(rep(0,5))
      design_row[t] <- 1
      temp_matrix <- rbind(model_matrix,design_row)
      
      ## add onto det_info
      det_t <- det(t(temp_matrix) %*% temp_matrix)
      det_info[t,] <- c(t,true_param[t],det_t)
      
    }
    
    ## finding the treatments which maximize a convex combo
    true_reward_scaled <- scale(true_param)
    
    sd_info <- ifelse(sd(det_info[,3])==0,1,sd(det_info[,3]))
    info_scaled <- (det_info[,3] - mean(det_info[,3]) ) / sd_info
    
    ## true convex combos
    non_dom <- comb(x=true_reward_scaled,y=info_scaled)
    
    ## save trial information
    trial[i,1] <- i
    ## new treatment for subject
    trial[i,2] <- sample(1:5,1,prob=c(p1,p2,p3,p4,p5))
    ## mean response for subject
    trial[i,3] <- mean_response(A=trial[i,2],
                                  mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,mu5=mu5)
    ## response for subject
    trial[i,4] <- rnorm(1,mean=trial[i,3],sd=sigma)
    
    ## was the choice non_dominated
    trial[i,5] <- ifelse(trial[i,2] %in% non_dom,1,0)
  }
  
  lst <- list(trial=trial,
              param_ests=param_ests,
              std_errors=std_errors,
              rand_probs=rand_probs)
  
  return(lst)
  
}


# test <- simple(train_set=train_set,burn_in=6,N=100,
#                       mu1=1,mu2=1,mu3=1,mu4=1.6,mu5=1.8,
#                       sigma=1)#B=500,c=0.5)
