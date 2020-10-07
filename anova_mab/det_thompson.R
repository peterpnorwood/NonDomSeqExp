## ----------------------------------------------------------------- ##
## det_thompson.R -------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: simulate an experiment with ---------------------------- ##
## restricted thompson sampling ------------------------------------ ##
## ----------------------------------------------------------------- ##

setwd("-------")
source("general_functions.R")


## det_thompson
## Purpose: run a trial with simple randomization
## param train_set: training set to base first patients off of
## param burn_in: number of subjects without adaptive randomization
## param N: total number of subjects
## param mu1-mu5: mean resposnes for each group
## param sigma: standard deviation of the response
## param c: dampening constant for thompson sampling
## param br: breaks for the grid search over alpha
## return lst: list of trial information
det_thompson <- function(train_set,burn_in,N,
                         mu1,mu2,mu3,mu4,mu5,
                         sigma,B,c,br=100){
  
  
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
    
    ## gather information on reward and information gain
    det_info <- matrix(NA,nrow=5,ncol=3)
    model_matrix <- model.matrix(fit)
    for(t in 1:5){
      ## adding extra row to design matrix
      design_row <- c(rep(0,5))
      design_row[t] <- 1
      temp_matrix <- rbind(model_matrix,design_row)
      
      ## add onto det_info
      det_info[t,] <- c(t,ests[t],det(t(temp_matrix) %*% temp_matrix))
    }
    
    ## finding the treatmetns which maximize a convex combo
    
    reward_scaled <- scale(det_info[,2])
    true_reward_scaled <- scale(true_param)
    sd_info <- ifelse(sd(det_info[,3])==0,1,sd(det_info[,3]))
    info_scaled <- (det_info[,3] - mean(det_info[,3]) ) / sd_info
    
    ## 
    max_info <- comb(x=reward_scaled,y= info_scaled)
    non_dom <- comb(x=true_reward_scaled,y=info_scaled)
    
    ## generate randomization probabilities from thompson sampling
    thomp <- thompson_probs(fit=fit,txt=c(max_info),
                            B=B,c=c)
    ## save probabilites
    rand_probs[i-burn_in,] <- c(i-burn_in,thomp$probs) 
    
    
    trial[i,1] <- i
    ## new treatment for subject
    trial[i,2] <- sample(1:5,1,prob=c(thomp$probs))
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


# test <- det_thompson(train_set=train_set,burn_in=6,N=100,
#                       mu1=1,mu2=1,mu3=1,mu4=1.6,mu5=1.8,
#                       sigma=1,B=500,c=0.5)
