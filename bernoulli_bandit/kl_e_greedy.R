## ----------------------------------------------------------------- ##
## e-greedy.R ------------------------------------------------------ ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: simulate an experiment with restricted ----------------- ##
## e-greedy randomization ------------------------------------------ ##
## ----------------------------------------------------------------- ##

setwd("-------")
source("general_functions.R")

library(dplyr)

## e_greedy
## Purpose: run a trial with e-greedy randomization
## param train_set: training set to base first patients off of
## param burn_in: number of subjects without adaptive randomization
## param N: total number of subjects
## param beta0-beta5: coefficients for the linear combo modeling the log odds
## param eps: randomize to greedy decision with 1-eps
## return lst: list three objects: the trial information, the param estimates
##             and the standard errors for the param estimates
kl_e_greedy <- function(train_set,burn_in,N,
                        beta0,beta1,beta2,beta3,beta4,beta5,eps){
  
  
  ## data we are interested in storing
  trial <- matrix(NA,nrow=N,ncol=10)
  true_param <- c(beta0,beta1,beta2,beta3,beta4,beta5)
  param_ests <- matrix(NA,nrow=N-burn_in,ncol=7)
  std_errors <- matrix(NA,nrow=N-burn_in,ncol=7)
  
  trial <- data.frame(trial)
  colnames(trial) <- c("subject","x1","A","log_odds","p","Y",
                       "non_dom","opt_A","regret_log_odds","regret_p")
  trial[1:burn_in,1:6] <- train_set[1:burn_in,]
  
  num_subjects <- burn_in
  i=1
  while(num_subjects<N){
    
    ## fit the anova model
    fit <- glm(Y~A+x1+I(x1**2)+x1:A+I(x1**2):A,data=na.omit(trial[,1:6]),
               family=binomial(link="logit"))
    
    X0 <- model.matrix(fit)
    Sigma0_est <- vcov(fit)
    predicted_vec <- predict(fit,type="response")
    p_vec <- na.omit(trial$p)
    Sigma0_true <- solve(t(X0) %*% diag(p_vec*(1-p_vec)) %*% X0)
    
    ## save estimates & standard errors
    ests <- fit$coefficients
    param_ests[i,] <- c(i,ests)
    ses <- diag(Sigma0_est)
    std_errors[i,] <- c(i,ses)
    
    ## bring in new subjects
    num_new <- sample(c(2,3,4),1)
    txt <- txt_choices[[num_new]]
    new_subjects <- train_set[1+num_subjects:(num_subjects+num_new-1),1:3]
    
    ## gather information on  information gain
    kl_info <- matrix(NA,nrow=2**num_new,ncol=5)
    for(t in 1:2**num_new){
      
      ## adding extra rows to design matrix
      new_subjects$A <- txt[,t]
      new_subjects$Y <- rnorm(num_new)
      X_add <- model.matrix(fit,data=new_subjects)
      X1 <- rbind(X0,X_add)
      
      ## estimated probabilites and cov mat
      log_hat <- X_add %*% ests
      p_hat <- get_prob(log_hat)
      p_hat_vec <- c(predicted_vec,p_hat)
      W_hat <- diag(p_hat_vec*(1-p_hat_vec))
      Sigma1_est <- solve( t(X1) %*% W_hat %*% X1)
      
      ## true probabilites
      log_true <- X_add %*% true_param
      p_true <- get_prob(log_true)
      p_true_vec <- c(p_vec,p_true)
      W_true <- diag(p_true_vec*(1-p_hat_vec))
      Sigma1_true <- solve( t(X1) %*% W_true %*% X1)
      
      ## add onto det_info
      kl_hat <- kl_div(mu0=ests,mu1=ests,
                       Sigma0=Sigma0_est,Sigma1=Sigma1_est,k=ncol(Sigma1_est))
      kl_true <- kl_div(mu0=true_param,mu1=true_param,
                        Sigma0=Sigma0_true,Sigma1=Sigma1_true,k=ncol(Sigma1_true)) 
      
      
      kl_info[t,] <- c(t,
                       ifelse(!is.na(mean(p_hat)),mean(p_hat),0),
                       ifelse(!is.na(mean(p_true)),mean(p_true),0),
                       ifelse(!is.na(kl_hat),kl_hat,0),
                       ifelse(!is.na(kl_true),kl_true,0))
      
    }
    
    ## finding the treatments which maximize a convex combo
    est_reward_scaled <- scale(kl_info[,2])
    true_reward_scaled <- scale(kl_info[,3])
    
    est_info_scaled <- scale(kl_info[,4])
    true_info_scaled <- scale(kl_info[,5])
    
    ## true convex combos
    est_non_dom <- comb(x=est_reward_scaled,y=est_info_scaled)
    true_non_dom <- comb(x=true_reward_scaled,y=true_info_scaled)

    
    ## choose treatment combo randomly
    num_options <- 2**num_new
    prob_vec <- rep(0,num_options)
    
    if(length(est_non_dom)==1){
      prob_vec[est_non_dom] = 1
    } else{
      prob_vec[c(est_non_dom)] <- eps/(length(est_non_dom)-1)
      prob_vec[which.max(kl_info[,2])] <- 1-eps
    }
    
    
    combo <- sample(num_options,1,prob=prob_vec)
    non_dom <- ifelse(combo %in% true_non_dom,1,0)
    
    ## additional columns to add
    new_subjects$A <- txt[,combo]
    new_subjects$log_odds <- NA
    new_subjects$p <- NA
    new_subjects$non_dom <- non_dom
    new_subjects$opt_A <- NA
    new_subjects$regret_log_odds <- NA
    new_subjects$regret_p <- NA
    
    for(j in 1:num_new){
      
      pat_info <- matrix(NA,nrow=2,ncol=3)
      for(A in c(0,1)){
        pat_info[A+1,1] <- A
        pat_info[A+1,2] <- mean_logit(A=A,x=new_subjects$x1[j],
                                      beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,
                                      beta4=beta4,beta5=beta5)
        #if(is.na(pat_info[A+1,2])){print(new_subjects[j,])}
        pat_info[A+1,3] <- get_prob(pat_info[A+1,2])
      }
      
      #print("pat_info simple")
      #print(pat_info[,2])
      #print(pat_info[,3])
      #print(which.max(pat_info[,3]))
      
      
      new_subjects$log_odds[j] <- pat_info[new_subjects$A[j]+1,2]
      new_subjects$p[j] <- pat_info[new_subjects$A[j]+1,3]
      new_subjects$Y[j] <- rbinom(1,1,new_subjects$p[j])
      new_subjects$opt_A[j] <- pat_info[which.max(pat_info[,3]),1]
      new_subjects$regret_log_odds[j] <- abs( max(pat_info[,2])-new_subjects$log_odds[j] )
      new_subjects$regret_p[j] <- abs( max(pat_info[,3])-new_subjects$p[j] )
      
    }
    
    ## save trial information
    new_subjects <- new_subjects %>% dplyr::select(subject,x1,A,log_odds,p,Y,
                                                   non_dom,opt_A,regret_log_odds,regret_p)
    trial[1+num_subjects:(num_subjects+num_new-1),] <- new_subjects
    num_subjects <- num_subjects+num_new
  }
  
  lst <- list(trial=trial,
              param_ests=param_ests,
              std_errors=std_errors)
  
  return(lst)
  
}


#test <- kl_e_greedy(train_set=train_set,burn_in=20,N=95,
#               beta0=0.1,beta1=0.1,beta2=0.1,beta3=0.1,beta4=0.1,beta5=0.1,
#               eps=0.25)  
