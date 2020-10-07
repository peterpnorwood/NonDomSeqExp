## ----------------------------------------------------------------- ##
## eig_simple.R ---------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: simulate an experiment with restricted simple ---------- ##
## randomization where the min eigenvalue is the info gain metric -- ##
## ----------------------------------------------------------------- ##

setwd("-----")
source("general_functions.R")


## eig_simple
## Purpose: run a trial with simple randomization
## param train_set: training set to base first patients off of
## param burn_in: number of subjects without adaptive randomization
## param N: total number of subjects
## param beta00-beta35: mean parameters
## param sigma: standard deviation of the response
## return lst: list with trial information
eig_simple <- function(train_set,burn_in,N,
                       beta00,beta10,beta20,beta30,
                       beta01,beta11,beta21,beta31,
                       beta02,beta12,beta22,beta32,
                       beta03,beta13,beta23,beta33,
                       beta04,beta14,beta24,beta34,
                       beta05,beta15,beta25,beta35,
                       sigma=0.1){
  
  
  ## data we are interested in storing
  trial <- matrix(NA,nrow=N,ncol=6)
  param_ests <- matrix(NA,nrow=N-burn_in,ncol=25)
  std_errors <- matrix(NA,nrow=N-burn_in,ncol=25)
  
  trial <- data.frame(trial)
  colnames(trial) <- colnames(train_set)
  trial <- train_set
  trial$non_dom <- NA
  trial$opt <- NA
  trial$min_eig <- NA
  trial$max_eig <- NA
  
  for(i in (burn_in+1):N){
  
    
    ## fit the reg model
    fit <- lm(Y ~ X1 + X2 + X1:X2 + I(X1^2) + I(X1^2):X2 + 
                as.factor(A) + as.factor(A):(X1 + X2 + X1:X2 + I(X1^2) + I(X1^2):X2),data=trial[1:i-1,])
    
    ## save estimates & standard errors
    ests <- fit$coefficients
    param_ests[i-burn_in,] <- c(i-burn_in,ests)
    ses <- diag(vcov(fit))
    std_errors[i-burn_in,] <- c(i-burn_in,ses)
    
    design <- model.matrix(fit)
    eigs_design <- eigen( t(design) %*% design)
    
    ## gather information on reward and information gain
    eig_info <- matrix(NA,nrow=4,ncol=4)
    for(a in 0:3){
      
      add <- data.frame(trial[i,1:3],a,1,1,1,1,1,1)
      colnames(add) <- colnames(trial)
      temp <- rbind(trial[1:i-1,],add)
      temp_matrix <- model.matrix(fit,data=temp)
      
      ## gather the eigenvalues
      XtX <- t(temp_matrix) %*% temp_matrix
      eigs <- eigen(XtX)
      
      ## gather the predicted response
      yhat <- predict(fit,add)
      mu <- mean_response(X1=trial[i,2],X2=trial[i,3],A=a,
                          beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                          beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                          beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                          beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                          beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                          beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35)
      
      ## add onto det_info
      eig_info[a+1,] <- c(a,yhat,mu,min(eigs$values))
      
    }
    
    ## finding the treatmetns which maximize a convex combo
    ## minus one to index these correctly
    max_info <- comb(x=eig_info[,2],y=eig_info[,4])-1
    true_non_dom <- comb(x=eig_info[,3],y=eig_info[,4])-1
    
    eig_info <- as.data.frame(eig_info)
    colnames(eig_info) <- c("A","yhat","mu","info_gain")
    eig_info$rand_prob <- ifelse(eig_info$A %in% max_info,1,0)/length(max_info)
    
    ## save patient information
    trial[i,1] <- i
    ## new treatment for subject
    trial[i,4] <- sample(0:3,1,prob=eig_info$rand_prob)
    ## mean response for subject
    trial[i,5] <- mean_response(X1=trial[i,2],
                                X2=trial[i,3],
                                A=trial[i,4],
                                beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35)
    ## response for subject
    trial[i,6] <- rnorm(1,mean=trial[i,5],sd=sigma)
    ## indicator if this is a truly non-dominated treatment
    trial[i,7] <- ifelse(trial[i,4] %in% true_non_dom,1,0)
    ## was the optimal treatment given
    trial[i,8] <- ifelse(trial[i,4]==eig_info$A[which.max(eig_info$mu)],1,0)
    trial[i,9] <- min(eigs_design$values)
    trial[i,10] <- max(eigs_design$values)
  }
  
  lst <- list(trial=trial,
              param_ests=param_ests,
              std_errors=std_errors)
  
  return(lst)
  
}


# beta00=0.5;beta01=0.45;beta02=0.05;beta03=-0.5;beta04=0.03;beta05=0.01
# beta10=0.35;beta11=-0.25;beta12=0.01;beta13=0.85;beta14=0.02;beta15=0.002
# beta20=0.61;beta21=-0.005;beta22=0.01;beta23=0.005;beta24=0.02;beta25=0.002
# beta30=0.45;beta31=-0.003;beta32=0.021;beta33=0.021;beta34=0.01;beta35=0.002
# 
# train_set <- gen_train(N=100,
#                        shape1=0.25,shape2=0.75,
#                        pX1=0.5,pA=0.5,
#                        beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
#                        beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
#                        beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
#                        beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
#                        beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
#                        beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
#                        sigma=0.1)
# 
# test <- eig_simple(train_set,burn_in=20,N=nrow(train_set),
#                beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
#                beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
#                beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
#                beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
#                beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
#                beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
#                sigma=0.1)
