## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: simulate a trial for each method and create plots ------ ##
## ----------------------------------------------------------------- ##

library(ggplot2)
library(gridExtra)
library(grid)

setwd("----------")
source("general_functions.R")

source("eig_simple.R")
source("det_simple.R")
source("tr_simple.R")

source("eig_ucb.R")
source("det_ucb.R")
source("tr_ucb.R")

source("eig_eps_greedy.R")
source("det_eps_greedy.R")
source("tr_eps_greedy.R")

source("eig_thompson.R")
source("det_thompson.R")
source("tr_thompson.R")

library(parallel)

## wrapper function
sim <- function(burn_in,N,
                shape1=1,shape2=1,
                beta00,beta10,beta20,beta30,
                beta01,beta11,beta21,beta31,
                beta02,beta12,beta22,beta32,
                beta03,beta13,beta23,beta33,
                beta04,beta14,beta24,beta34,
                beta05,beta15,beta25,beta35,
                sigma=0.1){
  
  
  
  ## generate training set
  train_set <- gen_train(N=N,
                         shape1=shape1,shape2=shape2,
                         beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                         beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                         beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                         beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                         beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                         beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                         sigma=sigma)
  
  
  ## EIGENVALUE AS INFORMATION GAIN METRIC
  eig_simple_sim <- try( eig_simple(train_set=train_set,burn_in=burn_in,N=N,
                            beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                            beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                            beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                            beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                            beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                            beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                            sigma=sigma) )
  
  eig_ucb_sim <- try( eig_ucb(train_set=train_set,burn_in=burn_in,N=N,
                             beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                             beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                             beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                             beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                             beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                             beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                             sigma=sigma,ci_level = 0.99) )
  
  eig_thompson_sim <- try( eig_thompson(train_set=train_set,burn_in=burn_in,N=N,
                                   beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                   beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                   beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                   beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                   beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                   beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                   sigma=sigma,B=500,c= 0.5) )
  
  eig_eps_greedy_sim <- try( eig_eps_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                          beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                          beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                          beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                          beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                          beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                          beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                          sigma=sigma,eps=0.1) )
  
  ## DET AS INFORMATION GAIN METRIC
  det_simple_sim <- try( det_simple(train_set=train_set,burn_in=burn_in,N=N,
                                    beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                    beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                    beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                    beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                    beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                    beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                    sigma=sigma) )
  
  det_ucb_sim <- try( det_ucb(train_set=train_set,burn_in=burn_in,N=N,
                              beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                              beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                              beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                              beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                              beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                              beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                              sigma=sigma,ci_level = 0.99) )
  
  det_thompson_sim <- try( det_thompson(train_set=train_set,burn_in=burn_in,N=N,
                                        beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                        beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                        beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                        beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                        beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                        beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                        sigma=sigma,B=500,c= 0.5) )
  
  det_eps_greedy_sim <- try( det_eps_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                            beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                            beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                            beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                            beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                            beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                            beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                            sigma=sigma,eps=0.1) )
  
  
  ## TR AS INFORMATION GAIN METRIC
  tr_simple_sim <- try( tr_simple(train_set=train_set,burn_in=burn_in,N=N,
                                    beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                    beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                    beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                    beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                    beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                    beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                    sigma=sigma) )
  
  tr_ucb_sim <- try( tr_ucb(train_set=train_set,burn_in=burn_in,N=N,
                              beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                              beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                              beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                              beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                              beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                              beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                              sigma=sigma,ci_level = 0.99) )
  
  tr_thompson_sim <- try( tr_thompson(train_set=train_set,burn_in=burn_in,N=N,
                                        beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                        beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                        beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                        beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                        beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                        beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                        sigma=sigma,B=500,c= 0.5) )
  
  tr_eps_greedy_sim <- try( tr_eps_greedy(train_set=train_set,burn_in=burn_in,N=N,
                                            beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
                                            beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
                                            beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
                                            beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
                                            beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
                                            beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
                                            sigma=sigma,eps=0.1) )
  
  output <- try( list(eig_simple=eig_simple_sim,
                      eig_thompson=eig_thompson_sim,
                      eig_ucb=eig_ucb_sim,
                      eig_eps_greedy=eig_eps_greedy_sim,
                      det_simple=det_simple_sim,
                      det_thompson=det_thompson_sim,
                      det_ucb=det_ucb_sim,
                      det_eps_greedy=det_eps_greedy_sim,
                      tr_simple=tr_simple_sim,
                      tr_thompson=tr_thompson_sim,
                      tr_ucb=tr_ucb_sim,
                      tr_eps_greedy=tr_eps_greedy_sim))
  return(output)
  
}


## fixed parameters
beta00=0.5;beta01=0.45;beta02=0.05;beta03=-0.5;beta04=0.03;beta05=0.01
beta10=0.35;beta11=-0.25;beta12=0.01;beta13=0.85;beta14=0.02;beta15=0.002
beta20=0.61;beta21=-0.005;beta22=0.01;beta23=0.005;beta24=0.02;beta25=0.002
beta30=0.45;beta31=-0.003;beta32=0.021;beta33=0.021;beta34=0.01;beta35=0.002

## sd param
sigma=0.1

## how many simulations to run
r <- 1

burn_in=75
N=2500

## run sim
sims <- sim(burn_in=burn_in,N=2500,
            shape1=1,shape2=1,
            beta00=beta00,beta10=beta10,beta20=beta20,beta30=beta30,
            beta01=beta01,beta11=beta11,beta21=beta21,beta31=beta31,
            beta02=beta02,beta12=beta12,beta22=beta22,beta32=beta32,
            beta03=beta03,beta13=beta13,beta23=beta23,beta33=beta33,
            beta04=beta04,beta14=beta14,beta24=beta24,beta34=beta34,
            beta05=beta05,beta15=beta15,beta25=beta25,beta35=beta35,
            sigma=sigma)

## just one replication here, so these loops are unnessary
## but if you wanted to run more, you can
good_sets <- c()
tick =1
for(i in 1:r){
  if(is.list(sims[[1]])){
    good_sets[tick]=i
    tick=tick+1
  }else{
    ## do nothin
  }
}

## organize results
trials <- data.frame()
for(i in good_sets){
  ## SIMPLE
  s_eig_simple <- sims$eig_simple[[1]]
  s_eig_simple$subject <- 1:nrow(s_eig_simple)
  s_eig_simple$method <- "eig_simple"
  s_eig_simple$info_gain <- "eig"
  s_eig_simple$method_char <- "RSR"
  s_eig_simple$rep <- i
  
  s_det_simple <- sims$det_simple[[1]]
  s_det_simple$subject <- 1:nrow(s_det_simple)
  s_det_simple$method <- "det_simple"
  s_det_simple$info_gain <- "det"
  s_det_simple$method_char <- "RSR"
  s_det_simple$rep <- i
  
  s_tr_simple <- sims$tr_simple[[1]]
  s_tr_simple$subject <- 1:nrow(s_tr_simple)
  s_tr_simple$method <- "tr_simple"
  s_tr_simple$info_gain <- "tr"
  s_tr_simple$method_char <- "RSR"
  s_tr_simple$rep <- i
  
  
  ## E GREEDY
  s_eig_eps_greedy <- sims$eig_eps_greedy[[1]]
  s_eig_eps_greedy$subject <- 1:nrow(s_eig_eps_greedy)
  s_eig_eps_greedy$method <- "eig_eps_greedy"
  s_eig_eps_greedy$info_gain <- "eig"
  s_eig_eps_greedy$method_char <- "REG(0.1)"
  s_eig_eps_greedy$rep <- i
  
  s_det_eps_greedy <- sims$det_eps_greedy[[1]]
  s_det_eps_greedy$subject <- 1:nrow(s_det_eps_greedy)
  s_det_eps_greedy$method <- "det_eps_greedy"
  s_det_eps_greedy$info_gain <- "det"
  s_det_eps_greedy$method_char <- "REG(0.1)"
  s_det_eps_greedy$rep <- i
  
  s_tr_eps_greedy <- sims$tr_eps_greedy[[1]]
  s_tr_eps_greedy$subject <- 1:nrow(s_tr_eps_greedy)
  s_tr_eps_greedy$method <- "tr_eps_greedy"
  s_tr_eps_greedy$info_gain <- "tr"
  s_tr_eps_greedy$method_char <- "REG(0.1)"
  s_tr_eps_greedy$rep <- i
  
  ## UCB
  ## UCB
  s_eig_ucb <- sims$eig_ucb[[1]]
  s_eig_ucb$subject <- 1:nrow(s_eig_ucb)
  s_eig_ucb$method <- "eig_ucb"
  s_eig_ucb$info_gain <- "eig"
  s_eig_ucb$method_char <- "RUCB(99)"
  s_eig_ucb$rep <- i
  
  s_det_ucb <- sims$det_ucb[[1]]
  s_det_ucb$subject <- 1:nrow(s_det_ucb)
  s_det_ucb$method <- "det_ucb"
  s_det_ucb$info_gain <- "det"
  s_det_ucb$method_char <- "RUCB(99)"
  s_det_ucb$rep <- i
  
  s_tr_ucb <- sims$tr_ucb[[1]]
  s_tr_ucb$subject <- 1:nrow(s_tr_ucb)
  s_tr_ucb$method <- "tr_ucb"
  s_tr_ucb$info_gain <- "tr"
  s_tr_ucb$method_char <- "RUCB(99)"
  s_tr_ucb$rep <- i
  
  ## THOMPSON SAMPLING
  ## thompson
  s_eig_thompson <- sims$eig_thompson[[1]]
  s_eig_thompson$subject <- 1:nrow(s_eig_thompson)
  s_eig_thompson$method <- "eig_thompson"
  s_eig_thompson$info_gain <- "eig"
  s_eig_thompson$method_char <- "RTH(0.5)"
  s_eig_thompson$rep <- i
  
  s_det_thompson <- sims$det_thompson[[1]]
  s_det_thompson$subject <- 1:nrow(s_det_thompson)
  s_det_thompson$method <- "det_thompson"
  s_det_thompson$info_gain <- "det"
  s_det_thompson$method_char <- "RTH(0.5)"
  s_det_thompson$rep <- i
  
  s_tr_thompson <- sims$tr_thompson[[1]]
  s_tr_thompson$subject <- 1:nrow(s_tr_thompson)
  s_tr_thompson$method <- "tr_thompson"
  s_tr_thompson$info_gain <- "tr"
  s_tr_thompson$method_char <- "RTH(0.5)"
  s_tr_thompson$rep <- i
  
  
  trials <- rbind(trials,
                  s_eig_eps_greedy,
                  s_eig_ucb,
                  s_eig_thompson,
                  s_eig_simple,
                  s_det_eps_greedy,
                  s_det_ucb,
                  s_det_thompson,
                  s_det_simple,
                  s_tr_eps_greedy,
                  s_tr_ucb,
                  s_tr_thompson,
                  s_tr_simple)
}

save(trials,file="-----")


## create eigenvalue ratio function
trials_post <- trials %>% filter(subject>burn_in) %>%
                          mutate(ratio=(log(max_eig)**1.1)/min_eig)



## create ratio plots
r1 <- ggplot(data=trials_post %>% filter(info_gain=="eig")) +
  geom_line(aes(x=subject,y=ratio,linetype=method_char)) +
  labs(x="Subject", y="Eigenvalue Ratio",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Min Eigenvalue") +
  theme_bw() +
  theme(legend.position="none")

r2 <- ggplot(data=trials_post %>% filter(info_gain=="det")) +
  geom_line(aes(x=subject,y=ratio,linetype=method_char)) +
  labs(x="Subject", y="Eigenvalue Ratio",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Determinant") +
  theme_bw() +
  theme(legend.position="none")

r3 <- ggplot(data=trials_post %>% filter(info_gain=="tr")) +
  geom_line(aes(x=subject,y=ratio,linetype=method_char)) +
  labs(x="Subject", y="Eigenvalue Ratio",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Trace") +
  theme_bw()

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

r_legend <- get_legend(r3)

r3 <- r3 + theme(legend.position="none")

grid.arrange(r2,r1,r3,r_legend,ncol=4,
             widths=c(2.3, 2.3,2.3, 0.8),
             top=textGrob("Ratio of Eigenvalues Over the Trial", 
                          gp=gpar(fontsize=17)))


## create eigenvalue growth plot
e1 <- ggplot(data=trials_post %>% filter(info_gain=="eig")) +
  geom_line(aes(x=subject,y=min_eig,linetype=method_char)) +
  labs(x="Subject", y="Minimum Eigenvalue",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Min Eigenvalue") +
  theme_bw() +
  theme(legend.position="none")

e2 <- ggplot(data=trials_post %>% filter(info_gain=="det")) +
  geom_line(aes(x=subject,y=min_eig,linetype=method_char)) +
  labs(x="Subject", y="Minimum Eigenvalue",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Determinant") +
  theme_bw() +
  theme(legend.position="none")

e3 <- ggplot(data=trials_post %>% filter(info_gain=="tr")) +
  geom_line(aes(x=subject,y=min_eig,linetype=method_char)) +
  labs(x="Subject", y="Minimum Eigenvalue",linetype="Method") +
  #ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Info Gain: Trace") +
  theme_bw()

e_legend <- get_legend(e3)

e3 <- e3 + theme(legend.position="none")

grid.arrange(e2,e1,e3,e_legend,ncol=4,
             widths=c(2.3, 2.3,2.3, 0.8),
             top=textGrob("Minimum Eigenvalues Over the Trial", 
                          gp=gpar(fontsize=17)))
