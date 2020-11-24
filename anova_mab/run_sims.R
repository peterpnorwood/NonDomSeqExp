## ----------------------------------------------------------------- ##
## run_sim.R ------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the two arm bandit ----------------- ##
## ----------------------------------------------------------------- ##

setwd("~/Research/NonDomSeqExp/NonDomSeqExp/bandit_inference")
source("sims.R")

library(parallel)
library(tidyverse)

## parameters 
burn_in=5
theta0=0
sigma=1
batch_size=100
batches=10
eps=0.1
clip=0.1

## reps
r=10

trials <- data.frame()
pvals <- data.frame()
  
wrapper <- function(theta1){
  
  run <- run_sim(burn_in=burn_in,theta0=theta0,
                 theta1=theta1,sigma=sigma,
                 batch_size=batch_size,batches=batches, 
                 eps=eps, clip=clip, r=r)
  
  lst <- list(trials=run$trials,pvals=run$pvals)
}

theta_list <- seq(0,0.5,by=0.25)

results <- mclapply(theta_list,wrapper,
                    mc.cores=1)


# pvals_null <- run_sim(burn_in=burn_in,theta0=theta0,
#                       theta1=0,sigma=sigma,
#                       batch_size=250,batches=10, 
#                       eps=eps, clip=clip,r=100)
# 
# pvals_null_2 <- run_sim(burn_in=burn_in,theta0=theta0,
#                       theta1=0,sigma=sigma,
#                       batch_size=25,batches=100, 
#                       eps=eps, clip=clip,r=100)
# 
# 
# null_trial <- pvals_null$trials
# null_pvals <- pvals_null$pval
# null_pvals$T_stat <- qnorm(1-null_pvals$pval)
# 
# 
# ggplot(data=null_pvals) +
#   geom_histogram(aes(x=pval),bins=10) +
#   facet_grid(cols=vars(method),rows=vars(rand_method))
# 
# ggplot(data=null_trial) +
#   geom_point(aes(x=subject,y=pi)) +
#   facet_wrap(vars(method)) +
#   labs(x="P(A=1)",y="") +
#   theme_bw()
# 
# 
# 
# ggplot(data=null_trial) + 
#   geom_histogram(aes(x=pi),bins=10) +
#   facet_wrap(vars(method)) +
#   scale_y_discrete()+
#   labs(x="P(A=1)",y="") +
#   theme_bw()
# 
# null_pvals %>% 
#   mutate(reject=ifelse(pval<0.05,1,0)) %>%
#   group_by(method,rand_method,theta1) %>%
#   summarise(power=mean(reject))
# 
# 
# 
# 
# 
# ## power analysis
# pval_summary <- pvals %>% 
#                 mutate(reject=ifelse(pval<0.05,1,0)) %>%
#                 group_by(method,rand_method,theta1) %>%
#                 summarise(power=mean(reject))
# 
# 
# ggplot(data=pval_summary %>% filter(rand_method %in% c("SR","RSR"))) +
#   geom_line(aes(x=theta1,y=power,color=rand_method)) +
#   geom_hline(yintercept = 0.05) +
#   facet_wrap(vars(method)) + 
#   theme_minimal()
# 
# ggplot(data=pval_summary %>% filter(rand_method %in% c("EG","REG"))) +
#   geom_line(aes(x=theta1,y=power,color=rand_method)) +
#   geom_hline(yintercept = 0.05) +
#   facet_wrap(vars(method)) + 
#   theme_minimal()
# 
# ggplot(data=pval_summary %>% filter(rand_method %in% c("TS","RTS"))) +
#   geom_line(aes(x=theta1,y=power,color=rand_method)) +
#   geom_hline(yintercept = 0.05) +
#   facet_wrap(vars(method)) + 
#   theme_minimal()
