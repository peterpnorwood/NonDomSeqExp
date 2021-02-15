## ----------------------------------------------------------------- ##
## sims.R ---------------------------------------------------------- ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: simulation experiments with different methods ---------- ##
## and store their results in a user-friendly format --------------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/NonDomSeqExp/NonDomSeqExp/contextual_bandit")
source("funcs.R")
source("simple.R")
source("ts.R")
source("nondom_ts.R")
source("ucb.R")
source("nondom_ucb.R")
source("nondom_simple.R")
source("e_greedy.R")
source("nondom_e_greedy.R")

## load packages
library(parallel)
library(tidyverse)
library(MASS)

## sim
## Purpose: simulate all methods on the same dataset
## param N: total sample size
## param p: number of features in the context (not including intercept)
## param K: number of arms
## param sd_X: standard deviation of the context
## param sd_Y: scaling factor sd(Y)=(sd_Y)*sqrt(p+1)
## param t0: time for greedy first to stay with 
##           greedy before checking
## param eps: P(explore) = eps in greedy first if it switches
## return output: list of trial information
sim <- function(N,p,K,sd_X,sd_Y,
                eps,percentile){
  
  ## standard deviation of Y
  sd_Y <- sd_Y*sqrt(p)
  ## burn in period is (# of parameters in full model)*3
  burn_in <- (p+1)*K*3
  ## vector of treatments
  A <- 1:K
  
  ## randomly sample the mean parameters
  theta <- rnorm((p+1)*K,0,1)
  
  ## generate training set
  train_set <- gen_data(N=N,p=p,sd_X=sd_X,A=1:K,sd_Y=sd_Y,theta=theta)
  
  ## run simple rand sim
  simple_sim <- try(simple(train_set=train_set,burn_in=burn_in,
                               sd_Y=sd_Y,A=A,theta=theta))
  
  ## run nondom simple sim
  nondom_simple_sim <- try(nondom_simple(train_set=train_set,burn_in=burn_in,
                                  sd_Y=sd_Y,A=A,theta=theta))
  
  ## run e_greedy experiment
  e_greedy_sim <- try(e_greedy(train_set=train_set,burn_in=burn_in,
                      sd_Y=sd_Y,A=A,theta=theta,eps=eps))
  
  ## run nondom_e_greedy experiment
  nondom_e_greedy_sim <- try(nondom_e_greedy(train_set=train_set,burn_in=burn_in,
                               sd_Y=sd_Y,A=A,theta=theta,eps=eps))
  
  ## run ucb
  ucb_sim <- try(ucb(train_set=train_set,burn_in=burn_in,
                  sd_Y=sd_Y,A=A,theta=theta,percentile=percentile))
  
  ## run nondom ucb
  nondom_ucb_sim <- try(nondom_ucb(train_set=train_set,burn_in=burn_in,
                        sd_Y=sd_Y,A=A,theta=theta,percentile=percentile))
  
  ## run ts
  ts_sim <- try(ts(train_set=train_set,burn_in=burn_in,
                     sd_Y=sd_Y,A=A,theta=theta))
  
  ## run nondom ts
  nondom_ts_sim <- try(nondom_ts(train_set=train_set,burn_in=burn_in,
                                   sd_Y=sd_Y,A=A,theta=theta))
  
  
  
  ## save output in a list
  output <- list(train_set=train_set,
                 simple=simple_sim,
                 nondom_simple=nondom_simple_sim,
                 e_greedy=e_greedy_sim,
                 nondom_e_greedy=nondom_e_greedy_sim,
                 ucb=ucb_sim,
                 nondom_ucb=nondom_ucb_sim,
                 ts=ts_sim,
                 nondom_ts=nondom_ts_sim)
  
  
  return(output)
  
}
