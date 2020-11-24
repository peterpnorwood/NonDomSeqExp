## ----------------------------------------------------------------- ##
## run_two.R ------------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: run simulations for the two arm bandit ----------------- ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("sims.R")

library(parallel)
library(tidyverse)

## parameters 
burn_in=5
theta0=0
theta1=0
sigma=1
batch_size=100
batches=50
eps=0.1
clip=0.1

## reps
r=5250

trials <- data.frame()
pvals <- data.frame()

results <- run_sim(burn_in=burn_in,theta0=theta0,
                   theta1=theta1,sigma=sigma,
                   batch_size=batch_size,batches=batches, 
                   eps=eps, clip=clip,r=r)

save(results,file="run_one.RData")