## ----------------------------------------------------------------- ##
## power_analysis.R ------------------------------------------------ ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: analyze power from the trials -------------------------- ##
## ----------------------------------------------------------------- ##

#setwd("~/Research/Non_Dominated/anova_mab")
load("5_25_trials.RData")

library(dplyr)
library(ggplot2)
library(lsmeans)

## load methods and different sample sizes to iterate through
methods <- c("det_simple","det_thompson_1","det_thompson_c",
             "simple","thompson_1","thompson_c")
Ns <- c(50,100,250,500)
R <- 5000

## loop through to calculate p-values
info <- matrix(NA,nrow=length(Ns)*length(methods)*R,ncol=8)
i = 1
for(r in 1:R){
  for(m in methods){
    for(n in Ns){
      
      ## store summary info
      info[i,1] <- m
      info[i,2] <- n
      info[i,3] <- r
      
      ## subset data
      temp <- trials %>% filter(rep==r,method==m & subject<=n)
      fit <- lm(Y~as.factor(A)-1,data=temp)
  
      ## calculate contrasts
      leastsquare <- lsmeans(fit, "A")
      contrasts <- list(c(-1,1,0,0,0),
                        c(-1,0,1,0,0),
                        c(-1,0,0,1,0),
                        c(-1,0,0,0,1),
                        c(0,0,0,-1,1))
      cont <- as.data.frame(contrast(leastsquare, contrasts))
      
      ## p-val mu2>mu1
      info[i,4] <- 1-pz(cont$t.ratio[1])
      ## p-val mu3>mu1
      info[i,5] <- 1-pz(cont$t.ratio[2])
      ## p-val mu4>mu1
      info[i,6] <- 1-pz(cont$t.ratio[3])
      ## p-val mu5>mu1
      info[i,7] <- 1-pz(cont$t.ratio[4])
      ## p-val mu5>mu4
      info[i,8] <- 1-pz(cont$t.ratio[5])
      
      i=i+1
      
    }
  }
}

info <- data.frame(info)

## analyze info
colnames(info) <- c("method","n","rep","p_mu2","p_mu3","p_mu4","p_mu5","p_mu5_mu4")
info_summary <- info %>%
                mutate(n=as.numeric(as.character(n)),
                       p_mu2=as.numeric(as.character(p_mu2)),
                       p_mu3=as.numeric(as.character(p_mu3)),
                       p_mu4=as.numeric(as.character(p_mu4)),
                       p_mu5=as.numeric(as.character(p_mu5)),
                       p_mu5_mu4=as.numeric(as.character(p_mu5_mu4)),
                       r_mu2=ifelse(p_mu2<0.05,1,0),
                       r_mu3=ifelse(p_mu3<0.05,1,0),
                       r_mu4=ifelse(p_mu4<0.05,1,0),
                       r_mu5=ifelse(p_mu5<0.05,1,0),
                       r_mu5_mu4=ifelse(p_mu5_mu4<0.05,1,0)) %>%
                group_by(method,n) %>%
                summarise(r_mu2=mean(r_mu2),
                          r_mu3=mean(r_mu3),
                          r_mu4=mean(r_mu4),
                          r_mu5=mean(r_mu5),
                          r_mu5_mu4=mean(r_mu5_mu4))

## save info
save(info_summary,file="5_25_info_summary.RData")
load("5_25_info_summary.RData")

