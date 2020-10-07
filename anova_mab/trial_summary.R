## ----------------------------------------------------------------- ##
## trial_summary.R ------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: analyze trial information for anova/mab scenario ------- ##
## ----------------------------------------------------------------- ##

setwd("-------")

## load the information
load("-------.RData")
load("-------.RData")
load("-------.RData")
load("-------.RData")


## load pacakges
library(ggplot2)
library(dplyr)
library(knitr)
library(gridExtra)
library(grid)

## create effective and superior treatment labels
## give new names to methods for plotting
trials <- trials %>% 
           mutate(A_sup=ifelse(A==5,1,0),
                  A_eff=ifelse(A %in% c(4,5),1,0),
                  method_char=case_when(
                   method=="simple"~"SRS",
                   method=="thompson_1"~"TH(1)",
                   method=="thompson_c"~"TH(0.5)",
                   method=="det_thompson_1"~"RTH(1)",
                   method=="det_thompson_c"~"RTH(0.5)",
                   method=="det_simple"~"RSR"))

## proportion of subjects randomized to effective or superior txt
trials_summary <- trials %>% 
                   group_by(method_char,rep) %>%
                   summarise(prob_sup=mean(A_sup),
                             prob_eff=mean(A_eff))
                             
trials_summary %>% group_by(method_char) %>% 
   summarise(mean_sup=mean(prob_sup),
             sd_sup=sd(prob_sup),
             five_sup=quantile(prob_sup,0.05),
             quarter_sup=quantile(prob_sup,0.25),
             mean_eff=mean(prob_eff),
             sd_eff=sd(prob_eff),
             five_eff=quantile(prob_eff,0.05),
             quarter_eff=quantile(prob_eff,0.25))

ggplot(data=trials %>% group_by(method,rep) %>% summarise(Y=mean(Y))) + 
   geom_boxplot(aes(x=method,y=Y))

## overall mean response by method
trials %>% group_by(method) %>% summarise(mean(Y))

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## looking at parameter estimates
colnames(param_ests) <- c("sub1","mu1","mu2","mu3","mu4","mu5","subject","method","rep")
param_ests$sub1 <- NULL
# 
t_mu1=1
t_mu2=1
t_mu3=1
t_mu4=1.8
t_mu5=2.0


## mean squared error
param_ests <- param_ests %>% 
                   mutate(mse1 = (mu1-t_mu1)**2,
                          mse2 = (mu2-t_mu2)**2,
                          mse3 = (mu3-t_mu3)**2,
                          mse4 = (mu4-t_mu4)**2,
                          mse5 = (mu5-t_mu5)**2,
                          mse_54=( (mu5-mu4)-(t_mu5-t_mu4) )**2)
param_summary <- param_ests %>% 
                  group_by(subject,method) %>%
                  summarise(mu1=mean(mu1),
                            mu2=mean(mu2),
                            mu3=mean(mu3),
                            mu4=mean(mu4),
                            mu5=mean(mu5),
                            mse1=mean(mse1),
                            mse2=mean(mse2),
                            mse3=mean(mse3),
                            mse4=mean(mse4),
                            mse5=mean(mse5),
                            mse_54=mean(mse_54))

## plot of mse throughout the trial for difference of trt 5 and 4
ggplot(data=param_summary) +
  geom_line(aes(x=subject,y=mse_54,color=method))


## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## looking at randomization probabilities
colnames(rand_probs) <- c("sub1","p1","p2","p3","p4","p5","subject","method","rep")
rand_probs$sub1 <- NULL

rand_summary <- rand_probs %>% 
                 group_by(method,subject) %>%
                 summarise(p1=mean(p1),
                           p2=mean(p2),
                           p3=mean(p3),
                           p4=mean(p4),
                           p5=mean(p5),
                           p45=p4+p5) %>%
                 mutate(method_char=case_when(
                   method=="simple"~"SRS",
                   method=="thompson_1"~"TH(1)",
                   method=="thompson_c"~"TH(0.5)",
                   method=="det_thompson_1"~"RTH(1)",
                   method=="det_thompson_c"~"RTH(0.5)",
                   method=="det_simple"~"RSR"))


## plots of average rand probabilities to trt 4
four_1 <- ggplot(data=rand_summary %>% filter(method %in% c("det_thompson_c","thompson_c"))) +
  geom_line(aes(x=subject,y=p4,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Average Randomization Probability",linetype="Method") +
  ylim(c(0,0.5)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

four_2  <- ggplot(data=rand_summary %>% filter(method %in% c("det_thompson_1","thompson_1"))) +
  geom_line(aes(x=subject,y=p4,linetype=method_char),size=1.25) +
  labs(x="Subject",linetype="Method",y="") +
  ylim(c(0,0.5)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

four_3 <- ggplot(data=rand_summary %>% filter(method %in% c("det_simple","simple"))) +
  geom_line(aes(x=subject,y=p4,linetype=method_char),size=1.25) +
  labs(x="Subject",linetype="Method",y="") +
  ylim(c(0,0.5)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

grid.arrange(four_1,four_2,four_3,ncol=3,top=textGrob("Average Randomization Probability to Treatment Four (Effective, Not Superior)", gp=gpar(fontsize=17)))


## plots of average randomization probabilites to trt 5
five_1 <- ggplot(data=rand_summary %>% filter(method %in% c("det_thompson_c","thompson_c"))) +
  geom_line(aes(x=subject,y=p5,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Average Randomization Probability",linetype="Method") +
  ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

five_2  <- ggplot(data=rand_summary %>% filter(method %in% c("det_thompson_1","thompson_1"))) +
  geom_line(aes(x=subject,y=p5,linetype=method_char),size=1.25) +
  labs(x="Subject", y="",linetype="Method") +
  ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

five_3 <- ggplot(data=rand_summary %>% filter(method %in% c("det_simple","simple"))) +
  geom_line(aes(x=subject,y=p5,linetype=method_char),size=1.25) +
  labs(x="Subject", y="",linetype="Method") +
  ylim(c(0.15,1)) + 
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

grid.arrange(five_1,five_2,five_3,ncol=3,
             top=textGrob("Average Randomization Probability to Treatment Five (Superior)", gp=gpar(fontsize=17)))


## average randomization probabilites to 4 or 5
ggplot(data=rand_summary) +
        geom_line(aes(x=subject,y=p45,color=method_char),size=1.25) +
        labs(x="Subject", y="Average Randomization Probability",color="Method") +
        ylim(c(0.15,1)) + 
        ggtitle("Average Randomization Probabilities to Effective Treatments") +
        theme_bw()

## average randomization probabilites to 5
ggplot(data=rand_summary) +
       geom_line(aes(x=subject,y=p5,color=method_char),size=1.25) +
       labs(x="Subject", y="",color="Method") +
       ggtitle("Average Randomization Probabilities to Superior Treatment") +
       ylim(c(0.15,1)) + 
       theme_bw()



## proportion of subject randomized to non-dominated treatments
trials %>% filter(subject>7) %>% 
  group_by(method_char,rep) %>% 
  summarise(prop_non_dom=mean(non_dom)) %>% 
  group_by(method_char) %>%
  summarise(mean_non_dom=mean(prop_non_dom),
            sd_non_dom=sd(prop_non_dom),
            five_non_dom=quantile(prop_non_dom,0.05),
            quarter_non_dom=quantile(prop_non_dom,0.25)) 




## ------------------------------------------------------------------ ##
## ----------------------------------------------------------------- ##

## looking at standard errors
colnames(std_errors) <- c("sub1","se1","se2","se3","se4","se5","subject","method","rep")


error_summary <- std_errors %>%
                 group_by(method,subject) %>%
                 summarise(se1=mean(se1),
                           se2=mean(se2),
                           se3=mean(se3),
                           se4=mean(se4),
                           se5=mean(se5))


ggplot(data=error_summary) +
  geom_line(aes(x=subject,y=se4,group=method,color=method))


## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## calculating cumulative regret

cumulative <- trials %>%
          mutate(regret=max(mu)-mu) %>%
          group_by(method_char,rep) %>%
          mutate(cumulative_regret=cumsum(regret)) %>%
          group_by(method_char,subject) %>%
          summarise(cumulative_regret=mean(cumulative_regret))

## regret plots
regret_1 <- ggplot(data=cumulative %>% filter(method_char %in% c("TH(0.5)","RTH(0.5)"))) +
  geom_line(aes(x=subject,y=cumulative_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

regret_2  <- ggplot(data=cumulative %>% filter(method_char %in% c("TH(1)","RTH(1)"))) +
  geom_line(aes(x=subject,y=cumulative_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

regret_3 <- ggplot(data=cumulative %>% filter(method_char %in% c("SRS","RSR"))) +
  geom_line(aes(x=subject,y=cumulative_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  scale_linetype_manual(values=c("dotted","solid")) +
  theme_bw()

grid.arrange(regret_1,regret_2,regret_3,ncol=3,
             top=textGrob("Average Cumulative Regret", gp=gpar(fontsize=17)))



