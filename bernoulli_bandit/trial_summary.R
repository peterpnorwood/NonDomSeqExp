## ----------------------------------------------------------------- ##
## trial_summary.R ------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: analyze trial information for bernoulli bandit --------- ##
## ----------------------------------------------------------------- ##


setwd("-------")

## load trial information
load("-------.RData")


## load packages
library(ggplot2)
library(dplyr)
library(knitr)
library(gridExtra)
library(grid)


## minor updates
trials <- trials %>%
          mutate(method_char=case_when(
                  method=="simple"~"SRS",
                  method=="e_greedy_10"~"EG(0.1)",
                  method=="e_greedy_25"~"EG(0.25)",
                  method=="kl_e_greedy_10"~"REG(0.1)",
                  method=="kl_e_greedy_25"~"REG(0.25)",
                  method=="kl_simple"~"RSR"),
                 correct=I(A==opt_A))


## proportion of patients (after burn in) randomized
## to optimal treatment
burn_in=25
trial_info <- trials %>%
 filter(subject>burn_in & subject < 500) %>%
 group_by(method_char,rep) %>%
 summarise(mean_Y1=mean(Y),
           mean_opt=mean(correct)) %>%
 group_by(method_char) %>%
 summarise(mean_Y=mean(mean_Y1),
           sd_Y=sd(mean_Y1),
           prop_opt=mean(mean_opt),
           sd_opt=sd(mean_opt),
           prop_5=quantile(mean_opt,0.05),
           prop_25=quantile(mean_opt,0.25))
 
save(trial_info,file="------.RData")

## proportion of patients randomized to non-dominated trts
non_dom_info <- trials %>%
  filter(subject>burn_in) %>%
  group_by(method_char,rep) %>%
  summarise(non_dom1=mean(non_dom)) %>%
  group_by(method_char) %>%
  summarise(non_dom=mean(non_dom1),
            sd_non_dom=sd(non_dom1),
            non_dom_5=quantile(non_dom1,0.05),
            non_dom_25=quantile(non_dom1,0.25))


save(non_dom_info,file="------.RData")


## cumulative regret
regret <- trials %>%
          filter(subject>burn_in) %>%
          group_by(method,rep) %>%
          mutate(cum_regret=cumsum(regret_p)) %>%
          group_by(method,subject) %>%
          summarise(mean_cum_regret=mean(cum_regret))

save(regret,file="regret.RData")

## cumulative proportion of correct trt assignment
cum_prop <- trials %>%
            filter(subject>burn_in) %>%
            group_by(method,rep) %>%
            mutate(cum_correct=cumsum(correct),
                   cum_prop=cum_correct/(subject-burn_in)) %>%
            group_by(method,subject) %>%
            summarise(mean_cum_prop=mean(cum_prop))

save(cum_prop,file="------.RData")


## tables for opt trt and non-dom trt
print(xtable(trial_info,digits=3),include.rownames = FALSE)
print(xtable(non_dom_info,digits=3),include.rownames = FALSE)

## regret plots
regret <- regret %>%
  mutate(method_char=case_when(
    method=="simple"~"SRS",
    method=="e_greedy_10"~"EG(0.1)",
    method=="e_greedy_25"~"EG(0.25)",
    method=="kl_e_greedy_10"~"REG(0.1)",
    method=="kl_e_greedy_25"~"REG(0.25)",
    method=="kl_simple"~"RSR"))

regret_1 <- ggplot(data=regret %>% filter(method_char %in% c("EG(0.25)","REG(0.25)"))) +
  geom_line(aes(x=subject,y=mean_cum_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) + 
  scale_linetype_manual(values=c("solid","dotted")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

regret_2  <- ggplot(data=regret %>% filter(method_char %in% c("EG(0.1)","REG(0.1)"))) +
  geom_line(aes(x=subject,y=mean_cum_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) + 
  scale_linetype_manual(values=c("solid","dotted")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

regret_3 <- ggplot(data=regret %>% filter(method_char %in% c("SRS","RSR"))) +
  geom_line(aes(x=subject,y=mean_cum_regret,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) + 
  scale_linetype_manual(values=c("solid","dotted")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

grid.arrange(regret_1,regret_2,regret_3,ncol=3,
             top=textGrob("Average Cumulative Regret", gp=gpar(fontsize=17)))

## cumulative proportion of assigning opt trt plot
cum_prop <- cum_prop %>%
  mutate(method_char=case_when(
    method=="simple"~"SRS",
    method=="e_greedy_10"~"EG(0.1)",
    method=="e_greedy_25"~"EG(0.25)",
    method=="kl_e_greedy_10"~"REG(0.1)",
    method=="kl_e_greedy_25"~"REG(0.25)",
    method=="kl_simple"~"RSR"))

opt_1 <- ggplot(data=cum_prop %>% filter(method_char %in% c("EG(0.25)","REG(0.25)"))) +
  geom_line(aes(x=subject,y=mean_cum_prop,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) + 
  ylim(0.45,1) +
  scale_linetype_manual(values=c("solid","dotted")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

opt_2  <- ggplot(data=cum_prop %>% filter(method_char %in% c("EG(0.1)","REG(0.1)"))) +
  geom_line(aes(x=subject,y=mean_cum_prop,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) +
  ylim(0.45,1) +
  scale_linetype_manual(values=c("solid","dotted")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

opt_3 <- ggplot(data=cum_prop %>% filter(method_char %in% c("SRS","RSR"))) +
  geom_line(aes(x=subject,y=mean_cum_prop,linetype=method_char),size=1.25) +
  labs(x="Subject", y="Cumulative Regret",linetype="Method") +
  xlim(c(25,500)) + 
  ylim(0.45,1) +
  scale_linetype_manual(values=c("dotted","solid")) +
  #ggtitle("Average Randomization Probabilities to Treatment 4 (Effective, Not Superior)") +
  theme_bw()

grid.arrange(opt_1,opt_2,opt_3,ncol=3,
             top=textGrob("Average Cumulative Probability of Assigning Optimal Treatment", gp=gpar(fontsize=17)))