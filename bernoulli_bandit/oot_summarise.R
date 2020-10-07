## ----------------------------------------------------------------- ##
## oot_summarise.R ------------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: summarise oot information ------------------------------ ##
## ----------------------------------------------------------------- ##

library(dplyr)
library(xtable)

setwd("-------")

## grab all of the datasets
load("-----.RData")
unif <- big_df
load("-----.RData")
left <- big_df
load("-----.RData")
right <- big_df

## combine the datasets
oot <- rbind(unif,left,right)

## analyze the rates of opt treatment
rates <- oot  %>%
  mutate(method_char=case_when(method=="simple"~"SRS",
                               method=="e_greedy_10"~"EG(0.1)",
                               method=="e_greedy_25"~"EG(0.25)",
                               method=="kl_e_greedy_10"~"REG(0.1)",
                               method=="kl_e_greedy_25"~"REG(0.25)",
                               method=="kl_simple"~"RSR")) %>%
        group_by(N,method_char,dist) %>%
        summarise(mean_y=mean(success),
                  sd_y=sd(success),
                  mean_correct=mean(correct),
                  sd_correct=sd(correct)) %>%
        arrange(N,method_char)

View(rates %>% filter(dist=="left"))
print(xtable(rates %>% filter(dist=="left") %>% dplyr::select(-dist),digits=3),include.rownames=FALSE)


## same distribution as before
View(rates %>% filter(dist=="unif"))
View(rates %>% filter(dist=="left")) ## most different than before


## Power analysis
load("------.RData")

power <- power_df %>%
  mutate(method_char=case_when(method=="simple"~"SRS",
    method=="e_greedy_10"~"EG(0.1)",
    method=="e_greedy_25"~"EG(0.25)",
    method=="kl_e_greedy_10"~"REG(0.1)",
    method=="kl_e_greedy_25"~"REG(0.25)",
    method=="kl_simple"~"RSR")) %>%
  group_by(method_char,N) %>%
  summarise(type_one=mean(r1),
            power_low=mean(r2),
            power_high=mean(r3)) %>%
  arrange(N,method_char) %>%
  dplyr::select(N,method_char,type_one,power_low,power_high)

print(xtable(power,digits=3),include.rownames = FALSE)

