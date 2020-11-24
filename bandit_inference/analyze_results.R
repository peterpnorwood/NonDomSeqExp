## ----------------------------------------------------------------- ##
## analyze_results.R ----------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: analyze simulation results ----------------------------- ##
## ----------------------------------------------------------------- ##

setwd("-----------------------------")
source("sims.R")

library(parallel)
library(tidyverse)
library(gridExtra)
library(grid)
## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## randomization plots
load("run_one.RData")
one_results <- results
one_null <- one_results
one_null_trial <- one_null$trials
one_null_trial <- one_null_trial
## re-level the method so it is in order on the plot
one_null_trial$method <- factor(one_null_trial$method,
                                levels=c("SR","RSR","TS","RTS"))
                                            
one_null_pvals <- one_null$pval
one_null_pvals$T_stat <- qnorm(1-one_null_pvals$pval)

## grab the last batch
## since all the rand probs are the same for the batch, we can just
## grab the first individual
last_batch_one <- one_null_trial %>% 
                  filter(rep<=5000 & batch==100 & subject==4951) %>%
                  na.omit()

## plotzz
rand_one <- ggplot(last_batch_one) +
  geom_histogram(aes(x=pi),bins=10) +
  facet_wrap(vars(method)) +
  scale_y_discrete()+
  labs(x="P(A=1)",y="") +
  ggtitle("Batch Size=50, Batches=100") +
  theme_bw()

## same thing for the second scenario
load("run_two.RData")
two_results <- results
two_null <- two_results
two_null_trial <- two_null$trials
two_null_trial$method <- factor(two_null_trial$method,levels=c("SR","RSR","TS","RTS"))

two_null_pvals <- two_null$pval
two_null_pvals$T_stat <- qnorm(1-two_null_pvals$pval)

last_batch_two <- two_null_trial %>% 
  filter(rep<=5000 & batch==50 & subject==4901) %>%
  na.omit()

rand_two <- ggplot(data=last_batch_two) +
  geom_histogram(aes(x=pi),bins=10) +
  facet_wrap(vars(method)) +
  scale_y_discrete()+
  labs(x="P(A=1)",y="") +
  ggtitle("Batch Size=100, Batches=50") +
theme_bw()

## arrange them in a grid
grid.arrange(rand_one,rand_two,ncol=2)

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

### type one error
one_pvals <- one_null_pvals %>% 
  filter(rep<=5000) %>%
  na.omit() %>%
  group_by(method,rand_method) %>%
  mutate(reject=ifelse(pval<0.05,1,0)) %>%
  summarise(type_one=mean(reject))

two_pvals <- two_null_pvals %>% 
  filter(rep<=5000) %>%
  na.omit() %>%
  group_by(method,rand_method) %>%
  mutate(reject=ifelse(pval<0.05,1,0)) %>%
  summarise(type_one=mean(reject))

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## power 
load("power_one.RData")
load("power_two.RData")

## scenario one
plot_one <- ggplot(data=power_one) +
  geom_line(aes(x=theta1,y=power,linetype=rand_method,color=rand_method),size=1) +
  scale_color_manual(values=c("#CC0000","#012169","#4B9CD3","#9E7E38")) +
  labs(x="Signal", y="Power",linetype="Method",color="Method") +
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Batch Size=100, Batches=50") +
  theme_bw()

## scenario two
plot_two <- ggplot(data=power_two) +
  geom_line(aes(x=theta1,y=power,linetype=rand_method,color=rand_method),size=1) +
  scale_color_manual(values=c("#CC0000","#012169","#4B9CD3","#9E7E38")) +
  labs(x="Signal", y="Power",linetype="Method",color="Method") +
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Batch Size=50, Batches=100") +
  theme_bw()

## put them all together in one plot

## Pull the legend function out
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

e_legend <- get_legend(plot_one)

## remove legends
plot_one <- plot_one + theme(legend.position="none")
plot_two <- plot_two + theme(legend.position="none")

## plot in a grid
grid.arrange(plot_one,plot_two,e_legend,ncol=3,
             widths=c(2.3, 2.3,0.8))

