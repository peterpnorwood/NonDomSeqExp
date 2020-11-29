## ----------------------------------------------------------------- ##
## analyze_results.R ----------------------------------------------- ##
## Author(s): Peter Norwood, NCSU ---------------------------------- ##
## Purpose: analyze simulation results ----------------------------- ##
## ----------------------------------------------------------------- ##


setwd("-------------------------")
library(parallel)
library(tidyverse)
library(gridExtra)
library(grid)

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## randomization probabilites under null plot
load("fixed_pi.RData")
fixed_pi$method <- factor(fixed_pi$method,
                          levels=c("SR","RSR","TS","RTS"))
                                            
fixed_pi_null <- fixed_pi %>% filter(theta1==0)

rand_fixed <- ggplot(fixed_pi_null) +
  geom_histogram(aes(x=pi),bins=10) +
  facet_wrap(vars(method)) +
  scale_y_discrete()+
  labs(x="P(A=1)",y="") +
  ggtitle("Fixed Batch Size") +
  theme_bw()

load("random_pi.RData")
random_pi$method <- factor(random_pi$method,levels=c("SR","RSR","TS","RTS"))

random_pi_null <- random_pi %>% filter(theta1==0)

rand_random <- ggplot(random_pi_null) +
  geom_histogram(aes(x=pi),bins=10) +
  facet_wrap(vars(method)) +
  scale_y_discrete()+
  labs(x="P(A=1)",y="") +
  ggtitle("Random Batch Size") +
  theme_bw()

## put plots together
grid.arrange(rand_fixed,rand_random,ncol=2)

## ----------------------------------------------------------------- ##
## ----------------------------------------------------------------- ##

## type one error analysis
load("fixed_power.RData")

fixed_type_one <- fixed_power %>% filter(theta1==0)
fixed_type_one

load("random_power.RData")
random_type_one <- random_power %>% filter(theta1==0)
random_type_one


## power 
plot_fixed <- ggplot(data=fixed_power %>% filter(method=="bols")) +
  geom_line(aes(x=theta1,y=power,linetype=rand_method,color=rand_method),size=1) +
  scale_color_manual(values=c("#CC0000","#012169","#4B9CD3","#9E7E38")) +
  labs(x="Signal", y="Power",linetype="Method",color="Method") +
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Fixed Batch Size") +
  theme_bw()

plot_random <- ggplot(data=random_power %>% filter(method=="bols")) +
  geom_line(aes(x=theta1,y=power,linetype=rand_method,color=rand_method),size=1) +
  scale_color_manual(values=c("#CC0000","#012169","#4B9CD3","#9E7E38")) +
  labs(x="Signal", y="Power",linetype="Method",color="Method") +
  scale_linetype_manual(values=c("dotted","solid","dashed","longdash")) +
  ggtitle("Random Batch Size") +
  theme_bw()

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

e_legend <- get_legend(plot_fixed)

plot_one <- plot_fixed + theme(legend.position="none")
plot_two <- plot_random + theme(legend.position="none")

grid.arrange(plot_one,plot_two,e_legend,ncol=3,
             widths=c(2.3, 2.3,0.8))


