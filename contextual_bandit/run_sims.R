## ----------------------------------------------------------------- ##
## run_sims.R ------------------------------------------------------ ##
## Author: Peter Norwood, NC State University ---------------------- ##
## Purpose: run a simulation experiment with certain inputs -------- ##
## ----------------------------------------------------------------- ##

## load functions
setwd("~/Research/NonDomSeqExp/NonDomSeqExp/contextual_bandit")
source("sims.R")
library(parallel)
library(gridExtra)
library(grid)

## run the simulations
## parameters
p=10
K=5

N=1000
sd_X=1
sd_Y=0.5
eps=0.05
percentile=0.95

## how many simulations to run
r <- 5
## Simulate the experiments
start <- Sys.time()
sims <- mclapply(X=1:r, 
                 function(X){sim(N=N,p=p,K=K,sd_X=sd_X,sd_Y=sd_Y,
                                 eps=eps,percentile=percentile)},
                 mc.cores = 1
)
end <- Sys.time()


## check bad sims
# good_sets <- c()
# tick=1
# for(i in 1:r){
#   check <- c()
#   for(j in 1:6){
#     check[j] <- is.data.frame(sims[[i]][[j]])
#   }
#   if(sum(check)==6){
#     good_sets[tick]=i
#     tick=tick+1
#   }else{
#     ## no nothing
#   }
# }

## collect experiment information into one dataframe
exps <- data.frame()
post <- data.frame()
tick <- 1
for(i in 1:r){
  
  ## e_greedy
  exp_e_greedy <- sims[[i]]$e_greedy
  exp_e_greedy$method <- "e_greedy"
  exp_e_greedy$rep <- tick
  
  ## nondom_e_greedy
  exp_nondom_e_greedy <- sims[[i]]$nondom_e_greedy
  exp_nondom_e_greedy$method <- "nondom_e_greedy"
  exp_nondom_e_greedy$rep <- tick
  
  ## simple
  exp_simple <- sims[[i]]$simple
  exp_simple$method <- "simple"
  exp_simple$rep <- tick
  
  ## nondom_simple
  exp_nondom_simple <- sims[[i]]$nondom_simple
  exp_nondom_simple$method <- "nondom_simple"
  exp_nondom_simple$rep <- tick
  
  ## ucb
  exp_ucb <- sims[[i]]$ucb
  exp_ucb$method <- "ucb"
  exp_ucb$rep <- tick
  
  ## nondom_ucb
  exp_nondom_ucb <- sims[[i]]$nondom_ucb
  exp_nondom_ucb$method <- "nondom_ucb"
  exp_nondom_ucb$rep <- tick
  
  ## ts
  exp_ts <- sims[[i]]$ts
  exp_ts$method <- "ts"
  exp_ts$rep <- tick
  
  ## nondom_ts
  exp_nondom_ts <- sims[[i]]$nondom_ts
  exp_nondom_ts$method <- "nondom_ts"
  exp_nondom_ts$rep <- tick
  
  
  
  ## create out of trial information
  ## attach these datasets to the bigger one
  exps <- rbind(exps,
                exp_e_greedy,exp_nondom_e_greedy,
                exp_simple,exp_nondom_simple,
                exp_ucb,exp_nondom_ucb,
                exp_ts,exp_nondom_ts)
  
  
  ## update the tick
  tick=tick+1
  
}

## save trial information
exps$K<- K
exps$p <- p

save(exps,file=paste0("exps_K=",K,"_p=",p,".RData"))
#save(post,file=paste0("post_K=",K,"_p=",p,".RData"))

burn_in=(p+1)*K*3

nondom_plot <- exps %>%
  filter(rep,sub>burn_in) %>%
  group_by(method,rep) %>%
  mutate(cum_regret=cumsum(regret),
         cum_nondom = cumsum(non_dom),
         nondom_prop = cum_nondom/(sub-burn_in),
         method_char = case_when(method=="simple"~"SR",
                                 method=="nondom_simple"~"RSR",
                                 method=="ucb"~"UCB",
                                 method=="nondom_ucb"~"RUCB",
                                 method=="ts"~"TS",
                                 method=="nondom_ts"~"RTS",
                                 method=="e_greedy"~"EG",
                                 method=="nondom_e_greedy"~"REG"),
         method_type = case_when(method %in% c("simple","nondom_simple")~"Simple",
                                 method %in% c("ucb","nondom_ucb")~"UCB",
                                 method %in% c("ts","nondom_ts")~"TS",
                                 method %in% c("e_greedy","nondom_e_greedy")~"E-Greedy")) %>%
  group_by(method_type,method_char,sub) %>%
  summarise(nondom_prop=mean(nondom_prop),
            cum_regret=mean(cum_regret),
            norm=mean(norm))

eg_plot <- ggplot(data=nondom_plot %>% filter(method_type=="E-Greedy")) +
  geom_line(aes(x=sub,y=nondom_prop,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#4B9CD3","#CC0000")) +
  labs(x="Subject",y="Proportion",color="Method") +
  theme_bw()

simple_plot <- ggplot(data=nondom_plot %>% filter(method_type=="Simple")) +
  geom_line(aes(x=sub,y=nondom_prop,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Proportion",color="Method") +
  theme_bw()

ucb_plot <- ggplot(data=nondom_plot %>% filter(method_type=="UCB")) +
  geom_line(aes(x=sub,y=nondom_prop,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Proportion",color="Method") +
  theme_bw()

ts_plot <- ggplot(data=nondom_plot %>% filter(method_type=="TS")) +
  geom_line(aes(x=sub,y=nondom_prop,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Proportion",color="Method") +
  theme_bw()

grid.arrange(eg_plot,simple_plot,ts_plot,ucb_plot, ncol=2, nrow=2,
             top = textGrob("Proportion Non-Dominated Assignment",gp=gpar(fontsize=20,font=1)))


## regret
eg_regret_plot <- ggplot(data=nondom_plot %>% filter(method_type=="E-Greedy")) +
  geom_line(aes(x=sub,y=cum_regret,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#4B9CD3","#CC0000")) +
  labs(x="Subject",y="Regret",color="Method") +
  theme_bw()

simple_regret_plot <- ggplot(data=nondom_plot %>% filter(method_type=="Simple")) +
  geom_line(aes(x=sub,y=cum_regret,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Regret",color="Method") +
  theme_bw()

ucb_regret_plot <- ggplot(data=nondom_plot %>% filter(method_type=="UCB")) +
  geom_line(aes(x=sub,y=cum_regret,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Regret",color="Method") +
  theme_bw()

ts_regret_plot <- ggplot(data=nondom_plot %>% filter(method_type=="TS")) +
  geom_line(aes(x=sub,y=cum_regret,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Regret",color="Method") +
  theme_bw()

grid.arrange(eg_regret_plot,simple_regret_plot,ts_regret_plot,ucb_regret_plot, ncol=2, nrow=2,
             top = textGrob("Cumulative Regret",gp=gpar(fontsize=20,font=1)))

## efficiency
eg_norm_plot <- ggplot(data=nondom_plot %>% filter(method_type=="E-Greedy")) +
  geom_line(aes(x=sub,y=norm,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#4B9CD3","#CC0000")) +
  labs(x="Subject",y="Norm",color="Method") +
  theme_bw()

simple_norm_plot <- ggplot(data=nondom_plot %>% filter(method_type=="Simple")) +
  geom_line(aes(x=sub,y=norm,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Norm",color="Method") +
  theme_bw()

ucb_norm_plot <- ggplot(data=nondom_plot %>% filter(method_type=="UCB")) +
  geom_line(aes(x=sub,y=norm,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Norm",color="Method") +
  theme_bw()

ts_norm_plot <- ggplot(data=nondom_plot %>% filter(method_type=="TS")) +
  geom_line(aes(x=sub,y=norm,group=method_char,color=method_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="Subject",y="Norm",color="Method") +
  theme_bw()

grid.arrange(eg_norm_plot,simple_norm_plot,ts_norm_plot,ucb_norm_plot, ncol=2, nrow=2,
             top = textGrob("Efficiency",gp=gpar(fontsize=20,font=1)))
