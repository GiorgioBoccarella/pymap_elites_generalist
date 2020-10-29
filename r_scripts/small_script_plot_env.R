library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("~/Documents/trial/#2")

#Trade-offs
#If many sim a sample is necessary 


env <- read.table("env_list.dat", quote="\"", comment.char="") 
env$k <- rep(0:1, each = 100)
env$l <- rep(1:100)
colnames(env) <- c("value", "env_pair", "sim","k",'l')

ggplot(env, aes(x = l, weight = value, fill = factor(k), alpha = 0.1)) +
  geom_density() +
  facet_wrap(env_pair~., ncol =5 )

  



