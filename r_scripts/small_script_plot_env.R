library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("~/Documents/clustered_sim/transfer_long_correct_inv_trial_env")

#Trade-offs

env <- read.table("env_list.dat", quote="\"", comment.char="") 

env <- gather(env, V61)

env$env_pair <- rep(1:3, each = 20)
env$env_k <- rep(0:1, each = 10)
env$env_l <- rep(1:20, each = 1)


ggplot(env, aes(x = env_l, y = env_pair, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~env_k, scales = "free_x")



