library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("~/Documents/clustered_sim/for_tuesday_yes")

#Trade-offs


archive <- read.table("archive_sim_.dat", quote="\"", comment.char="")
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 
                       "f03avg", "f03A", "f03B", "f13avg", "f13A", "f13B", 'p0', 'p1',  "ps", "gen", 
                       "sim",  "transf", "transf_n", "inv_with",
                       "inv_T", 'p')

archive$p <- archive$p0 + archive$p1


setwd("~/Documents/clustered_sim/transfer_19_10")

#Trade-offs


archive_1 <- read.table("archive_sim_.dat", quote="\"", comment.char="")
colnames(archive_1) <- c("env_pair", "fitness", "modularity","trade_off", 
                       "f03avg", "f03A", "f03B", "f13avg", "f13A", "f13B", 'p0', 'p1',  "ps", "gen", 
                       "sim",  "transf", "transf_n", "inv_with",
                       "inv_T", 'p')

archive_1$p <- archive_1$p0 + archive_1$p1


archive <- rbind(archive, archive_1)

summary(archive)

#Transfer
#archive_t <- archive[archive$transf == 0.3 | archive$transf == 1.3,]

# Evolution
archive_0t <- archive[archive$transf != 0,]
#archive_0t <- archive_0t[archive_0t$ps != 50,]

archive_0t <- archive_0t[archive_0t$env_pair == 1.1,]
#archive_0t <- archive_0t[archive_0t$inv_T == 1,]

ggplot() +
  geom_line(data = archive_0t, aes(x = gen, y = fitness, group = interaction(transf, transf_n, ps, env_pair, inv_T),
                                   color = factor(env_pair), linetype = factor(inv_T)), stat="smooth",method = 'gam', alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  #xlab("Mutation steps") +
  #ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(transf~ps, scale = 'free_x') 



  ggplot() +
  geom_line(data = archive, aes(x = gen, y = trade_off, group = interaction(transf, transf_n, inv_T, ps, env_pair),
                                color = factor(env_pair), linetype=factor(ps)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(transf~inv_T, ncol = 2, scale = 'free') 



ggplot() +
  geom_line(data = archive, aes(x = gen, y = modularity, group = interaction(transf, transf_n, inv_T, ps, env_pair),
                                color = factor(env_pair), linetype=factor(ps)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(transf~inv_T, ncol = 2, scale = 'free') 

ggplot() +
  geom_line(data = archive, aes(x = gen, y = p, group = interaction(transf, transf_n, inv_T, ps, env_pair),
                                color = factor(env_pair), linetype=factor(ps)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Regulators (p)") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(transf~inv_T, ncol = 2, scale = 'free') 




ggplot() +
  geom_line(data = archive, aes(x = gen, y = trade_off, group = interaction(transf, transf_n, inv_T, ps, env_pair, sim),
                                color = factor(env_pair), linetype=factor(ps)), stat="smoo
            th",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness change") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(sim~., ncol = 2) 



ggplot() +
  geom_line(data = archive, aes(x = gen, y = modularity, group = interaction(transf, transf_n, inv_T, ps, env_pair, sim),
                                color = factor(env_pair), linetype=factor(ps)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness change") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(sim~., ncol = 2) 


ggplot() +
  geom_line(data = archive, aes(x = gen, y = fitness, group = interaction(transf, transf_n, inv_T),
                                color = factor(transf), linetype=factor(inv_T)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) + 
  xlab("Mutation steps") +
  ylab("Fitness change") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(env_pair~ps + invasion_rate, ncol = 4) 


ggplot() +
  geom_line(data = archive, aes(x = gen, y = modularity, group = interaction(transf, transf_n, inv_T),
                                color = factor(transf), linetype=factor(inv_T)), stat="smooth",method = 'gam', se = FALSE, alpha = 0.7,
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness change") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(env_pair~ps + invasion_rate) 














ggplot() +
  geom_line(data = archive, aes(x = gen, y = fitness, group = interaction(inv_T, env_pair)
                                   , linetype= factor(inv_T),  color = factor(env_pair),
                                   alpha = I(10)), size = 1.2, stat="smooth", method = "gam", se = TRUE) +
  facet_wrap(inv_T~ps)


x <- archive_0t %>% gather(env, fit, 5:10)


ggplot() +
  geom_line(data = x, aes(x = gen, y = fit, group = interaction(env, inv_T, env_pair), color = factor(env_pair),
                          linetype = factor(inv_T), alpha = I(10)), size = 1.2, stat="smooth", method = "gam", se = TRUE) +
  facet_wrap(env~.) +
  scale_color_viridis_d()


x_ <- x[x$env == '_9',]

ggplot() +
  geom_line(data = x_, aes(x = gen, y = fit, group = interaction(env_pair), color = factor(env_pair),
                          alpha = I(10)), size = 1.2, stat="smooth", method = "gam", se = TRUE) +
  scale_color_viridis_d()

