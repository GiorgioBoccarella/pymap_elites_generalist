library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("~/Documents/trial/#2")


# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="")
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1',  "ps", "gen", "sim",  "transfer_in", "transfer_n", "inv_with",
                       "inv_T", 'p')
# Count p sum
archive$p <- archive$p0 + archive$p1


# Evolution
archive_evo <- archive[archive$transfer_in == 0,]


ggplot() +
  geom_line(data = archive_evo, aes(x = gen, y = trade_off, group = interaction(transfer_in, transfer_n, ps, env_pair, inv_T),
                                    color = factor(env_pair), linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  #xlab("Mutation steps") +
  #ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion"))



archive_1.3 <- archive[archive$env_pair == 1.3,]

#Transfer
archive_transfer <- archive[archive$transfer_in != 0,]
archive_transfer_1.3 <- archive[archive$transfer_in == 1.3,]


ggplot() +
  geom_line(data = archive_evo, aes(x = gen, y = fitness, group = interaction(transfer_in, transfer_n, ps, env_pair, inv_T),
                                   color = factor(env_pair), linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  #xlab("Mutation steps") +
  #ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion"))



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

