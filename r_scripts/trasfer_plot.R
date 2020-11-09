library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

setwd("~/Documents/trial/#S25")


# 1 -----------------------------------------------------------------------


# Load_data
archive_1 <- read.table("archive.dat", quote="\"", comment.char="", nrows = 150300)
colnames(archive_1) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B','f3avg', 'f3A', 'f3B',
                       "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive_1$p <- archive_1$p0 + archive_1$p1

f <- ggplot() +
  geom_line(data = archive_1, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                    transfer_n, env_pair, inv_T), color = factor(env_pair),
                                    linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_1, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                    linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_1, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                                 transfer_n, env_pair, inv_T), color = factor(env_pair),
                                    linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_1, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                    linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_1.png", plot=r, height = 10, width = 14,  dpi=300)


# 2 -----------------------------------------------------------------------

# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="", skip = 150300)
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                       "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive$p <- archive$p0 + archive$p1

archive_2 <- archive[archive$inv_with == 2,]

f <- ggplot() +
  geom_line(data = archive_2, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_2, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                               transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_2, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                      transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_2, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_2.png", plot=r, height = 10, width = 14,  dpi=300)




# 3 -----------------------------------------------------------------------

# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="", skip = 150300)
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                         'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                         "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive$p <- archive$p0 + archive$p1

archive_3 <- archive[archive$inv_with == 3,]

f <- ggplot() +
  geom_line(data = archive_3, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_3, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                               transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_3, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                      transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_3, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_3.png", plot=r, height = 10, width = 14,  dpi=300)

# 4 -----------------------------------------------------------------------

# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="", skip = 150300)
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                       "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive$p <- archive$p0 + archive$p1

archive_4 <- archive[archive$inv_with == 4,]

f <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                               transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                      transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_4.png", plot=r, height = 10, width = 14,  dpi=300)

# 5 -----------------------------------------------------------------------

# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="", skip = 150300)
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                       "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive$p <- archive$p0 + archive$p1

archive_4 <- archive[archive$inv_with == 5,]

f <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                               transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                      transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_4, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_5.png", plot=r, height = 10, width = 14,  dpi=300)

# 6 -----------------------------------------------------------------------

# Load_data
archive <- read.table("archive.dat", quote="\"", comment.char="", skip = 150300)
colnames(archive) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                       'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                       "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive$p <- archive$p0 + archive$p1

archive_6 <- archive[archive$inv_with == 6,]

f <- ggplot() +
  geom_line(data = archive_6, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

m <- ggplot() +
  geom_line(data = archive_6, aes(x = gen, y = modularity, group = interaction(transfer_in,
                                                                               transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


p <- ggplot() +
  geom_line(data = archive_6, aes(x = gen, y = p, group = interaction(transfer_in,
                                                                      transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("p regulators") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) 


t <- ggplot() +
  geom_line(data = archive_6, aes(x = gen, y = trade_off, group = interaction(transfer_in,
                                                                              transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Trade-off") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none")


r <- grid.arrange(f, m, t, p, nrow = 2)

ggsave("my_graph_6.png", plot=r, height = 10, width = 14,  dpi=300)



# Invasion ----------------------------------------------------------------


inv_d <- read.table("invasion_id.dat", quote = "\"", comment.char = "")
colnames(inv_d) <- c("Invader", "Wild_type", "gen", "sim", "inv_rate", 'rate', 'psi', 'psw', 'ip0', 'ip1', 'wp0', 'wp1')

inv_d$ip <- inv_d$ip0 + inv_d$ip1
inv_d$wp <- inv_d$wp0 + inv_d$wp1

ggplot(data=inv_d, aes(x=as.factor(Invader), fill = factor(Invader))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  facet_wrap(inv_rate~., ncol = 2, scale = 'free_x')




# Transfer ----------------------------------------------------------------

# Load_data
archive_t <- read.table("archive.dat", quote="\"", comment.char="")
colnames(archive_t) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                         'p1', "ps", 'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B',
                         "gen", "sim",  "transfer_in", "transfer_n", "inv_with", 
                         "rate", "inv_T", 'p')

archive_t$p <- archive_t$p0 + archive_t$p1



ggplot() +
  geom_line(data = archive_t, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                  transfer_n, env_pair, inv_T), color = factor(env_pair),
                                  linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") 

ggplot() +
  geom_line(data = archive_t, aes(x = gen, y = fitness, group = interaction(transf, factor(transf_n), inv_T, size = 0.3),
                                  color = factor(transf), linetype=factor(inv_T)), stat="smooth", 
            method = 'lm', se = FALSE, alpha = 0.7, arrow = arrow(length = unit(0.25, "cm"),
                                                                  ends="last", type = "closed"), size = 0.9) +
  guides(col=guide_legend("Transfered in Δ E"), linetype = guide_legend("") ) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Modularity change") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  facet_wrap(env_pair~., ncol = 2) 

