# Invasion ----------------------------------------------------------------
setwd("~/Documents/result_wrap/test_inv_rate_new")

inv_d <- read.table("invasion_id.dat", quote = "\"", comment.char = "")
colnames(inv_d) <- c("Invader", "Wild_type", "gen", "sim", "inv_rate", 'rate', 'psi', 'psw', 'ip0', 'ip1', 'wp0', 'wp1', 'e1', 'e2')

inv_d$ip <- inv_d$ip0 + inv_d$ip1
inv_d$wp <- inv_d$wp0 + inv_d$wp1

inv_d$env_l <- with(inv_d, paste0(e1, e2))

inv_d_1 <- inv_d[inv_d$inv_rate == 1,] 

ggplot(data=inv_d, aes(x=as.factor(Invader), fill = factor(Invader))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  facet_wrap(env_l~inv_rate, ncol = 2, scale = 'free_x')

ggplot(data=inv_d, aes(x=as.factor(Wild_type), fill = factor(Invader))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  facet_wrap(inv_rate~., ncol = 2, scale = 'free_x')

# 1 -----------------------------------------------------------------------


# Load_data
archive_1 <- read.table("archive.dat", quote="\"", comment.char="")
colnames(archive_1) <- c("env_pair", "fitness", "modularity","trade_off", 'p0',
                         'p1', "ps", "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p', 'e1', 'e2')

archive_1$p <- archive_1$p0 + archive_1$p1
archive_1$env_l <- with(archive_1, paste0(e1, e2))
archive_1 <- archive_1[archive_1$env_l == "[0.1,1.3]",]

f <- ggplot() +
  geom_line(data = archive_1, aes(x = gen, y = fitness, group = interaction(transfer_in,
                                                                            transfer_n, env_pair, inv_T),
                                  color = factor(env_pair), linetype = factor(inv_T)),
            stat="smooth",method = 'gam', alpha = 0.7, 
            arrow = arrow(length = unit(0.25, "cm"), ends="last", type = "closed"), size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("No Invasion", "Invasion")) +
  theme(legend.position = "none") +
  facet_wrap(env_l~inv_with)

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
  theme(legend.position = "none")+
  facet_wrap(env_l~inv_with)


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
  theme_gray(base_size = 15) +
  facet_wrap(env_l~inv_with)



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
  theme(legend.position = "none") +
  facet_wrap(env_l~inv_with)


r <- grid.arrange(f, m, t, p, nrow = 2)
