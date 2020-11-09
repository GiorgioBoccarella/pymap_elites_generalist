library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(viridis)

setwd("~/Documents/result_wrap/base_line_p25")


# 1 -----------------------------------------------------------------------

# Load_data
archive_b <- read.table("archive.dat", quote="\"", comment.char="")
colnames(archive_b) <- c("env_pair", "fitness", "modularity","trade_off", 'p0','p1', "ps",
                         'f1avg', 'f1A', 'f1B', 'f2avg', 'f2A', 'f2B','f3avg', 'f3A', 'f3B',
                         'f4avg', 'f4A', 'f4B', 'f5avg', 'f5A', 'f5B','f6avg', 'f6A', 'f6B',
                         "gen", "sim",  "transfer_in", "transfer_n", "inv_with", "rate", "inv_T", 'p')

archive_b$p <- archive_b$p0 + archive_b$p1

long_archive <- archive_b %>% gather(env_fit, fit, 8:25)

# Filter for average

archive_avg <- dplyr::filter(long_archive, grepl("avg",env_fit))

e_names <- c(
  `f1avg` = "Fitness in 0.1 ",
  `f2avg` = "Fitness in 0.3 ",
  `f3avg` = "Fitness in 0.9 ",
  `f4avg` = "Fitness in 1.1 ",
  `f5avg` = "Fitness in 1.2 ",
  `f6avg` = "Fitness in 1.3 "
)



ggplot() +
  geom_line(data = archive_avg, aes(x = gen, y = fit, group = interaction(env_pair, env_fit),
                                    color = factor(env_pair)),
            stat="smooth",method = 'gam', alpha = 0.7, size = 0.9) +
  theme_gray(base_size = 15) +
  xlab("Mutation steps") +
  ylab("Fitness") +
  theme_gray(base_size = 15) +
  scale_color_viridis_d() +
  facet_wrap(~env_fit, labeller = as_labeller(e_names)) +
  guides(col=guide_legend("Δ E"), linetype = guide_legend(""), scale = "free" ) +
  xlim(0,250) +
  ylim(-0.85, -0.6)



archive_A <- dplyr::filter(long_archive, grepl("A",env_fit))
colnames(archive_A)[17] <- "fitA"

archive_B <- dplyr::filter(long_archive, grepl("B",env_fit))
colnames(archive_B)[17] <- "fitB"

archive_A$fitB <- archive_B$fitB

e_names_AB <- c(
  `f1A` = "0.1 ",
  `f2A` = "0.3 ",
  `f3A` = "0.9 ",
  `f4A` = "1.1 ",
  `f5A` = "1.2 ",
  `f6A` = "1.3 "
)

archive_A_1 <- archive_A[archive_A$env_pair == 1.3,] 
archive_A_1 <- archive_A[archive_A$sim == 0,] 


ggplot() +
  geom_point(data = archive_A_1, aes(x = fitA, y = fitB, group = interaction(env_pair, env_fit),
                                    color = factor(env_pair))) +
  theme_gray(base_size = 15) +
  xlab("Fit in A") +
  ylab("Fit in B") +
  theme_gray(base_size = 15) +
  scale_color_viridis_d() +
  facet_wrap(~factor(env_fit), labeller = as_labeller(e_names_AB), ncol = 3) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)

ggplot() +
  geom_line(data = archive_A, aes(x = fitA, y = fitB, group = interaction(env_pair, env_fit),
                                  color = factor(env_pair)),
            stat="smooth",method = 'lm', alpha = 0.7, size = 0.9)+
  theme_gray(base_size = 15) +
  xlab("Fit in A") +
  ylab("Fit in B") +
  theme_gray(base_size = 15) +
  scale_color_viridis_d() +
  facet_wrap(~factor(env_fit), labeller = as_labeller(e_names_AB), ncol = 3) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)



library(plotly)

archive_A_9 <- archive_A[archive_A$env_fit == "f5A",] 
plot_ly(x=archive_A_9$fitA, y=archive_A_9$fitB, z=archive_A_9$gen,
        type="scatter3d", mode="markers", color=archive_A_9$env_pair, size=I(8))

ggsave("my_graph_1.png", plot=f, height = 10, width = 14,  dpi=300)
