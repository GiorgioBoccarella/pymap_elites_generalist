library(ggplot2)
library(gtools)
library(dplyr)
library(ggpubr)

###########
#Import data
file_list <- list.files("/home/giorgio/Documents/results", pattern=".dat", full.names=T)
file_list <- mixedsort(sort(file_list))
data_m <- do.call("rbind", lapply(file_list, read.table))
colnames(data_m) <- c("fitness", "envID", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", 'St', "gen")
#colnames(data_m) <- c("fitness", "envID", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10","E11", "E12", "E13", "E14", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10","T11", "T12", "T13", "T14", 'St', "gen")

#Plot number of specialist and generalist

df2 <- data_m %>% group_by(gen) %>% count(St)

sg <- ggplot(df2, aes(x = gen, y = n, color = as.factor(St))) +
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  xlim(2000, 190000) +
  ylab("Number of individual") +
  theme(legend.position = "none")


#Plot fitness and standard deviation of Id

data_m$sd_t <- apply(data_m[13:22], 1, sd)

fit <- ggplot(data_m, aes(x = gen, y = fitness, group = envID, color = as.factor(St), alpha = I(0.2))) +
  geom_line() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "red")) +
  xlim(2000, 160000) 


st <- ggplot(data_m, aes(x = gen, y = sd_t, group = envID, color = as.factor(St), alpha = I(0.2))) +
  geom_line() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "red"))
#scale_color_manual(name= "Start as: ",labels = c("Generalist", "Specialist"), values = c("blue", "red"))

gst <- ggplot(data_m, aes(x = gen, y = sd_t, color = as.factor(St), alpha = I(0.2))) +
  geom_line() +
  scale_color_manual(name= "Start as: ",labels = c("Generalist", "Specialist"), values = c("blue", "red"))

#Merge plots
ggarrange(sg, fit, st, gst)
