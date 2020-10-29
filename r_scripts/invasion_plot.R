# Invasion ----------------------------------------------------------------


inv_d <- read.table("invasion_id.dat", quote = "\"", comment.char = "")
colnames(inv_d) <- c("Invader", "Wild_type", "gen", "sim", "inv_rate", 'rate', 'psi', 'psw', 'ip0', 'ip1', 'wp0', 'wp1')

inv_d$ip <- inv_d$ip0 + inv_d$ip1
inv_d$wp <- inv_d$wp0 + inv_d$wp1

ggplot(data=inv_d, aes(x=as.factor(Invader), fill = factor(Invader))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  facet_wrap(inv_rate~., ncol = 2, scale = 'free_x')
