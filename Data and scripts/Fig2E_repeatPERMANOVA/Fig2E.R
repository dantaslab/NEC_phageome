library(ggplot2)
library(readxl)
library(dplyr)
set.seed(129)

###### Load and merge data #######

data_all <-read.csv("PERMANOVA_result_eachDBO_NEC_all.csv", header = T)
data_all$Group = "All"
data_early <-read.csv("PERMANOVA_result_eachDBO_NEC_early.csv", header = T)
data_early$Group = "Early"
data_late <-read.csv("PERMANOVA_result_eachDBO_NEC_late.csv", header = T)
data_late$Group = "Late"

Data = rbind(data_all, data_early, data_late)
Data$negLog10_P_value = -log10(Data$Pr..F.)
# Data$SigLabel[Data$Pr..F.<0.1] = "º"
Data$SigLabel[Data$Pr..F.<0.05] = "*"
Data$SigLabel[Data$Pr..F.<0.01] = "**"
Data$SigLabel[Data$Pr..F.<=0.001] = "***"

###### Plot  #######

plot1 <- ggplot(Data1, aes(x=Days_before_onset, y=R2)) + 
  # geom_line() +
  geom_smooth(span=2, alpha=0.1, color = "#004c4c") + #method = loess
  geom_point(alpha=0.7, size = 2, color = "#004c4c") +
  xlab("Days before NEC onset") +
  ylab("Effect size") + 
  scale_color_brewer(palette = "Spectral") +
  theme(legend.position="none", legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"),
        panel.grid.minor.y = element_blank()) +
  scale_x_reverse(limits = c(14, 0),breaks=c(0,2,4,6,8,10,12,14)) +
  ylim(-0.07,0.26) +
  geom_text(aes(label=SigLabel),nudge_x=0, nudge_y=0.02,check_overlap=F,size=8,color="#004c4c") +
  facet_grid(Group~Dataset, scale="fixed") +
  theme(strip.text = element_text(size = 12))
plot1
ggsave("PERMANOVA_NEC_early_late_line.pdf", plot1, width=15, height=8)





