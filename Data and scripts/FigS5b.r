library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)

setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/DataS')
sample_df = read_excel('Data_S2_v1.xlsx')
sample_df$NEC_onset = sample_df$DOL+sample_df$Days_before_NEC_onset
p_onset <- ggplot(sample_df, aes(x=NEC_onset)) + 
  geom_histogram(binwidth=5, fill="light gray", color="dark gray", alpha=1) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),legend.position="none") +
  ylab("Count") +
  xlab("NEC onset DOL") +
  # scale_y_log10() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept=mean(sample_df$NEC_onset), lwd=0.75, linetype=2, color="black") +
  geom_text(aes(label=round(mean(sample_df$NEC_onset),1), angle=90, y=270, x=mean(sample_df$NEC_onset)), vjust=1.2, col="black", size=4) +
  geom_vline(xintercept=median(sample_df$NEC_onset), lwd=0.75, linetype=2, color="#2293b3") +
  geom_text(aes(label=round(median(sample_df$NEC_onset),1), angle=90, y=270, x=median(sample_df$NEC_onset)), vjust=1.2, col="#2293b3", size=4)
p_onset
ggsave("250826_histogram_sample_NEConset_DOL_count.pdf", p_onset, width=4, height=4)
