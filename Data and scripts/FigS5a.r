library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)

setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/DataS')
subject_df = read_excel('Data_S1_v1.xlsx')

p_onset <- ggplot(subject_df, aes(x=NEC_onset)) + 
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
  geom_vline(xintercept=mean(subject_df$NEC_onset), lwd=0.75, linetype=2, color="black") +
  geom_text(aes(label=round(mean(subject_df$NEC_onset),1), angle=90, y=23, x=mean(subject_df$NEC_onset)), vjust=1.3, col="black", size=4) +
  geom_vline(xintercept=median(subject_df$NEC_onset), lwd=0.75, linetype=2, color="#2293b3") +
  geom_text(aes(label=round(median(subject_df$NEC_onset),1), angle=90, y=23, x=median(subject_df$NEC_onset)), vjust=1.3, col="#2293b3", size=4)
p_onset
ggsave("250826_histogram_infant_NEConset_DOL_count.pdf", p_onset, width=4, height=4)

