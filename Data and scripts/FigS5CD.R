library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggstatsplot)
set.seed(125)

setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 2/Fig2C_extension')

metadata<-read.delim('241004_nec_metadata_144.txt', header = TRUE)
metadata = filter(metadata, !grepl("204-01|2156-01|415-01|2132-01|2234-01|234-01|
                                 421-01|2218-01|2238-02|49-01|257-01|2022-01|
                                 2077-01|2190-01|2228-01", Patient))
metadata = filter(metadata, Patient!="421-01" & Patient!="2077-01")
metadata$nec_stage[is.na(metadata$nec_stage)] = "Control"
# row.names(metadata)<-metadata[,1]
# metadata<-metadata[,-1]

metadata$NEC_onset2[metadata$NEC_onset>40] = "Late"
metadata$NEC_onset2[metadata$NEC_onset<=40] = "Early"

metadata_case = filter(metadata, NEC_status == 1)


##--------------------------------------------------------------------------------------------------------------------------------
metadata$NEC_status[metadata$NEC_status==1]="Case"
metadata$NEC_status[metadata$NEC_status==0]="Control"
plot_NEC_onset2_birthweight <- 
  ggplot(metadata, aes(x=as.factor(NEC_onset2), y=Birthweight, fill=as.factor(NEC_onset2))) +
  geom_boxplot() +
  geom_jitter(color="gray", size=2, alpha=1, width = 0.25) + #color="gray",
  # geom_jitter(color="gray", size=0.8, alpha=1, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
  # scale_fill_manual(values=c("Early"="#e6809d","Late"="#6788a6")) +
  scale_fill_manual(values=c("Early"="#a32b32","Late"="#d4edf8")) +
  xlab("NEC_onset2") +
  theme(legend.position="none", legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c("Early", "Late")),
              map_signif_level = F, textsize=4) +
# stat_compare_means(method = "wilcox.test", label ="p.format", label.x = 2) #p.format p.signif
# facet_wrap(~nec_stage, scale="fixed", nrow = 1)
  facet_wrap(~NEC_status, scale="fixed", nrow = 1)
plot_NEC_onset2_birthweight
ggsave("241005_Box_NEC_onset2_birthweight.pdf", plot_NEC_onset2_birthweight, width=4, height=5)

plot_NEC_onset2_GAbirth <- 
  ggplot(metadata, aes(x=as.factor(NEC_onset2), y=GA_birth, fill=as.factor(NEC_onset2))) +
  geom_boxplot() +
  geom_jitter(color="gray", size=2, alpha=1, width = 0.25) + #color="gray",
  # geom_jitter(color="gray", size=0.8, alpha=1, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
  # scale_fill_manual(values=c("Early"="#e6809d","Late"="#6788a6")) +
  scale_fill_manual(values=c("Early"="#a32b32","Late"="#d4edf8")) +
  xlab("NEC_onset2") +
  theme(legend.position="none", legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c("Early", "Late")),
              map_signif_level = F, textsize=4) +
# stat_compare_means(method = "wilcox.test", label ="p.format", label.x = 2) #p.format p.signif
  facet_wrap(~NEC_status, scale="fixed", nrow = 1)
plot_NEC_onset2_GAbirth
ggsave("241005_Box_NEC_onset2_GAbrith.pdf", plot_NEC_onset2_GAbirth, width=4, height=5)


##--------------------------------------------------------------------------------------------------------------------------------
scatter_BW_GA <-
  ggplot(metadata, aes(x=Birthweight, y=GA_birth, color=NEC_onset2)) +
  geom_point() +
  scale_color_manual(values=c("Early"="#e6809d","Late"="#6788a6")) +
  scale_fill_manual(values=c("Early"="#e6809d","Late"="#6788a6")) +
  xlab("NEC_onset2") +
  theme(legend.position=c(0.85,0.15), legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method=lm , color="light gray", fill="light gray", se=T, alpha = 0.2) +
  # geom_smooth(method=lm , aes(color=NEC_onset2, fill=NEC_onset2), se=T, alpha = 0.2) +
  stat_cor(label.y = 31) + 
  stat_regline_equation(label.y = 32)

scatter_BW_GA
ggsave("241005_Scatter_GAbrith.pdf", scatter_BW_GA, width=5, height=5)


