library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggpmisc)
library(reshape2)
set.seed(125)

data_ImpVir = read_excel("9ImpVir_phage_info.xlsx", sheet = "9ImpVir_wNA")
data_ImpVir = filter(data_ImpVir, !is.na(norm_tpmean_RNA))
data_ImpVir = filter(data_ImpVir, NEC_onset == "Early")

data = filter(data_ImpVir, Days_before_onset<=8)
# data$norm_tpmean_DNA[data$norm_tpmean_DNA==0] = 0.00000000000001

plot_box3.1 <- ggplot(data, aes(x=as.factor(Days_before_onset), y=norm_tpmean_DNA, fill=NEC_status)) +
  geom_boxplot(outlier.size=0) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  # ylab("DNA") +
  xlab("Days before NEC onset") +
  # scale_y_log10() +
  scale_x_discrete(limits = c("8","7","6","5","4","3","2","1","0"), labels = c(-8:0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # geom_signif(comparisons = list(c("case", "control")),
  #             map_signif_level = F, textsize=5)
  #stat_compare_means() +
  # ylim(0,0.5) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7) +
  facet_wrap(~Phage_host, scale="free", nrow = 3)
plot_box3.1
# ggsave("250424_box_compare_DNA_ImpVir_early_DBO_v1.pdf", plot_box3.1, width=12, height=7)

##--------------------------------------------------------------------------------------------------------------
# data1 = filter(data_ImpVir, Phage_host == "Dickeya" | Phage_host == "Enterobacter" | Phage_host == "Klebsiella" | Phage_host == "Salmonella")
data1 = filter(data_ImpVir, Phage_host == "Klebsiella" | Phage_host == "Salmonella")
data2 = filter(data1, Days_before_onset<=8)

plot_box3.2 <- ggplot(data2, aes(x=as.factor(Days_before_onset), y=norm_tpmean_DNA, fill=NEC_status)) +
  geom_boxplot(outlier.size=0.05) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("DNA") +
  xlab("Days before NEC onset") +
  # scale_y_log10() +
  scale_x_discrete(limits = c("3","2","1","0"), labels = c(-3:0)) +
  # scale_x_discrete(limits = c("8","7","6","5","4","3","2","1","0"), labels = c(-8:0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # geom_signif(comparisons = list(c("case", "control")),
  #             map_signif_level = F, textsize=5)
  #stat_compare_means() +
  ylim(0,40000) +
  stat_compare_means(method = "wilcox.test", label ="p.format", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7) +
  facet_wrap(~Phage_host, scale="fixed", nrow = 1)
plot_box3.2
ggsave("250424_box_compare_DNA_ImpVir_early_DBO_v3.pdf", plot_box3.2, width=6.2, height=3)
