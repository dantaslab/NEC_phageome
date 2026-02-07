library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggpmisc)
library(reshape2)
library(ggrepel)
set.seed(125)

data_ImpVir = read_excel("9ImpVir_phage_info.xlsx", sheet = "9ImpVir_wNA")
data_ImpVir = filter(data_ImpVir, !is.na(norm_tpmean_RNA))
data_ImpVir = filter(data_ImpVir, NEC_onset == "Early")

##--------------------------------------------------------------------------------------------------------------
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth>50] = "50-"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth<=50] = "40-50"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth<=40] = "30-40"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth<=30] = "20-30"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth<=20] = "10-20"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth<=10] = "1-10"
data_ImpVir$DOL_group[data_ImpVir$Date_after_birth==0] = "0"

plot_box1.2 <- ggplot(data_ImpVir, aes(x=DOL_group, y=RNA_to_DNA, fill=NEC_status)) +
  geom_boxplot(outlier.size=0.25) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("LOG(RNA/DNA)") +
  xlab("DOL") +
  # scale_y_log10() +
  # scale_x_discrete(limits = c("50-60","40-50","30-40","20-30","10-20","1-10","0")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # geom_signif(comparisons = list(c("case", "control")),
  #             map_signif_level = F, textsize=5)
  #stat_compare_means() +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7) +
  # facet_wrap(~Phage_host, scale="free", nrow = 2)
  facet_wrap(~Phage_host, scale="free", nrow = 3)
plot_box1.2
ggsave("250421_box_compare_RNA_DNA_ImpVir_early_DOL_v1.pdf", plot_box1.2, width=10, height=7)


# ##--------------------------------------------------------------------------------------------------------------
# 
# # data = filter(data_ImpVir, Days_before_onset<=7)
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset<=10] = "9-10"
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset<=8] = "7-8"
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset<=6] = "5-6"
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset<=4] = "3-4"
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset<=2] = "1-2"
# data_ImpVir$DBO_group2[data_ImpVir$Days_before_onset==0] = "0"
# 
# plot_box2.1 <- ggplot(data_ImpVir, aes(x=DBO_group2, y=RNA_to_DNA, fill=NEC_status)) +
#   geom_boxplot() +
#   #scale_fill_brewer(palette="RdBu") +
#   scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
#   #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#   #geom_jitter(color="gray", size=0.8, alpha=0.9) +
#   #theme_ipsum() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
#   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
#   ylab("RNA/DNA") +
#   xlab("Days before NEC onset") +
#   scale_y_log10() +
#   scale_x_discrete(limits = c("9-10","7-8","5-6","3-4","1-2","0")) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   # geom_signif(comparisons = list(c("case", "control")),
#   #             map_signif_level = F, textsize=5)
#   #stat_compare_means() +
#   stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.95,  hide.ns = TRUE, size = 5) +
#   facet_wrap(~Phage_host, scale="free", nrow = 3)
# plot_box2.1
# ggsave("250421_box_compare_RNA_DNA_ImpVir_early_DBO_v1.pdf", plot_box2.1, width=14, height=14)

##--------------------------------------------------------------------------------------------------------------
data = filter(data_ImpVir, Days_before_onset<=8)
# data$RNA_to_DNA[data$RNA_to_DNA==0] = 0.00000000000001

plot_box3.1 <- ggplot(data, aes(x=as.factor(Days_before_onset), y=RNA_to_DNA, fill=NEC_status)) +
  geom_boxplot(outlier.size=0) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  # ylab("LOG(RNA/DNA)") +
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
  ylim(0,0.5) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7) +
  facet_wrap(~Phage_host, scale="free", nrow = 3)
plot_box3.1
ggsave("250424_box_compare_RNA_DNA_ImpVir_early_DBO_v2.pdf", plot_box3.1, width=12, height=7)

##--------------------------------------------------------------------------------------------------------------
# data1 = filter(data_ImpVir, Phage_host == "Dickeya" | Phage_host == "Enterobacter" | Phage_host == "Klebsiella" | Phage_host == "Salmonella")
data1 = filter(data_ImpVir, Phage_host == "Klebsiella" | Phage_host == "Salmonella")
data2 = filter(data1, Days_before_onset<=8)

plot_box3.2 <- ggplot(data2, aes(x=as.factor(Days_before_onset), y=RNA_to_DNA, fill=NEC_status)) +
  geom_boxplot(outlier.size=0.05) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("RNA/DNA") +
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
  ylim(0,0.5) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7) +
  facet_wrap(~Phage_host, scale="fixed", nrow = 1)
plot_box3.2
ggsave("250424_box_compare_RNA_DNA_ImpVir_early_DBO_v3.pdf", plot_box3.2, width=6, height=3)

##--------------------------------------------------------------------------------------------------------------
#data_sample = filter(data_sample, tpmean_DNA_filtered>0)
plot_smooth1 <- ggplot(data_ImpVir, aes(x=Days_before_onset, y=RNA_to_DNA, color=NEC_status)) + #geom_point() +
  geom_smooth(method = loess, alpha=0.2) +
  #geom_smooth(alpha=0.2, span = 0.5) +
  # ylab("Phage RNA/DNA") +
  xlab("Days before NEC onset") +   # Postnatal age [days]
  # scale_color_discrete(name="NEC status") +
  scale_color_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  theme(legend.position="none", legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_reverse(limits = c(3, 0)) +
  facet_wrap(~Phage_host, scale="free", nrow = 3)
plot_smooth1
ggsave("250421_smooth_compare_RNA_DNA_ImpVir_early.pdf", plot_smooth1, width=14, height=7)


# ##--------------------------------------------------------------------------------------------------------------
# data_sample2= gather(data_sample, "TestGroup", "Norm_tpmean", norm_tpmean_RNA, norm_tpmean_DNA)
# plot2 <- ggplot(data_sample2, aes(x=Days_before_onset, y=Norm_tpmean, color=NEC_status)) + #geom_point() +
#   geom_smooth(method = loess, alpha=0.2) +
#   #geom_smooth(alpha=0.2, span = 0.5) +
#   ylab("Normallized coverage") +
#   xlab("Days before NEC onset") +   # Postnatal age [days]
#   # scale_color_discrete(name="NEC status") +
#   scale_color_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
#   theme(legend.position=c(0.8, 0.2), legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
#   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   scale_x_reverse(limits = c(60, 0)) +
#   facet_wrap(~TestGroup, scale="free", nrow = 1)
# plot2
# ggsave("241005_RA_Klebsiella_phage_bac.pdf", plot2, width=7, height=4)

##--------------------------------------------------------------------------------------------------------------
data_sample = filter(data2, Days_before_onset<=1)
plot_scatter <-
  ggplot(data_sample, aes(x=norm_tpmean_DNA, y=norm_tpmean_RNA)) +
  geom_point(aes(color = NEC_status), size = 3) +
  scale_color_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  # geom_smooth(method = loess, alpha=0.2) +
  # geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  # stat_poly_line(alpha=0.2, se = F, size = 1.5) +
  # stat_poly_eq(use_label(c("eq","adj.R2", "p")), label.x = 0.2) +
  # stat_poly_eq(use_label(c("eq", "adj.R2", "f", "p", "n"))) +
  # scale_color_manual(values = c("#293B5F")) +
  # geom_abline(intercept = 0, slope = 1, color = "orange", size = 1) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position=c(0.85,0.85)) +
  ylab("Phage RNA") +
  xlab("Phage DNA") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # coord_fixed( ratio=1) +
  facet_grid(Days_before_onset~Phage_host, scale="free") +
  geom_text_repel(
    aes(label = Contig),
    color = "black",
    size = 9/.pt, # font size 9 pt
    point.padding = 0.1, 
    box.padding = 0.6,
    min.segment.length = 2200,
    max.overlaps = 5,
    seed = 7654) # For reproducibility reasons
plot_scatter
ggsave("250424_Scatter_RNA_DNA_DBO0-1.pdf", plot_scatter, width=8, height=8)

