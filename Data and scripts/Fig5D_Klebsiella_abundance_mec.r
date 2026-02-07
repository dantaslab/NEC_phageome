library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr) 

data_mec = read_excel("virMec_phage_info.xlsx", sheet="virMec_sample")
data_bactaxa = read.csv('NEC_bacteria_genus_RA.csv', header = T)
data_bactaxa = select(data_bactaxa, SampleID, Klebsiella, Enterobacter, Serratia,  Atlantibacter)

data_mec_bactaxa = merge(data_bactaxa, data_mec, by.x="SampleID", by.y="Sample")
row.names(data_mec_bactaxa) = data_mec_bactaxa$SampleID
data_mec_bactaxa$SampleID = NULL

data_mec_bactaxa_early = filter(data_mec_bactaxa, NEC_onset2=="Early" & virMec_gene!="No9vir")
data_mec_bactaxa_early2 = filter(data_mec_bactaxa_early, Klebsiella_virMec!="NoVir")


plot_box <- ggplot(data_mec_bactaxa_early2, aes(x=Klebsiella_virMec, y=Klebsiella, fill=Klebsiella_virMec)) +
  geom_boxplot(outlier.size=0.25) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("Yes"="#ffcd58","No"="#345c7c"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="top") +
  ylab("Relative abundance") +
  xlab("Mec gene") +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size=4)+
  # scale_y_log10() +
  # scale_x_discrete(limits = as.factor(c(0:8)), labels = c(-8:0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = F, textsize=3)

plot_box
ggsave("250423_box_compare_Klebsiella_RA_mec_early.pdf", plot_box3.2, width=2.8, height=5)










