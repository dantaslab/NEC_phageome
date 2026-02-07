library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(reshape2)

data_mec = read_excel("virMec_phage_info.xlsx", sheet="virMec_sample")
data_irep_genus = read.csv('mmc15_iRep_genus.csv', header = T, row.names = 1)
data_irep_genus = select(data_irep_genus, Klebsiella, Enterobacter, Serratia, Shigella, Atlantibacter)
data_irep_species = read.csv('mmc15_iRep_species.csv', header = T, row.names = 1)
data_irep_species = select(data_irep_species, Klebsiella_grimontii)

data_irep = merge(data_irep_genus,data_irep_species, by="row.names")

data_mec_irep = merge(data_irep, data_mec, by.x="Row.names", by.y="Sample")
row.names(data_mec_irep) = data_mec_irep$Row.names
data_mec_irep$Row.names = NULL

data_mec_irep_early = filter(data_mec_irep, NEC_onset2=="Early" & virMec_gene!="No9vir")
data_mec_irep_early2 = filter(data_mec_irep_early, Klebsiella_virMec!="NoVir")


plot_box <- ggplot(data_mec_irep_early2, aes(x=Klebsiella_virMec, y=Klebsiella, fill=Klebsiella_virMec)) +
  geom_boxplot(outlier.size=0.25) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("Yes"="#ffcd58","No"="#345c7c"))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="gray", size=0.8, alpha=0.9) +
  #theme_ipsum() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="top") +
  ylab("iRep") +
  xlab("Mec gene") +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5), size=4)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = F, textsize=3)

plot_box
ggsave("250423_box_compare_Klebsiella_irep_mec_early.pdf", plot_box, width=2.8, height=5)





