library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(reshape2)

data_mec = read_excel("virMec_phage_info.xlsx", sheet="virMec_sample")
data_mec$virMec_gene[data_mec$virMec_gene!="Yes"] = "No"

##----------------------------------------------------------------------------------------------------------------------------------------------------------------
data_mec.early.case = filter(data_mec, NEC_onset2=="Early" & NEC_status=="case")
mec_gene_DBO.early.case = as.data.frame.array(table(data_mec.early.case$Days_before_NEC_onset, data_mec.early.case$virMec_gene))
mec_gene_DBO.early.case$DBO = rownames(mec_gene_DBO.early.case)
mec_gene_DBO.early.case$sample_count = rowSums(mec_gene_DBO.early.case[,1:2])
mec_gene_DBO.early.case$prevalence = mec_gene_DBO.early.case$Yes/mec_gene_DBO.early.case$sample_count*100
mec_gene_DBO.early.case$NEC_status = "case"

data_mec.early.control = filter(data_mec, NEC_onset2=="Early" & NEC_status=="control")
mec_gene_DBO.early.control = as.data.frame.array(table(data_mec.early.control$Days_before_NEC_onset, data_mec.early.control$virMec_gene))
mec_gene_DBO.early.control$DBO = rownames(mec_gene_DBO.early.control)
mec_gene_DBO.early.control$sample_count = rowSums(mec_gene_DBO.early.control[,1:2])
mec_gene_DBO.early.control$prevalence = mec_gene_DBO.early.control$Yes/mec_gene_DBO.early.control$sample_count*100
mec_gene_DBO.early.control$NEC_status = "control"

mec_gene_DBO.early = rbind(mec_gene_DBO.early.case, mec_gene_DBO.early.control)
mec_gene_DBO.early$DBO = as.numeric(mec_gene_DBO.early$DBO)
mec_gene_DBO.early$NEC_onset = "Early"

##----------------------------------------------------------------------------------------------------------------------------------------------------------------
mec_gene_DBO.early$DBO_group[mec_gene_DBO.early$DBO>16] = "16-"
mec_gene_DBO.early$DBO_group[mec_gene_DBO.early$DBO<=16] = "09-16"
mec_gene_DBO.early$DBO_group[mec_gene_DBO.early$DBO<=8] = "01-08"
mec_gene_DBO.early$DBO_group[mec_gene_DBO.early$DBO==0] = "0"

plot_box <- ggplot(mec_gene_DBO.early, aes(x=DBO_group, y=prevalence, fill=NEC_status)) +
  geom_boxplot(outlier.size=0.05) +
  #scale_fill_brewer(palette="RdBu") +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  # geom_jitter(color="gray", size=0.8, alpha=0.9) +
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.75), size=3)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),legend.position="none") +
  ylab("Prevalence [% of cohort]") +
  xlab("Days before NEC onset") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # scale_x_discrete(limits = c("21-","11-20","1-10","0")) +
  scale_x_discrete(limits = c("16-","09-16","01-08","0")) +
  # scale_x_discrete(limits = c("21-","08-14","01-07","0")) +
  # geom_signif(comparisons = list(c("case", "control")),
  #             map_signif_level = F, textsize=5)
  #stat_compare_means() +
  # ylim(0,0.5) +
  stat_compare_means(method = "wilcox.test", label ="p.signif", label.x = 2, label.y.npc=0.9,  hide.ns = TRUE, size = 7)
plot_box
ggsave("250610_box_compare_mec_prevalence_NECstatus_early_DBO_v3.pdf", plot_box, width=4, height=4)

