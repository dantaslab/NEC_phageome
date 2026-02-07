library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')

# Inputs are the results from Machine learning > Fig4 > Feature selection > MaAsLin2
test = c("0_1", "1_2", "2_3", "3_4", "4_5", "5_6", "6_7", "7_8", "8_9")
data_all_feature = c()
for (DBO in 0:8){

data_DBO_genus = read.delim(paste("250416_Phage_HostGenus_drop_LM_LOG_0.01_early_DBO", test[DBO+1], "/all_results.tsv", sep = ""),
                            header = TRUE, sep="\t")
data_DBO_genus = filter(data_DBO_genus, metadata == "NEC_status")
data_DBO_genus$Days_before_NEC_onset = DBO

data_DBO_species = read.delim(paste("250416_Phage_HostSpecies_drop_LM_LOG_0.001_early_DBO", test[DBO+1], "/all_results.tsv", sep = ""),
                              header = TRUE, sep="\t")
data_DBO_species = filter(data_DBO_species, metadata == "NEC_status")
data_DBO_species$Days_before_NEC_onset = DBO

data_all_feature = rbind(data_all_feature, data_DBO_genus, data_DBO_species)
}
##------------------------------------------------------------------------------------------------------------------------------
data_all_feature$onset_group = "Early"
data = data_all_feature
data = filter(data, feature =='Klebsiella' | feature =='Enterobacter' | feature =='Enterococcus_faecalis' | feature =='Serratia' |
                feature =='Dickeya' | feature =='Salmonella' | feature =='Shigella' |feature =='Klebsiella_grimontii' | feature =='Atlantibacter') # 9 selected viral features

# ref is control
data$sig_label[data$qval<0.25 & data$coef>0] = "all_control"
data$sig_label[data$qval<0.25 & data$coef<0] = "all_case"
data$sig_label[data$qval>=0.25] = "not_significant"

##----------------------------------------------------------------------------------------------------------------------------------------------------------------

plot_coef <- ggplot(data, aes(x=feature, y=-coef, color=sig_label)) +
  geom_errorbar(aes(ymin=-coef-stderr, ymax=-coef+stderr, group=onset_group), width=.2, color="black", position=position_dodge(0.8)) + 
  geom_jitter(
    aes(shape = onset_group),
    position = position_jitterdodge(jitter.width = 0),
    size = 3, color = "black"
  ) +
  geom_jitter(
    aes(shape = onset_group, color = sig_label),
    position = position_jitterdodge(jitter.width = 0),
    size = 2.5
  ) +
  # geom_point(size = 3, aes(shape = onset_group)) +
  # geom_point(aes(fill=sig_label,shape = onset_group), size = 3, color = "black")+
  # geom_point(aes(shape = onset_group),size = 2.4) +
  coord_flip() +
  scale_color_manual(values=c("all_case"="#f28a5f","not_significant"="#DDDDDD","all_control"="#3291a4"))+
  ylab("Coefficient") +
  xlab("Phage host") +
  # scale_color_discrete(name="NEC status") +
  theme(legend.position="none", legend.title = element_text(size=10), legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) +
  theme(axis.text.y=element_text(face="italic")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(limits=c('Atlantibacter','Klebsiella_grimontii','Shigella','Salmonella','Dickeya','Serratia','Enterococcus_faecalis','Enterobacter','Klebsiella')) +
  geom_hline(yintercept=-0, linetype="dashed", color = "black", size = 0.5) +
  facet_grid(onset_group~Days_before_NEC_onset, scale="fixed")
plot_coef
ggsave("250417_Plot_masslin_phageHost_early_coef_RA.pdf", plot_coef, width=14, height=3.1)

