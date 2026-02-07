library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')

setwd("~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 3/Fig3A_extension/Fig3A_eachDBO")
set.seed(129)

test = c("0_1", "1_2", "2_3", "3_4", "4_5", "5_6", "6_7", "7_8", "8_9")
data_all_feature = c()
for (i in 1:9){
# DBO_larger1 sample
data_DBO_genus = read.delim(paste("MaAsLin2_Bac_Resistome_RA/250416_Bac_Resisome_class_drop_LM_LOG_0_late_DBO", test[i], "/significant_results.tsv", sep = ""),
                      header = TRUE, sep="\t")
data_DBO_genus = filter(data_DBO_genus, metadata == "NEC_status")
if(dim(data_DBO_genus)[1]>0){data_DBO_genus$Days_before_NEC_onset = test[i]}

data_DBO_species = read.delim(paste("MaAsLin2_Bac_Resistome_RA/250416_Bac_Resisome_gene_drop_LM_LOG_0_late_DBO", test[i], "/significant_results.tsv", sep = ""),
                            header = TRUE, sep="\t")
data_DBO_species = filter(data_DBO_species, metadata == "NEC_status")
data_DBO_species$Days_before_NEC_onset = test[i]

data_all_feature = rbind(data_all_feature, data_DBO_genus, data_DBO_species)
}

write_csv(data_all_feature,"250416_maaslin2_resistome_sig_features_late.csv")


