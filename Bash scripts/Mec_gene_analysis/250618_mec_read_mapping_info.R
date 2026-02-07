library(readxl)

setwd("~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 5/9ImpVir/mec_gene")

data_mec_ref = read_excel("250615_phage_mec_protein_189.xlsx", sheet = "Sheet4")

data_sample_DNA = read_excel("250615_phage_mec_protein_189.xlsx", sheet = "DNA")

data = merge(data_sample_DNA, data_mec_ref, by = "Patient")

write.csv(data, "250618_mec_DNA_mapping_info_3035.csv", row.names = FALSE)



data_mec_ref = read_excel("250615_phage_mec_protein_189.xlsx", sheet = "Sheet4")

data_sample_RNA = read_excel("250615_phage_mec_protein_189.xlsx", sheet = "RNA")

data = merge(data_sample_RNA, data_mec_ref, by = "Patient")

write.csv(data, "250618_mec_RNA_mapping_info_2911.csv", row.names = FALSE)
