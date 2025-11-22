library("doParallel")
library(ggplot2)
library(reshape)
library(labdsv)
library(vegan)
library(ggpubr)
library(permute)
library(dplyr)
library(rsample)
library(purrr)
library(tidyr)
library(ggpubr)
library(readxl)
set.seed(129)

###### Load and edit DF #######
NEC_DF<-read.csv('Phage_relAb_vOTU.csv',
                 header = T,
                 row.names = 1)
NEC_DF = subset(NEC_DF, select = -c(NEC_onset2, NEC_status, DOL, Days_before_NEC_onset, Patient))
k <- rowSums(NEC_DF) > 0
NEC_DF <- NEC_DF[k,]

metadata<-read.delim('nec_metadata.txt', header = TRUE)
row.names(metadata)<-metadata[,1]
metadata<-metadata[,-1]
NEC_DF_meta = merge(NEC_DF,metadata, by = "row.names")
row.names(NEC_DF_meta) <-NEC_DF_meta$Row.names
metadata = NEC_DF_meta[,(ncol(NEC_DF_meta)-83+2):ncol(NEC_DF_meta)] 
metadata = select(metadata, where(function(x) all(!is.na(x))))
metadata = select(metadata, -Specimen_ID, -DOB, -Date_of_isolation, -Rankorder, -Dummy, ABX_recent_abx,
                  -ABX_recent_class, -Other_ANTIM_recent_antim, -Other_ANTIM_recent_antim_class, -mat_inpatient_abx,
                  -mat_inpatient_abx_class)

NEC_DF = NEC_DF_meta[,2:(ncol(NEC_DF_meta)-83+1)]

NEC_DF_BC <-vegdist(NEC_DF, method = "bray")

subject<-metadata$Patient

subject_data<-read.delim("subject_Birthweight.txt", header = TRUE,check.names = FALSE, row.names=1)

##--------------------------------------------------------------------------------------
source("repeated_measures_PERMANOVA.R")
df_abx = read_excel("sample_abx_diet.xlsx")
sample_row_name = df_abx$DNA_Sample_ID
df_abx = df_abx[,-1]
row.names(df_abx) = sample_row_name

for (i in 1:dim(df_abx)[2]){
  sample_data<- df_abx[,i]
  row.names(sample_data) = sample_row_name
  subject_data = filter(subject_data, row.names(subject_data) %in% subject)
  sample_data = filter(sample_data, row.names(sample_data) %in% row.names(NEC_DF))
  sample_data$names = rownames(sample_data)
  library(dplyr)
  target <- rownames(as.matrix(NEC_DF_BC)) 
  sample_data_tmp = sample_data %>% arrange(factor(names, levels = target))
  row.names(sample_data_tmp) = sample_data_tmp$names
  sample_data = select(sample_data_tmp, -names)
  row.names(sample_data) = sample_data_tmp$names
  
  #Sample_data
  PERMANOVA_result = 
    PERMANOVA_repeat_measures(
      NEC_DF_BC,
      subject,
      sample_data = sample_data,
      metadata_order = c(names(sample_data)),
      permutations=999, ncores=18)
  
  write.table(PERMANOVA_result[1,], "PERMANOVA_result_phage_vOTU_abx_diet.csv", 
              append = TRUE, col.names = FALSE, sep = " ")
}

