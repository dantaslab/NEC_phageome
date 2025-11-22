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
source("repeated_measures_PERMANOVA.R")
df_subject = read_excel("subject_ALL.xlsx")
subject_row_name = df_subject$Patient
df_subject = df_subject[,-1]
row.names(df_subject) = subject_row_name

for (i in 13:dim(df_subject)[2]){
  subject_data<- df_subject[,i]
  row.names(subject_data) =  subject_row_name
  na_patient = subject_row_name[subject_data[,1]=="NA"]
  subject_data = filter(subject_data, subject_data[,1]!="NA")
  
  ##########################################
  
  NEC_DF<-read.csv('Phage_relAb_vOTU.csv',
                   header = T,
                   row.names = 1)
  
  NEC_DF = filter(NEC_DF, !grepl(paste(na_patient, collapse="|"), Patient))
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
  subject_data = filter(subject_data, row.names(subject_data) %in% subject)
  
  ##########################################
  #subject_data
  PERMANOVA_result = 
    PERMANOVA_repeat_measures(
      NEC_DF_BC,
      subject, subject_data = subject_data,
      metadata_order = c(names(subject_data)),
      permutations=999, ncores=18)
  
  write.table(PERMANOVA_result[1,], "PERMANOVA_result_vOTU_base_mat_NA.csv", 
              append = TRUE, col.names = FALSE, sep = " ")
}

