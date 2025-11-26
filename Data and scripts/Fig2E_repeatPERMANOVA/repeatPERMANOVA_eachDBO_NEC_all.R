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
set.seed(56)

###### PERMANOVA function #######
source("repeated_measures_PERMANOVA_bary.R")

input_name = c('Phage_relAbundance_host_genus.csv',
               'NEC_bacteria_genus_RA.csv',
               'Phage_relAbundance_species_all.csv',
               'NEC_bacteria_species_RA.csv',
               'Phage_relAb_vOTU.csv')

output_group = c('Phage_Host_Genus',
                'Bac_Genus_meta4',
                'Phage_Host_Species',
                'Bac_Species_meta4',
                'Phage_vOTU')

###### Write file with column names #######
write.table(cbind("Dataset","Days_before_onset","sample#","infant#","Df","SumOfSqs","R2","F","Pr(>F)"), 
            "result/PERMANOVA_result_eachDBO_NEC_all.csv", 
            col.names = FALSE, row.names = FALSE, sep = ",")

for (i in 1:5){
  print(i)
  ###### Loop with DBO #######
  for (DBO in 0:14){ # DBO: Day_before_onset 
    
    ###### Load and edit DF #######
   
    NEC_DF<-read.csv(input_name[i],
                     header = T,
                     row.names = 1)
    
    NEC_DF = subset(NEC_DF, select = -c(NEC_onset2, NEC_status, DOL, Days_before_NEC_onset, Patient))
        k <- rowSums(NEC_DF) > 0
    NEC_DF <- NEC_DF[k,]
    
    ###### Load metadata and merge with DF #######
    metadata<-read.delim('nec_metadata.txt', header = TRUE)
    row.names(metadata)<-metadata[,1]
    metadata<-metadata[,-1]
    metadata = filter(metadata,NEC_onset-DOL==DBO)
    NEC_DF_meta = merge(NEC_DF,metadata, by = "row.names")
    row.names(NEC_DF_meta) <-NEC_DF_meta$Row.names
    metadata = NEC_DF_meta[,(ncol(NEC_DF_meta)-83+2):ncol(NEC_DF_meta)] 
    metadata = select(metadata, where(function(x) all(!is.na(x))))
    metadata = select(metadata, -Specimen_ID, -DOB, -Date_of_isolation, -Rankorder, -Dummy, ABX_recent_abx,
                      -ABX_recent_class, -Other_ANTIM_recent_antim, -Other_ANTIM_recent_antim_class, -mat_inpatient_abx,
                      -mat_inpatient_abx_class)
    NEC_DF = NEC_DF_meta[,2:(ncol(NEC_DF_meta)-83+1)]
    
    ###### BC calculation #######
    NEC_DF_BC <-vegdist(NEC_DF, method = "bray")
    
    ###### Assign subject as patient ID and subject data as NEC status #######
    subject<-metadata$Patient
    subject_data<-read.delim("subject_NEC.txt", header = TRUE,check.names = FALSE, row.names=1)
    subject_data = filter(subject_data, row.names(subject_data) %in% subject)
    
    ###### PERMOVONA test with subject data #######
    PERMANOVA_result =
      PERMANOVA_repeat_measures(
      NEC_DF_BC,
      subject, subject_data = subject_data,
      metadata_order = c(names(subject_data)),
      permutations=999, ncores=6)
    
    write.table(c(output_group[i],DBO,dim(NEC_DF)[1],dim(subject_data)[1],PERMANOVA_result[1,]), "result/PERMANOVA_result_eachDBO_NEC_all.csv", 
                append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
  }
}
