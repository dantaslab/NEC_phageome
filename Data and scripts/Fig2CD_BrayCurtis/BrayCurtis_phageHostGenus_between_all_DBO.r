library("compositions")
library(readxl)
library(xlsx)
library(dplyr)
library("vegan")
library("Matrix")
library("phyloseq")
library(data.table)
library("reshape2")
library("ggplot2")
set.seed(129)

SM_DF_all = read.csv('Phage_relAbundance_host_genus.csv',
                     header = T)

k <- rowSums(SM_DF_all[,2:154]) > 0
SM_DF_all <- SM_DF_all[k,]
SM_DF_all[,2:154] = SM_DF_all[,2:154]/rowSums(SM_DF_all[,2:154]) * 100

taxmat <- as.data.frame(colnames(SM_DF_all))
colnames(taxmat)[1] <- "phage_hostGenus"

RowName <- SM_DF_all$SampleID
metadata <- as.data.frame(SM_DF_all$Days_before_NEC_onset)
rownames(metadata) = RowName
colnames(metadata)[1] <- "Group"
SM_DF_all <- select(SM_DF_all, -SampleID, -NEC_onset2, -NEC_status, -DOL, -Days_before_NEC_onset, -Patient)
rownames(SM_DF_all) = RowName
t_SM_DF_all = t(SM_DF_all)

setdiff(rownames(t_SM_DF_all),taxmat$phage_hostGenus)

rownames(taxmat) <- taxmat$phage_hostGenus
taxmat =  as.matrix(taxmat)

OTU = otu_table(t_SM_DF_all, taxa_are_rows = TRUE)
sampledata = sample_data(metadata)
TAX = tax_table(taxmat)
physeq1 = phyloseq(OTU, TAX, sampledata)

physeq1
head(otu_table(physeq1)[,1:6])

######################

abrel_bray <- phyloseq::distance(physeq1, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]

all_BC <- melt(abrel_bray)
all_BC_pair <- all_BC[all_BC$Var1!=all_BC$Var2,]

SM_sample <- read.csv('Phage_relAb_vOTU.csv',
                  header = T)
A = merge(all_BC_pair, SM_sample[,c(1,1718:1722)], by.x = "Var1", by.y = "SampleID")
colnames(A)[c(3:8)] = c("BrayCurtis_dissimilarity", "NEC_onset_patient1", "NEC_status_patient1", "DOL_patient1", "DBO_patient1", "Patient1")
B = merge(A, SM_sample[,c(1,1718:1722)], by.x = "Var2", by.y = "SampleID")
colnames(B)[c(9:13)] = c("NEC_onset_patient2", "NEC_status_patient2", "DOL_patient2", "DBO_patient2", "Patient2")
B$DOL_diff = abs(B$DOL_patient2-B$DOL_patient1)
B$BC_pair_group <- ifelse(B$Patient1==B$Patient2, 'Within', 'Between')
head(B)
BC_between_case_control_DBO = filter(B, BC_pair_group=="Between" & DBO_patient1==DBO_patient2)

write.csv(BC_between_case_control_DBO,"BC_phyloseq_hostGenus_between_all_DBO.csv")

