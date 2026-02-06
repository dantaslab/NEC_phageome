# RANDOM FOREST CLASSIFIER ANALYSES ############################################
# For those only interested in this analysis: 
#save(list=c('taxa.filt'), file='RData/ACS_AD/230109_STMRevision/230111_bactaxafilt.RData')

## Load packages and data-------------------------------------------------------
# Packages specific to these analyses: 
library(dplyr)     #v1.0.8    
library(tidyr)     #v1.2.0
library(ggplot2)   #v3.3.5 
#library(phyloseq)  #v1.38.0
library(Boruta)    #v7.0.0
library(caret)     #v6.0.86
library(VIM)       #v6.1.0
library(stringr)   #v1.4.0
library(egg)       #v0.4.5
library(scales)    #v1.1.1
library(readxl)
library(ggpubr) 
library(forcats)

## Load data----------------------------------------------------
#setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta')

## Phage host data
dataset_all_genus_phage = read.csv('Input/Phage_relAbundance_genus_all.csv', header = T, row.names = 1)
dataset_all_genus_phage = select(dataset_all_genus_phage, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_genus_phage = dataset_all_genus_phage*100
dataset_all_genus_phage[dataset_all_genus_phage < 0.01] <- 0
dataset_all_genus_phage = dataset_all_genus_phage[rowSums(dataset_all_genus_phage) > 0, ]
dataset_all_genus_phage <- dataset_all_genus_phage/rowSums(dataset_all_genus_phage) * 100
colnames(dataset_all_genus_phage)<-paste(colnames(dataset_all_genus_phage),"vir",sep="_")
dataset_all_species_phage = read.csv('Input/Phage_relAbundance_species_all.csv', header = T, row.names = 1)
dataset_all_species_phage = select(dataset_all_species_phage, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_species_phage = dataset_all_species_phage*100
dataset_all_species_phage[dataset_all_species_phage < 0.001] <- 0
dataset_all_species_phage = dataset_all_species_phage[rowSums(dataset_all_species_phage) > 0, ]
dataset_all_species_phage <- dataset_all_species_phage/rowSums(dataset_all_species_phage) * 100
colnames(dataset_all_species_phage)<-paste(colnames(dataset_all_species_phage),"vir",sep="_")
dataset_all_phage = merge(dataset_all_genus_phage,dataset_all_species_phage,by="row.names")
rownames(dataset_all_phage) = dataset_all_phage$Row.names
dataset_all_phage$Row.names = NULL

## Bacterial function data
dataset_all_resistome_class = read.csv('Input/241104_NEC_ShortBRED_byclass.csv', header = T, row.names = 1) #241104_NEC_ShortBRED_byclass.csv
dataset_all_resistome_class = select(dataset_all_resistome_class, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_resistome_class = dataset_all_resistome_class[rowSums(dataset_all_resistome_class) > 0, ]
dataset_all_resistome_class <- dataset_all_resistome_class/rowSums(dataset_all_resistome_class) * 100
dataset_all_resistome_gene = read.csv('Input/240831_NEC_ShortBRED_V2.csv', header = T, row.names = 1) #241104_NEC_ShortBRED_byclass.csv
dataset_all_resistome_gene = select(dataset_all_resistome_gene, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_resistome_gene = dataset_all_resistome_gene[rowSums(dataset_all_resistome_gene) > 0, ]
dataset_all_resistome_gene <- dataset_all_resistome_gene/rowSums(dataset_all_resistome_gene) * 100
dataset_all_resistome = merge(dataset_all_resistome_class,dataset_all_resistome_gene,by="row.names")
rownames(dataset_all_resistome) = dataset_all_resistome$Row.names
dataset_all_resistome$Row.names = NULL
colnames(dataset_all_resistome)<-paste(colnames(dataset_all_resistome),"rest",sep="_")

# dataset_all_virulence = read.csv('mmc17_MAG_virulence.csv', header = T, row.names = 1)
# dataset_all_virulence = select(dataset_all_virulence, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
# dataset_all_virulence = dataset_all_virulence*100

dataset_all_pathway = read.csv('Input/mmc12_MAG_Pathway_RA.csv', header = T, row.names = 1)
dataset_all_pathway = select(dataset_all_pathway, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_pathway = dataset_all_pathway*100
colnames(dataset_all_pathway)<-paste(colnames(dataset_all_pathway),"path",sep="_")

dataset_all_iRep_genus = read.csv('Input/mmc15_iRep_genus.csv', header = T, row.names = 1)
dataset_all_iRep_genus = select(dataset_all_iRep_genus, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
colnames(dataset_all_iRep_genus)<-paste(colnames(dataset_all_iRep_genus),"iRep",sep="_")
dataset_all_iRep_species = read.csv('Input/mmc15_iRep_species.csv', header = T, row.names = 1)
dataset_all_iRep_species = select(dataset_all_iRep_species, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
colnames(dataset_all_iRep_species)<-paste(colnames(dataset_all_iRep_species),"iRep",sep="_")
dataset_all_iRep = merge(dataset_all_iRep_genus,dataset_all_iRep_species,by="row.names")
rownames(dataset_all_iRep) = dataset_all_iRep$Row.names
dataset_all_iRep$Row.names = NULL
dataset_imp_Virtaxairep  = merge(dataset_all_resistome,dataset_all_iRep,by="row.names")
rownames(dataset_imp_Virtaxairep) = dataset_imp_Virtaxairep$Row.names
dataset_imp_Virtaxairep$Row.names = NULL

## MAX data
dataset_all_MAX_taxa = read.csv('Input/mmc13_MAX_Species_RA.csv', header = T, row.names = 1)
dataset_all_MAX_taxa = select(dataset_all_MAX_taxa, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_MAX_taxa = dataset_all_MAX_taxa[rowSums(dataset_all_MAX_taxa) > 0, ]
# dataset_all_MAX_taxa <- dataset_all_MAX_taxa/rowSums(dataset_all_MAX_taxa) * 100
dataset_all_MAX_path = read.csv('Input/mmc14_MAX_Pathway_RA.csv', header = T, row.names = 1) #241104_NEC_ShortBRED_byclass.csv
dataset_all_MAX_path = select(dataset_all_MAX_path, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_MAX_path = dataset_all_MAX_path[rowSums(dataset_all_MAX_path) > 0, ]
# dataset_all_MAX_path <- dataset_all_MAX_path/rowSums(dataset_all_MAX_path) * 100
dataset_all_MAX = merge(dataset_all_MAX_taxa,dataset_all_MAX_path,by="row.names")
rownames(dataset_all_MAX) = dataset_all_MAX$Row.names
dataset_all_MAX$Row.names = NULL
colnames(dataset_all_MAX)<-paste(colnames(dataset_all_MAX),"MAX",sep="_")
dataset_imp_VirtaxairepMAX = merge(dataset_imp_Virtaxairep,dataset_all_MAX,by="row.names")


## Bactieral taxonomy data
dataset_all_genus_bac = read.csv('Input/240510_NEC_bacteria_genus_RA_meta4.csv', header = T, row.names = 1)
dataset_all_genus_bac = select(dataset_all_genus_bac, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_genus_bac[dataset_all_genus_bac < 0.01] <- 0
dataset_all_genus_bac = dataset_all_genus_bac[rowSums(dataset_all_genus_bac) > 0, ]
dataset_all_genus_bac <- dataset_all_genus_bac/rowSums(dataset_all_genus_bac) * 100
colnames(dataset_all_genus_bac)<-paste(colnames(dataset_all_genus_bac),"bac",sep="_")
dataset_all_species_bac = read.csv('Input/240510_NEC_bacteria_species_RA_meta4.csv', header = T, row.names = 1)
dataset_all_meta = select(dataset_all_species_bac, Days_before_NEC_onset,Patient,NEC_status,NEC_onset2)
dataset_all_species_bac = select(dataset_all_species_bac, -NEC_onset2,-NEC_status,-DOL,-Days_before_NEC_onset,-Patient)
dataset_all_species_bac[dataset_all_species_bac < 0.001] <- 0
Index = rowSums(dataset_all_species_bac) > 0
dataset_all_species_bac = dataset_all_species_bac[Index, ]
dataset_all_species_bac <- dataset_all_species_bac/rowSums(dataset_all_species_bac) * 100
colnames(dataset_all_species_bac)<-paste(colnames(dataset_all_species_bac),"bac",sep="_")
dataset_all_species_bac = merge(dataset_all_species_bac,dataset_all_meta[Index, ],by="row.names")
rownames(dataset_all_species_bac) = dataset_all_species_bac$Row.names
dataset_all_species_bac$Row.names = NULL
dataset_all_bac = merge(dataset_all_genus_bac,dataset_all_species_bac,by="row.names")
row.names(dataset_all_bac) = dataset_all_bac$Row.names
dataset_all_bac = dataset_all_bac[,-1]

# Combine data
dataset_all1 = merge(dataset_all_phage,dataset_all_bac,by="row.names")
row.names(dataset_all1) = dataset_all1$Row.names
dataset_all1$Row.names = NULL
dataset_all2 = merge(dataset_all1,dataset_all_pathway,by="row.names")
dataset_all  = merge(dataset_all2,dataset_imp_VirtaxairepMAX,by="Row.names")
row.names(dataset_all) = dataset_all$Row.names
dataset_all = filter(dataset_all, !grepl("204-01|2156-01|415-01|2132-01|2234-01|234-01|
                                   421-01|2218-01|2238-02|49-01|257-01|2022-01|
                                   2077-01|2190-01|2228-01", Patient))
dataset_all = filter(dataset_all, Patient!="421-01" & Patient!="2077-01")
dataset_all = filter(dataset_all, NEC_onset2=="Early")

# Access clinical metadata
metadata_sample1 = read_excel("Input/240404_sample_ALL.xlsx")
metadata_sample2 = read_excel("Input/241103_sample_medications.xlsx")
metadata_sample  = merge(metadata_sample1, metadata_sample2, by="DNA_Sample_ID")
metadata_subject = read_excel("Input/241104_subject_ALL_noMat.xlsx")
metadata         = merge(metadata_sample, metadata_subject, by="Patient")
# metadata_stage   = read.delim('subject_nec_stage2.txt', header = TRUE)
# metadata         = merge(metadata, metadata_stage, by="Patient")
row.names(metadata) = metadata$DNA_Sample_ID
metadata = filter(metadata, !grepl("204-01|2156-01|415-01|2132-01|2234-01|234-01|
                                 421-01|2218-01|2238-02|49-01|257-01|2022-01|
                                 2077-01|2190-01|2228-01", Patient))
metadata = filter(metadata, Patient!="421-01" & Patient!="2077-01")

ImpVar.df   = read.csv('Input/241128_ML_sallTaxa_MAX_VarImpTest2_early_v1.2.csv', header = T)
ImpVar.all  = unique(ImpVar.df$Feature)
ImpVar.vir  = unique(ImpVar.df$Feature[ImpVar.df$DataSet=="VirTaxa"])
# Boruta.df   = read.csv('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_Boruta/250415_NEC_Boruta_early_sum.csv', header = T)
Boruta.df   = read.csv('Input/250415_NEC_Boruta_early_sum.csv', header = T)
Boruta.vir  = unique(Boruta.df$Feature[Boruta.df$Dataset=="VirTaxa"])
# Maaslin.df  = read.csv('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 3/Fig3A_extension/Fig3A_eachDBO/250416_maaslin2_phage_host_sig_features_early.csv', header = T)
Maaslin.df  = read.csv('Input/250416_maaslin2_phage_host_sig_features_early.csv', header = T)
Maaslin.vir =  paste(unique(Maaslin.df$feature) ,"vir",sep="_")

write.table(cbind('Days_before_onset','Control.Control', 'Control.Case', 'Case.Control', 'Case.Case', 
                  'Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyNull', 'AccuracyPValue', 
                  'McnemarPValue', 'Sensitivity', 'Specificity', 'Pos.Pred.Value', 'Neg.Pred.Value', 
                  'Precision', 'Recall', 'F1', 'Prevalence', 'Detection.Rate', 'Detection.Prevalence', 
                  'Balanced.Accuracy', 'Data', 'DataSet', 'Taxa', 'Seed'), 
            "250416_ML_ImpVIr_early_pred_v1.2.csv", 
            col.names = FALSE, row.names = FALSE, sep = ",")

write.table(cbind('Days_before_onset','DataSet', 'Taxa', 'Feature', 'Feature.MeanImp', 'Feature.SDImp', 'Feature.NULLMeanImp', 'Feature.NULLSDImp'), 
            "250416_ML_ImpVIr_early_VarImpSumm_v1.2.csv", 
            col.names = FALSE, row.names = FALSE, sep = ",")

for (DBO in 0:14){
  GM_all_case = filter(dataset_all,  NEC_status=="case" & (Days_before_NEC_onset==DBO | Days_before_NEC_onset==DBO+1))
  GM_all_control= filter(dataset_all, NEC_status=="control" & Days_before_NEC_onset==DBO)
  GM_all = rbind(GM_all_case,GM_all_control)
  GM_all = GM_all[,-1]
  
  # Transform to counts, preserving precision to 0.00001%
  GM.all <- select(GM_all, -Days_before_NEC_onset, -Patient, -NEC_status, -NEC_onset2)
  GM.filt <- GM.all[, colSums(GM.all) > 0]
  GM.filt.int <- round(GM.filt*10000000)
  GM.caret <- data.frame(GM.filt.int)
  
  meta.caret <- metadata # filter(metadata, Days_before_NEC_onset==DBO)

  # VIM also removes variable class information (factors, etc.). Fix here by
  # converting categorical variables to factors - v1
  factor.vars <- c("NEC_status",
                   "Recent_penicillin", "Recent_aminoglycoside", "Recent_lincosamide", "Recent_glycopeptide",
                   "Recent_carbapenem", "Recent_macrolide", "Recent_1st_gen_cephalosporin", "Recent_2nd_gen_cephalosporin",
                   "Recent_3rd_gen_cephalosporin", "Recent_4th_gen_cephalosporin",
                   "Other_ANTIM_recent", "feeding_current", "milk_current",
                   "Acetaminophen_recent", "Caffeine_recent", "Calcium_gluconate_recent", "Chlorathiazide_recent",
                   "Cholecalciferol_recent", "Dexamethasone_recent", "Dopamine_recent", "Dobutamine_recent",
                   "Famotidine_recent", "Fentanyl_recent", "Hydrocortisone_recent", "Indomethacin_recent", "Insulin_recent",
                   "Immunoglobulin_recent", "Iron_recent", "Lasix_recent", "Midazolam_recent", "Morphine_recent",
                   "Multi_vitamin_with_iron_recent", "Multi_vitamin_no_iron_recent", "Pancuronium_or_Vecurionium_recent",
                   "Phenobarbitol_recent", "Potassium_chloride_recent", "Sodium_chloride_recent", "Sodium_bicarbonate_recent",
                   "Surfactant_recent", "Vitamin_A_recent", "Sex", "route", "multiple", "Site")

  meta.caret[, factor.vars] <- lapply(meta.caret[, factor.vars], factor)

  GM.caret <- merge(GM.caret, meta.caret[21], by='row.names')
  rownames(GM.caret) <- GM.caret$Row.names
  GM.caret$Row.names <- NULL
  
  
  # Subset taxa abundance data to these all existing taxa (for entire cohort including
  # train and validation sets -- will be re-partitioned later using same index). 
  fs.virtaxa.all     = colnames(dataset_all_phage)
  fs.virtaxa.imp     = ImpVar.vir
  fs.virtaxa.boruta  = Boruta.vir
  fs.virtaxa.maaslin = Maaslin.vir
  fs.all.impvari1    = c("Klebsiella_vir", "Enterobacter_vir", "Cronobacter_vir", "Serratia_vir", "Dickeya_vir", "Bariatricus_vir", "Salmonella_vir", "Enterococcus_faecalis_vir")
  fs.all.impvari2    = c('Shigella_vir', 'Enterobacter_vir', 'Salmonella_vir', 'Klebsiella_vir', 'Klebsiella_grimontii_vir', 'Atlantibacter_vir', 'Dickeya_vir', 'Serratia_vir', 'Enterococcus_faecalis_vir')
  fs.all.impvari3    = c('Clostridioides_vir', 'Clostridium_vir', 'Cronobacter_vir', 'Enterobacter_vir', 'Klebsiella_vir', 'Dickeya_vir', 'Enterobacter_oligotrophica_vir', 'Klebsiella_grimontii_vir', 'Escherichia_vir', 'Salmonella_vir', 'Serratia_vir', 'Citrobacter_werkmanii_vir', 'Enterococcus_faecalis_vir', 'Escherichia_albertii_vir', 'Citrobacter_vir', 'Terrisporobacter_vir', 'Atlantibacter_vir', 'Shigella_vir')
  fs.all.impvari4    = c('Cronobacter_vir', 'Dickeya_vir', 'Enterobacter_vir', 'Enterobacter_oligotrophica_vir', 'Klebsiella_grimontii_vir', 'Escherichia_vir', 'Klebsiella_vir', 'Salmonella_vir', 'Serratia_vir', 'Citrobacter_werkmanii_vir', 'Enterococcus_faecalis_vir', 'Escherichia_albertii_vir', 'Citrobacter_vir', 'Clostridioides_vir', 'Atlantibacter_vir', 'Shigella_vir')
  fs.all.impvari5    = c('Dickeya_vir', 'Enterobacter_vir', 'Klebsiella_vir', 'Serratia_vir')
  
  virtaxa.all     = GM.caret[ , colnames(GM.caret) %in% c(fs.virtaxa.all, "NEC_status")] 
  virtaxa.imp     = GM.caret[ , colnames(GM.caret) %in% c(fs.virtaxa.imp, "NEC_status")] 
  virtaxa.boruta  = GM.caret[ , colnames(GM.caret) %in% c(fs.virtaxa.boruta, "NEC_status")] 
  virtaxa.maaslin = GM.caret[ , colnames(GM.caret) %in% c(fs.virtaxa.maaslin, "NEC_status")] 
  all.impvari1    = GM.caret[ , colnames(GM.caret) %in% c(fs.all.impvari1, "NEC_status")]
  all.impvari2    = GM.caret[ , colnames(GM.caret) %in% c(fs.all.impvari2, "NEC_status")]
  all.impvari3    = GM.caret[ , colnames(GM.caret) %in% c(fs.all.impvari3, "NEC_status")]
  all.impvari4    = GM.caret[ , colnames(GM.caret) %in% c(fs.all.impvari4, "NEC_status")]
  all.impvari5    = GM.caret[ , colnames(GM.caret) %in% c(fs.all.impvari5, "NEC_status")]
  
  ## Create data subsets (i.e. omit biomarker categories)-------------------------
  
  # FUNCTION: createSubsets - supply clinical metadata features to omit, then 
  # merge remaining clinical metadata features with previously selected taxonomic 
  # features.
  # Required packages: none
  # Arguments:
  #  RFmetadata (dataframe)    = base.df (imputed metadata)
  #  nullvars (chr list)       = list of vars to exclude from the model
  #  Taxa (num matrix)         = taxa.caret.boruta (abundances of selected taxa)
  #
  # Return:
  #  Named list: Base = values for selected metadata features, WithTaxa = values
  #  for selected metadata features as well as feature-selected taxa. This will
  #  enable testing of improvements in model performances with addition of
  #  taxonomic feature data.
  
  createSubsets <- function(RFmetadata, nullvars, Taxa){
    if (!all(nullvars %in% colnames(RFmetadata))){
      stop("At least one variable name not in dataframe")
    }
    
    #Non-taxonomic features only  
    metadata <- RFmetadata
    metadata[ , nullvars] <- list(NULL)
    
    #With taxonomic features
    metadata.wtax <- merge(metadata, Taxa, all = TRUE, by='row.names')
    rownames(metadata.wtax) <- metadata.wtax$Row.names
    metadata.wtax$Row.names <- NULL
    
    #RETURN
    all.out <- list("Base"= metadata, 
                    "WithTaxa"= metadata.wtax)
  }
  
  
  
  # Define metadata subsets (biomarker exclusion lists):
  
  # Main text models: 
  #  No exclusions
  all.nullvars <- list()
  
  # Create Data Subsets
  #  Main Text models:
  all.virtaxa     <- createSubsets(virtaxa.all, all.nullvars, Taxa=virtaxa.all)
  imp.virtaxa     <- createSubsets(virtaxa.imp, all.nullvars, Taxa=virtaxa.imp)
  boruta.virtaxa  <- createSubsets(virtaxa.boruta, all.nullvars, Taxa=virtaxa.boruta)
  maaslin.virtaxa <- createSubsets(virtaxa.maaslin, all.nullvars, Taxa=virtaxa.maaslin)
  imp.allvari1    <- createSubsets(all.impvari1, all.nullvars, Taxa=all.impvari1)
  imp.allvari2    <- createSubsets(all.impvari2, all.nullvars, Taxa=all.impvari2)
  imp.allvari3    <- createSubsets(all.impvari3, all.nullvars, Taxa=all.impvari3)
  imp.allvari4    <- createSubsets(all.impvari4, all.nullvars, Taxa=all.impvari4)
  imp.allvari5    <- createSubsets(all.impvari5, all.nullvars, Taxa=all.impvari5)
  
 
  ## Train Random Forest classifiers----------------------------------------------
  
  # Define the test harness. Within the training cohort, will train using 10-fold
  # cross-validation. 
  cv10 <- trainControl(method='cv', number=10, classProbs=T, savePredictions =T)
  
  
  # Define categorical variables that should not be normalized / scaled.
  varsnot2norm <- factor.vars
  
  
  # FUNCTION: train_rf_models trains a classifier on the provided training data  
  # and test harness. It does this iteratively (100x) on random 80:20 partitions of
  # the training cohort, testing on the entire validation cohort at each iteration. 
  # Predictive results are collated.
  # Required packages: caret
  # Arguments:
  #  data (dataframe)         = data subset generated by createSubsets()
  #  control.harness          = cv10. Control harness generated  by caret::trainControl()
  #  data.name.string (str)   = Identifier for the data subset provided
  #  varsNot2Norm (chr list)  = List of categorical variables that shouldn't be
  #                             normalized / scaled.
  #  shuffle.class (logical)  = TRUE if class labels should be shuffled during
  #                             model training to generate null performance 
  #                             parameter distributions.
  #
  # Return: 
  #  list(prediction results [on training cohort], prediction results [on 
  #            validation cohort], variable importances)
  
  train_rf_models <- function(data, control.harness, data.name.string, 
                              varsNot2Norm, shuffle.class) {
    #Cross Validation (within training cohort) 
    out <- list()
    varimportance <- list() 
    
    #Validation Set 
    out.val <- list()
    
    #Separate into train/test and validation subsets, using same index as before.
    set.seed(123)
    
    
    for (i in 1:100) {
      #Create random partition of training cohort (80:20). 
      set.seed(i)
      
      data1_idx <- createDataPartition(data$NEC_status, p = 0.75, list=FALSE)
      data1 <- data[data1_idx, ]
      data.val <- data[-data1_idx, ]
      
      # train_idx <- createDataPartition(data1$NEC_status, p = 0.8, list=FALSE)
      
      data.train <- data1 #[train_idx, ]
      data.train <- select(data.train, where(function(x) length(levels(factor(x)))!=1))
      # data.test <- data1[-train_idx, ]
      
      # Optionally shuffle class labels in training data (data.train).
      if (shuffle.class == TRUE) {
        data.train$NEC_status <- sample(data.train$NEC_status)
      }
      
      #Pre-process data (center and scale)
      preprocessrule <- preProcess(data.train[, !(colnames(data.train) %in% varsNot2Norm)], 
                                   method = c('center', 'scale'))
      data.train.p <- predict(preprocessrule, data.train)
      # data.test.p <- predict(preprocessrule, data.test)
      
      data.val.p <- predict(preprocessrule, data.val)
      
      
      #train model
      set.seed(123)
      fit.rf <- train(NEC_status~., data = data.train.p, method = 'rf',
                      metric = 'Accuracy', trControl = control.harness)
      
      #Make predictions for this iteration's test set and store performance measures.
      # predictions.rf <- predict(fit.rf, data.test.p)
      # cM <- confusionMatrix(predictions.rf, data.test.p$NEC_status)
      # pred.results <- c('Control-Control'=cM$table[1,1], 
      #                   'Control-Case'=cM$table[2,1], 
      #                   'Case-Control'=cM$table[1,2], 
      #                   'Case-Case'= cM$table[2,2], 
      #                   cM$overall, 
      #                   cM$byClass,
      #                   'Data' = data.name.string, 
      #                   'Seed' = i)
      # 
      # out[[i]] <- pred.results
      
      #Make predictions for retained VALIDATION set and store performance measures
      predictions.val.rf <- predict(fit.rf, data.val.p)
      cM.val <- confusionMatrix(predictions.val.rf, data.val.p$NEC_status)
      pred.results.val <- c('Control-Control'=cM.val$table[1,1], 
                            'Control-Case'=cM.val$table[2,1], 
                            'Case-Control'=cM.val$table[1,2], 
                            'Case-Case'= cM.val$table[2,2], 
                            cM.val$overall, 
                            cM.val$byClass,
                            'Data' = data.name.string, 
                            'Seed' = i)
      
      out.val[[i]] <- pred.results.val
      
      
      #Find important vars
      fit.rf.importance <- varImp(fit.rf, scale=FALSE)
      fit.rf.importance.df <- fit.rf.importance$importance
      var.importance <- fit.rf.importance.df$Overall
      names(var.importance) <- row.names(fit.rf.importance.df)
      
      varimportance[[i]] <- var.importance
    }
    
    #Return 
    # out.df <- data.frame(do.call('rbind', out))
    out.val.df <- data.frame(do.call('rbind', out.val))
    varimportance.df <- data.frame(do.call('rbind', varimportance))
    
    allout <- list(#'Pred.Results.CV'=out.df,
      'Pred.Results.Val'=out.val.df, 
      'Var.Importance'=varimportance.df)
  }
  
  
  
  # Train the models (without shuffling class labels)
  require(doParallel)
  cl = makeCluster(24)
  registerDoParallel(cl)
  
  # all (all biomarkers)

  rf.all_Virtaxa     <- train_rf_models(all.virtaxa$Base, cv10, 'VirTaxa', 
                                        varsnot2norm, FALSE)
  rf.imp_Virtaxa     <- train_rf_models(imp.virtaxa$Base, cv10, 'impVirTaxa', 
                                        varsnot2norm, FALSE)
  rf.boruta_Virtaxa  <- train_rf_models(boruta.virtaxa$Base, cv10, 'BorutaVirTaxa', 
                                        varsnot2norm, FALSE)
  rf.maaslin_Virtaxa <- train_rf_models(maaslin.virtaxa$Base, cv10, 'MaaslinVirTaxa', 
                                        varsnot2norm, FALSE)
  rf.imp_Allvari1    <- train_rf_models(imp.allvari1$Base, cv10, 'impAllVari8', 
                                        varsnot2norm, FALSE)
  rf.imp_Allvari2    <- train_rf_models(imp.allvari2$Base, cv10, 'impAllVari9', 
                                        varsnot2norm, FALSE)
  rf.imp_Allvari3    <- train_rf_models(imp.allvari3$Base, cv10, 'impAllVari18', 
                                        varsnot2norm, FALSE)
  rf.imp_Allvari4    <- train_rf_models(imp.allvari4$Base, cv10, 'impAllVari16', 
                                        varsnot2norm, FALSE)
  rf.imp_Allvari5    <- train_rf_models(imp.allvari5$Base, cv10, 'impAllVari4', 
                                        varsnot2norm, FALSE)

  # # Beta (Excluding amyloid)
  # rf.beta <- train_rf_models(beta.data$Base, cv10, 'Beta', 
  #                            varsnot2norm, FALSE)
  # rf.beta_T <- train_rf_models(beta.data$WithTaxa, cv10, 'Beta plus Taxa', 
  #                              varsnot2norm, FALSE)
  
  # Collate results: 
  
  all.pred.Val <- rbind(rf.all_Virtaxa$Pred.Results.Val,
                        rf.imp_Virtaxa$Pred.Results.Val,
                        rf.boruta_Virtaxa$Pred.Results.Val,
                        rf.maaslin_Virtaxa$Pred.Results.Val,
                        rf.imp_Allvari1$Pred.Results.Val,
                        rf.imp_Allvari2$Pred.Results.Val,
                        rf.imp_Allvari3$Pred.Results.Val,
                        rf.imp_Allvari4$Pred.Results.Val,
                        rf.imp_Allvari5$Pred.Results.Val)
  
  rf.all_Virtaxa$Var.Importance$Dataset     = "virtaxa"
  rf.imp_Virtaxa$Var.Importance$Dataset     = "impvirtaxa"
  rf.boruta_Virtaxa$Var.Importance$Dataset  = "impvirboruta"
  rf.maaslin_Virtaxa$Var.Importance$Dataset = "impvirmaaslin"
  rf.imp_Allvari1$Var.Importance$Dataset    = "impallvari8"
  rf.imp_Allvari2$Var.Importance$Dataset    = "impallvari9"
  rf.imp_Allvari3$Var.Importance$Dataset    = "impallvari18"
  rf.imp_Allvari4$Var.Importance$Dataset    = "impallvari16"
  rf.imp_Allvari5$Var.Importance$Dataset    = "impallvari4"
  
  all.VarImp.list <- list(rf.all_Virtaxa$Var.Importance,
                          rf.imp_Virtaxa$Var.Importance,
                          rf.boruta_Virtaxa$Var.Importance,
                          rf.maaslin_Virtaxa$Var.Importance,
                          rf.imp_Allvari1$Var.Importance,
                          rf.imp_Allvari2$Var.Importance,
                          rf.imp_Allvari3$Var.Importance,
                          rf.imp_Allvari4$Var.Importance,
                          rf.imp_Allvari5$Var.Importance)
  
  all.VarImp <- all.VarImp.list %>% Reduce(function(d1, d2) full_join(d1, d2), .)
  all.VarImp$Dataset = NULL
  
  ## Shuffle class labels and retrain for null distributions----------------------
  
  # We'd like to compare importance of individual features to the performance of 
  # the predictive models against their importance when class labels have been
  # shuffled during model training (to generate random null distributions for 
  # the variable importances). 
  
  # Alpha (all biomarkers)

  rfN.all_Virtaxa     <- train_rf_models(all.virtaxa$Base, cv10, 'VirTaxa', 
                                         varsnot2norm, TRUE)
  rfN.imp_Virtaxa     <- train_rf_models(imp.virtaxa$Base, cv10, 'impVirTaxa', 
                                         varsnot2norm, TRUE)
  rfN.boruta_Virtaxa  <- train_rf_models(boruta.virtaxa$Base, cv10, 'BorutaVirTaxa', 
                                         varsnot2norm, TRUE)
  rfN.maaslin_Virtaxa <- train_rf_models(maaslin.virtaxa$Base, cv10, 'MaaslinVirTaxa', 
                                         varsnot2norm, TRUE)
  rfN.imp_Allvari1    <- train_rf_models(imp.allvari1$Base, cv10, 'impAllVari8', 
                                         varsnot2norm, TRUE)
  rfN.imp_Allvari2    <- train_rf_models(imp.allvari2$Base, cv10, 'impAllVari9', 
                                         varsnot2norm, TRUE)
  rfN.imp_Allvari3    <- train_rf_models(imp.allvari3$Base, cv10, 'impAllVari18', 
                                         varsnot2norm, TRUE)
  rfN.imp_Allvari4    <- train_rf_models(imp.allvari4$Base, cv10, 'impAllVari16', 
                                         varsnot2norm, TRUE)
  rfN.imp_Allvari5    <- train_rf_models(imp.allvari5$Base, cv10, 'impAllVari4', 
                                         varsnot2norm, TRUE)
  
  # # Beta (Excluding amyloid)
  # rfN.beta <- train_rf_models(beta.data$Base, cv10, 'Beta', 
  #                            varsnot2norm, FALSE)
  # rfN.beta_T <- train_rf_models(beta.data$WithTaxa, cv10, 'Beta plus Taxa', 
  #                              varsnot2norm, FALSE)
  
  # Collate results: 
  
  NULL.all.pred.Val <- rbind(rfN.all_Virtaxa$Pred.Results.Val,
                             rfN.imp_Virtaxa$Pred.Results.Val,
                             rfN.boruta_Virtaxa$Pred.Results.Val,
                             rfN.maaslin_Virtaxa$Pred.Results.Val,
                             rfN.imp_Allvari1$Pred.Results.Val,
                             rfN.imp_Allvari2$Pred.Results.Val,
                             rfN.imp_Allvari3$Pred.Results.Val,
                             rfN.imp_Allvari4$Pred.Results.Val,
                             rfN.imp_Allvari5$Pred.Results.Val)
  
  rfN.all_Virtaxa$Var.Importance$Dataset     = "virtaxa"
  rfN.imp_Virtaxa$Var.Importance$Dataset     = "impvirtaxa"
  rfN.boruta_Virtaxa$Var.Importance$Dataset  = "impvirboruta"
  rfN.maaslin_Virtaxa$Var.Importance$Dataset = "impvirmaaslin"
  rfN.imp_Allvari1$Var.Importance$Dataset    = "impallvari8"
  rfN.imp_Allvari2$Var.Importance$Dataset    = "impallvari9"
  rfN.imp_Allvari3$Var.Importance$Dataset    = "impallvari18"
  rfN.imp_Allvari4$Var.Importance$Dataset    = "impallvari16"
  rfN.imp_Allvari5$Var.Importance$Dataset    = "impallvari4"
  
  NULL.all.VarImp.list <- list(rfN.all_Virtaxa$Var.Importance,
                               rfN.imp_Virtaxa$Var.Importance,
                               rfN.boruta_Virtaxa$Var.Importance,
                               rfN.maaslin_Virtaxa$Var.Importance,
                               rfN.imp_Allvari1$Var.Importance,
                               rfN.imp_Allvari2$Var.Importance,
                               rfN.imp_Allvari3$Var.Importance,
                               rfN.imp_Allvari4$Var.Importance,
                               rfN.imp_Allvari5$Var.Importance)

  NULL.all.VarImp <- NULL.all.VarImp.list %>% Reduce(function(d1, d2) full_join(d1, d2), .)
  NULL.all.VarImp$Dataset = NULL
  
  stopCluster(cl)
  registerDoSEQ()
  
  ## Plot predictive performance metrics------------------------------------------ 
  
  # We'd like to plot the distributions of accuracy, sensitivity, and specificity
  # of predictions on the validation cohort for each of the trained models. Also
  # plot null distributions of these parameters from the class label shuffling. 
  
  # For ease of faceting in plots, add columns DataSet indicating the base model
  # (e.g. all, Beta ..), and Taxa indicating if feature-selected taxa were 
  # included as features in the model, to prediction results (all.pred.Val). 
  
  all.pred.Val2 <- all.pred.Val %>% 
    separate(Data, into=c('DataSet', 'Taxa'), sep = ' plus ', 
             remove = FALSE, fill = 'right') %>% 
    mutate(Taxa = ifelse(is.na(Taxa), 'all VirTaxa', 
                         'Plus others'))
  write.table(c(DBO,all.pred.Val2), 
              "250416_ML_ImpVIr_early_pred_v1.2.csv", 
              append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
  
  ## Plot empirical and null variable importances---------------------------------
  
  # Add DataSet and Taxa variables to the VARIMP dataframes. These are inherited
  # from all.pred.Val2 since the data were generated simultaneously in 
  # train_rf_models(). 
  all.VarImp2 <- all.VarImp
  all.VarImp2$DataSet <- all.pred.Val2$Data
  all.VarImp2$Taxa <- all.pred.Val2$Taxa
  # all.VarImp2$X <- NULL
  
  all.VarImp_NULL2 <- NULL.all.VarImp
  all.VarImp_NULL2$DataSet <- all.pred.Val2$Data
  all.VarImp_NULL2$Taxa <- all.pred.Val2$Taxa
  # all.VarImp_NULL2$X<- NULL
  
  
  # Summarise Variable Importances by Model (DataSet), +/- GM features (Taxa)
  all.VarImp.mean <- all.VarImp2 %>% group_by(DataSet, Taxa) %>%
    summarise_all(mean) 
  
  all.VarImp.sd <- all.VarImp2 %>% group_by(DataSet, Taxa) %>%
    summarise_all(sd) 
  
  all.VarImpNULL.mean <- all.VarImp_NULL2 %>% group_by(DataSet, Taxa) %>%
    summarise_all(mean)
  
  all.VarImpNULL.sd <- all.VarImp_NULL2 %>% group_by(DataSet, Taxa) %>%
    summarise_all(sd)
  
  
  
  # Convert to long and merge into a single data frame
  VarImp.mean.long <- gather(all.VarImp.mean, Feature, Feature.MeanImp, 
                             colnames(all.VarImp.mean)[3]:colnames(all.VarImp.mean)[ncol(all.VarImp.mean)],
                             factor_key=TRUE)
  
  VarImp.sd.long <- gather(all.VarImp.sd, Feature, Feature.SDImp, 
                           colnames(all.VarImp.sd)[3]:colnames(all.VarImp.sd)[ncol(all.VarImp.sd)],
                           factor_key=TRUE)
  
  VarImp.NULL.mean.long <- gather(all.VarImpNULL.mean, Feature, Feature.NULLMeanImp, 
                                  colnames(all.VarImpNULL.mean)[3]:colnames(all.VarImpNULL.mean)[ncol(all.VarImpNULL.mean)],
                                  factor_key=TRUE)
  
  VarImp.NULL.sd.long <- gather(all.VarImpNULL.sd, Feature, Feature.NULLSDImp, 
                                colnames(all.VarImpNULL.sd)[3]:colnames(all.VarImpNULL.sd)[ncol(all.VarImpNULL.sd)],
                                factor_key=TRUE)
  
  
  VarImp.summ <- merge(VarImp.mean.long, VarImp.sd.long)
  VarImp.summ <- merge(VarImp.summ, VarImp.NULL.mean.long)
  VarImp.summ <- merge(VarImp.summ, VarImp.NULL.sd.long)
  
  
  # Fix feature names for plotting
  
  VarImp.summ <- VarImp.summ %>%
    mutate(Feature = str_remove(Feature, "_recent")) %>%
    mutate(Feature = str_remove(Feature, "Recent_")) %>%
    mutate(Feature = if_else(
      substring(Feature, nchar(Feature), nchar(Feature)) == "1", 
      substring(Feature, 1, nchar(Feature) - 1), 
      Feature
    ))
  varimp.labels <- unique(VarImp.summ$Feature)
  
  
  # FUNCTION plotVarImp takes long-format summary of variable importances (mean
  # and SD) for a given model, both the empiricial variable importances and
  # those from training after class label shuffling (to find null distributions of
  # variable importances). It plots variable importances (bar plots with SD error 
  # bars) ranked from most to least important, and overlays null mean and SD variable
  # importances from the permutational analysis (class label shuffling).
  # Required packages: ggplot2
  # Arguments:
  #  VarImp.df (dataframe)       = Mean and SD of variable importances for a model,
  #                                and optionally multiple models. Should include
  #                                column names: 'Feature.MeanImp', 'Feature.SDImp',
  #                                'Feature.NULLMeanImp' and 'Feature.NULLSDImp',  
  #                                as well as 'DataSet' (model names).
  #  Data.Set (str)              = Model (e.g. 'Alpha') to which to subset, if 
  #                                data collated for more than one model are provided.
  #  var.labels (named chr list) = Feature labels for plotting.
  # Return:
  #  Plot of ranked variable importances for a given model, both empirical and 
  #  null from class label shuffling. 
  
  plotVarImp <- function(VarImp.df, Data.Set, var.labels){
    
    p_varimp <- ggplot(subset(VarImp.df, DataSet== Data.Set & !is.na(Feature.MeanImp)), 
                       aes(x=reorder(Feature, Feature.MeanImp), y=Feature.MeanImp))+
      geom_col(fill='#66b2b2')+
      geom_point(color='black')+
      geom_errorbar(aes(ymin = Feature.MeanImp - Feature.SDImp, 
                        ymax = Feature.MeanImp + Feature.SDImp), width=0.7, size=0.8)+
      geom_point(aes(x=reorder(Feature, Feature.MeanImp), y=Feature.NULLMeanImp),
                 color='cyan', alpha=0.80)+
      geom_errorbar(aes(ymin = Feature.NULLMeanImp - Feature.NULLSDImp,
                        ymax = Feature.NULLMeanImp + Feature.NULLSDImp), width=0.5, size=0.6,
                    color='cyan', alpha=0.80)+
      ylab("Feature importance")+
      theme_classic()+
      theme(axis.text.y = element_text(size=10),
            axis.text.x = element_text(size=10),
            axis.title.x = element_text(size=10),
            axis.title.y = element_blank(),
            plot.title = element_blank(),
            panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
      facet_wrap(~DataSet, nrow=2)+
      # scale_x_discrete(labels = var.labels)+
      coord_flip()
  }
  
  
  # Main text models
  # pvar_bactaxa <- plotVarImp(VarImp.summ, 'BacTaxa', varimp.labels)
  # pvar_virtaxa <- plotVarImp(VarImp.summ, 'plus VirTaxa', varimp.labels)
  # pvar_impvirtaxa <- plotVarImp(VarImp.summ, 'plus impVirTaxa', varimp.labels)
  # pvar_feature <- plotVarImp(VarImp.summ, 'plus impAllVari', varimp.labels)
  # pvarimp_main <- ggarrange(pvar_all, pvar_bactaxa, 
  #                           pvar_virtaxa, pvar_impvirtaxa, nrow=2)
  
  write.table(c(DBO,VarImp.summ), 
              "250416_ML_ImpVIr_early_VarImpSumm_v1.2.csv", 
              append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
  
  ## Compare empirical and null variable importance------------------------------
  
  # Combine empirical and null variable importance data frames, subsetting to 
  # models that include taxonomic features
  all.VarImp3 <- subset(all.VarImp2, Taxa == 'all VirTaxa')
  all.VarImp3$Group <- 'Empirical'
  
  all.VarImp_NULL3 <- subset(all.VarImp_NULL2, Taxa == 'all VirTaxa')
  all.VarImp_NULL3$Group <- 'Null'
  
  all.VarImp.combined <- rbind(all.VarImp3, all.VarImp_NULL3)
  all.VarImp.combined$Group <- as.factor(all.VarImp.combined$Group)
  
  #all.VarImp.combined$RaceOther <- NULL
  # write.table(c(DBO,all.VarImp.combined), 
  #             paste("250416_ML_ImpVIr_early_VarImpCombined_DBO", DBO, "_v1.2.csv", sep = ""),
  #             append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

  # For each model (e.g. Alpha, Beta..), for each feature, perform t-test for
  # differences in importance by whether it was null (class label-shuffled) or 
  # empirical version of the model. 
  
  # Set up empty data frame to collate t-test results.
  varimp.tests <- data.frame(matrix(ncol=5, nrow=1000))
  colnames(varimp.tests) <- c('DataSet', 'Feature', 'P.val', 'Mean.Emp', 'Mean.Null')
  
  # Iterate through each model, and features included therein 
  models <- unique(all.VarImp.combined$DataSet)
  count = 0
  for (i in 1:length(models)){
    varimp.data <- subset(all.VarImp.combined, DataSet == models[i])
    
    varimp.data.long <- gather(varimp.data, Feature, Feature.Importance, 
                               colnames(varimp.data)[1]:colnames(varimp.data)[ncol(varimp.data)-3],
                               factor_key=TRUE)
    varimp.data.long <- na.omit(varimp.data.long)
    
    features <- as.character(unique(varimp.data.long$Feature))
    
    for (j in 1:length(features)) {
      count = count + 1
      print(count)
      students.t <- t.test(Feature.Importance ~ Group, 
                           data = subset(varimp.data.long, Feature == features[j]))
      
      newrow <- c(DataSet = models[i], Feature = features[j], 
                  P.val = students.t$p.value, Mean.Emp = students.t$estimate[[1]],
                  Mean.Null = students.t$estimate[[2]])
      
      varimp.tests[count, ] <- newrow
    }
  }
  
  
  # Subset to comparisons in which mean Empirical is greater than mean Null
  numvars <- c('P.val', 'Mean.Emp', 'Mean.Null')
  
  varimp.tests[, numvars] <- lapply(varimp.tests[,numvars], as.numeric)
  
  varimp.tests <- varimp.tests %>%
    mutate(Diff = Mean.Emp - Mean.Null)
  
  varimp.tests2 <- subset(varimp.tests, Diff > 0)
  
  # Apply BH (FDR) correction for multiple hypothesis testing and subset to 
  # significant comparisons. Add star nomenclature.
  varimp.tests2$P.val.adj <- p.adjust(varimp.tests2$P.val, method='BH')
  varimp.tests2 <- subset(varimp.tests2, P.val.adj < 0.05)
  
  varimp.tests2 <- varimp.tests2 %>% 
    mutate(Adj_Stars = if_else(P.val.adj <= 0.05 & P.val.adj > 0.01, '*',
                               if_else(P.val.adj <= 0.01 & P.val.adj > 0.001, '**', 
                                       if_else(P.val.adj <= 0.001, '***', 'NS'))))
  write.table(c(DBO,varimp.tests), 
              paste("250416_ML_ImpVIr_early_VarImpTest_DBO", DBO, "_v1.2.csv", sep = ""),
              append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")
}
