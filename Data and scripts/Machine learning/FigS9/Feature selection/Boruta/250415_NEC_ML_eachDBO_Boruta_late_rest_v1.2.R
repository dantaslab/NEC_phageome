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
dataset_all_Bacrestirep  = merge(dataset_all_resistome,dataset_all_iRep,by="row.names")
rownames(dataset_all_Bacrestirep) = dataset_all_Bacrestirep$Row.names
dataset_all_Bacrestirep$Row.names = NULL

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
dataset_all_BacrestirepMAX = merge(dataset_all_Bacrestirep,dataset_all_MAX,by="row.names")


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
dataset_all  = merge(dataset_all2,dataset_all_BacrestirepMAX,by="Row.names")
row.names(dataset_all) = dataset_all$Row.names
dataset_all = filter(dataset_all, !grepl("204-01|2156-01|415-01|2132-01|2234-01|234-01|
                                   421-01|2218-01|2238-02|49-01|257-01|2022-01|
                                   2077-01|2190-01|2228-01", Patient))
dataset_all = filter(dataset_all, Patient!="421-01" & Patient!="2077-01")
dataset_all = filter(dataset_all, NEC_onset2=="Late")

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

# write.table(cbind('Days_before_onset','Control.Control', 'Control.Case', 'Case.Control', 'Case.Case', 
#                   'Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyNull', 'AccuracyPValue', 
#                   'McnemarPValue', 'Sensitivity', 'Specificity', 'Pos.Pred.Value', 'Neg.Pred.Value', 
#                   'Precision', 'Recall', 'F1', 'Prevalence', 'Detection.Rate', 'Detection.Prevalence', 
#                   'Balanced.Accuracy', 'Data', 'DataSet', 'Taxa', 'Seed'), 
#             "241127_ML_allTaxa_MAX_pred_v1.2.csv", 
#             col.names = FALSE, row.names = FALSE, sep = ",")
# 
# write.table(cbind('Days_before_onset','DataSet', 'Taxa', 'Feature', 'Feature.MeanImp', 'Feature.SDImp', 'Feature.NULLMeanImp', 'Feature.NULLSDImp'), 
#             "241127_ML_allTaxa_MAX_VarImpSumm_v1.2.csv", 
#             col.names = FALSE, row.names = FALSE, sep = ",")

for (DBO in 0:8){
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
  fs.bactaxa.all = colnames(dataset_all_bac)[1:(length(colnames(dataset_all_bac))-3)]
  fs.virtaxa.all = colnames(dataset_all_phage)
  fs.bacrest.all = colnames(dataset_all_resistome)
  fs.bacirep.all = colnames(dataset_all_iRep)
  fs.bacpath.all = colnames(dataset_all_pathway)
  fs.bacMAX.all  = colnames(dataset_all_MAX)
  
  bactaxa.all = GM.caret[ , colnames(GM.caret) %in% fs.bactaxa.all]
  virtaxa.all = GM.caret[ , colnames(GM.caret) %in% fs.virtaxa.all] 
  bacrest.all = GM.caret[ , colnames(GM.caret) %in% c(fs.bacrest.all, "NEC_status")] 
  base.df = bacrest.all
  bacirep.all = GM.caret[ , colnames(GM.caret) %in% fs.bacirep.all]
  bacpath.all = GM.caret[ , colnames(GM.caret) %in% fs.bacpath.all]
  bacMAX.all  = GM.caret[ , colnames(GM.caret) %in% fs.bacMAX.all]
  feature.all = GM.caret[ , !colnames(GM.caret) %in% c(fs.bacrest.all, "NEC_status")] 
  
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
  
  # #  Exclude amyloid only
  # beta.nullvars <- c('PET_amyloid', 'CSF.ratio.ab42.ab40')
  
  
  # Create Data Subsets
  #  Main Text models:
  all.bactaxa <- createSubsets(base.df, all.nullvars, Taxa=bactaxa.all)
  all.virtaxa <- createSubsets(base.df, all.nullvars, Taxa=virtaxa.all)
  # all.bacrest <- createSubsets(base.df, all.nullvars, Taxa=bacrest.all)
  all.bacirep <- createSubsets(base.df, all.nullvars, Taxa=bacirep.all)
  all.bacpath <- createSubsets(base.df, all.nullvars, Taxa=bacpath.all)
  all.bacMAX  <- createSubsets(base.df, all.nullvars, Taxa=bacMAX.all)
  all.feature <- createSubsets(base.df, all.nullvars, Taxa=feature.all)
  
  
  # ## Partition out training cohort for feature selection--------------------------
  # set.seed(42)
  # #train_idx <- createDataPartition(meta.caret.im$amyloid.positive.AF, p = 0.6, list=FALSE)
  # train_idx <- createDataPartition(GM.caret$NEC_status, p = 0.75, list=FALSE)
  # 
  # # Access metadata (biomarkers) and taxanomic abundances for training cohort.
  # #meta.train <- meta.caret.im[train_idx, ]
  # GM.NEC.train <- GM.NEC[train_idx, ]
  
  # ## Iterative taxonomic feature selection on training cohort---------------------
  # 
  # # Define categorical variables that should not be normalized / scaled.
  # varsnot2norm <- factor.vars
  # 
  # #Pre-process data (center and scale)
  # preprocessrule <- preProcess(GM.NEC.train[, !(colnames(GM.NEC.train) %in% varsnot2norm)], 
  #                              method = c('center', 'scale'))
  # GM.NEC.train.p <- predict(preprocessrule, GM.NEC.train)
  
  
  # FUNCTION: callBoruta helper function (called by runFeatureSelection() below)
  # Required packages: Boruta, stringr
  # Arguments:
  #  taxa.data (dataframe)  = GM.NEC.train (tax abundances for training cohort)
  #  seed.Boruta (int)      = will be passed iteratively in defined range
  #
  # Return:
  #  list of names of feature-selected taxa for current iteration
  
  callBoruta <- function(taxa.data, seed.Boruta){
    set.seed(seed.Boruta)
    
    taxa.boruta <- Boruta(NEC_status~., data = taxa.data, maxRuns=500, doTrace=0)
    taxa.boruta.fix <- TentativeRoughFix(taxa.boruta)
    
    taxa.boruta.df <- data.frame('boruta' = taxa.boruta.fix$finalDecision)
    
    taxanames <- str_replace_all(rownames(subset(taxa.boruta.df, 
                                                 boruta=='Confirmed')), "`", "")
    
  }
  
  
  # FUNCTION: runFeatureSelection - iterative feature selection function
  # Required packages: Boruta, stringr
  # Arguments:
  #  taxa.train.data (dataframe)   = GM.NEC.train (passed to callBoruta())
  #  seed.range (num list)         = range for iteration, e.g. 1:100.
  #
  # Return:
  #  dataframe summarizing frequency at which unique taxa were feature-selected
  #  across all iterations (random seeds).
  
  runFeatureSelection <- function(taxa.train.data, seed.range){
    
    ## Partition out training cohort for feature selection--------------------------
    set.seed(123)
    #train_idx <- createDataPartition(meta.caret.im$amyloid.positive.AF, p = 0.6, list=FALSE)
    train_idx <- createDataPartition(taxa.train.data$NEC_status, p = 0.75, list=FALSE)
    
    # Access metadata (biomarkers) and taxanomic abundances for training cohort.
    #meta.train <- meta.caret.im[train_idx, ]
    GM.NEC.train <- taxa.train.data[train_idx, ]
    
    ## Iterative taxonomic feature selection on training cohort---------------------
    # Define categorical variables that should not be normalized / scaled.
    varsnot2norm <- factor.vars
    
    #Pre-process data (center and scale)
    preprocessrule <- preProcess(GM.NEC.train[, !(colnames(GM.NEC.train) %in% varsnot2norm)], 
                                 method = c('center', 'scale'))
    GM.NEC.train.p <- predict(preprocessrule, GM.NEC.train)
    
    fs.taxa.train <- vector('list', length(seed.range))
    for (i in seed.range){
      fs.taxa.train[[i]] <- callBoruta(GM.NEC.train.p, i)
    }
    
    fs.taxa.train.summ <- data.frame(table(unlist(fs.taxa.train)))
  } 
  
  
  require(doParallel)
  cl = makeCluster(24)
  registerDoParallel(cl)
  
  
  # Carry out iterative feature selection (may take a little while..~20 min)
  # fs.taxa <- runFeatureSelection(GM.NEC.train.p, 1:100)
  # fs.all.virtaxa <- runFeatureSelection(all.virtaxa$WithTaxa, 1:100)
  # fs.all.bactaxa <- runFeatureSelection(all.bactaxa$WithTaxa, 1:100)
  fs.all.bacrest <- runFeatureSelection(all.bactaxa$Base, 1:100)
  # fs.all.bacirep <- runFeatureSelection(all.bacirep$WithTaxa, 1:100)
  # fs.all.bacpath <- runFeatureSelection(all.bacpath$WithTaxa, 1:100)
  # fs.all.bacMAX  <- runFeatureSelection(all.bacMAX$WithTaxa, 1:100)
  fs.all.feature <- runFeatureSelection(all.feature$WithTaxa, 1:100)
  
  stopCluster(cl)
  registerDoSEQ()
  
  # Recommended to save workspace at this point.
  
  # Filter for taxa selected in > 25% of iterations.
  # fs.taxa.top <- subset(fs.taxa, Freq > 25, select='Var1')
  # fs.taxa.top
  # 
  # fs.all.virtaxa.top <- subset(fs.all.virtaxa, Freq > 25) #, select='Var1')
  # fs.all.virtaxa.top$Test_group = "VirTaxa"
  # 
  # fs.all.bactaxa.top <- subset(fs.all.bactaxa, Freq > 25) #, select='Var1')
  # fs.all.bactaxa.top$Test_group = "BacTaxa"
  
  fs.all.bacrest.top <- subset(fs.all.bacrest, Freq > 25) #, select='Var1')
  fs.all.bacrest.top$Test_group = "BacRest"
  
  # fs.all.bacirep.top <- subset(fs.all.bacirep, Freq > 25) #, select='Var1')
  # fs.all.bacirep.top$Test_group = "BaciRep"
  # 
  # fs.all.bacpath.top <- subset(fs.all.bacpath, Freq > 25) #, select='Var1')
  # fs.all.bacpath.top$Test_group = "BacPath"
  # 
  # fs.all.bacMAX.top  <- subset(fs.all.bacMAX, Freq > 25) #, select='Var1')
  # fs.all.bacMAX.top$Test_group = "BacMAX"
  
  fs.all.feature.top <- subset(fs.all.feature, Freq > 25) #, select='Var1')
  fs.all.feature.top$Test_group = "AllVars"
  
  fs.all.top = rbind(# fs.all.virtaxa.top,
                     # fs.all.bactaxa.top,
                     fs.all.bacrest.top,
                     # fs.all.bacirep.top,
                     # fs.all.bacpath.top,
                     # fs.all.bacMAX.top,
                     fs.all.feature.top)
  
  
  ## Plot relative abundances of feature-selected taxa----------------------------
  
  # Subset melted phyloseq dataframe to the feature-selected taxa.
  # GM.NEC.virtaxa.all <- subset(GM.caret, select = c(as.matrix(fs.all.top$Var1), "NEC_status"))
  # GM.NEC.virtaxa.all.df = gather(GM.NEC.virtaxa.all, "Feature", "Abundance", -NEC_status)

  write.table(c(DBO,fs.all.top), 
              paste("250415_NEC_boruta_late_rest_DBO", DBO, ".csv", sep = ""),
              append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")
  
}
