

setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_Boruta/Input')

boruta_df = read.csv("Boruta/250215_NEC_Boruta_sum.csv",header = 1)
boruta_Virtaxa = filter(boruta_df, Dataset == "VirTaxa")
boruta_Virtaxa_count  = as.data.frame(table(boruta_Virtaxa$Feature, boruta_Virtaxa$Dataset))
boruta_Virtaxa_count_desc = arrange(boruta_Virtaxa_count, desc(Freq))

# boruta_all = filter(boruta_df, Dataset == "AllVars")
# boruta_all_count  = as.data.frame(table(boruta_all$Feature, boruta_all$Dataset))
# boruta_all_count_desc = arrange(boruta_all_count, desc(Freq))

colnames(boruta_Virtaxa_count) = c("boruta_","DataSet", "Count")



VarImp_df = read.csv("241127_ML_sallTaxa_MAXVarImpTest2_v1.2.csv",header = 1)
VarImp_Virtaxa = filter(VarImp_df, DataSet == "VirTaxa")
VarImp_Virtaxa_count  = as.data.frame(table(VarImp_Virtaxa$Feature, VarImp_Virtaxa$DataSet))
VarImp_Virtaxa_count_desc = arrange(VarImp_Virtaxa_count, desc(Freq))

VarImp_Bactaxa = filter(VarImp_df, DataSet == "Plus BacTaxa")
VarImp_Bactaxa_count  = as.data.frame(table(VarImp_Bactaxa$Feature, VarImp_Bactaxa$DataSet))
VarImp_Bactaxa_count_desc = arrange(VarImp_Bactaxa_count, desc(Freq))

# colnames(boruta_Virtaxa_count) = c("boruta_","DataSet", "Count")


write.table(VarImp_count_desc, "250223_ML_allmeta_VarImpTest2SUM_v1.2.csv", col.names = TRUE, row.names = FALSE, sep = ",")