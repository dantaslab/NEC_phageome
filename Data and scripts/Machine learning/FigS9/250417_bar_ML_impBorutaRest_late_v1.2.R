library(ggplot2)
library(dplyr) 
library(forcats)
library(tidyr)
library(ggpubr) 
library(rstatix)
library(stringr)

setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_ImpVar/250417_ML_ImpVIr_late_rest_Test3_re')

## Plot empirical and null variable importances---------------------------------

# Add DataSet and Taxa variables to the VARIMP dataframes. These are inherited
# from all.pred.Val2 since the data were generated simultaneously in 
# train_rf_models(). 

## Load data
VarImp.summ <- read.csv('250417_ML_ImpVIr_late_rest_VarImpSumm_v1.2.csv', header = T)
unique(VarImp.summ$DataSet)
# VarImp.summ <- filter(VarImp.summ, Days_before_onset==DBO)
varimp.labels <- unique(VarImp.summ$Feature)

# VarImp.summ.combine <- VarImp.summ %>%
#                        group_by(DataSet,Feature) %>%
#                        summarise(Feature.MeanImp=mean(Feature.MeanImp),
#                                  Feature.SDImp=mean(Feature.SDImp),
#                                  Feature.NULLMeanImp=mean(Feature.NULLMeanImp),
#                                  Feature.NULLSDImp=mean(Feature.NULLSDImp))


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
               color='pink', alpha=1,shape=18,size=4)+
    geom_errorbar(aes(ymin = Feature.NULLMeanImp - Feature.NULLSDImp,
                      ymax = Feature.NULLMeanImp + Feature.NULLSDImp), width=0.5, size=0.6,
                  color='pink', alpha=0.80)+
    ylab("Feature importance")+
    theme_classic()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
    facet_wrap(~Days_before_onset, nrow=2)+
    # facet_wrap(~DataSet, nrow=2)+
    # scale_x_discrete(labels = var.labels)+
    coord_flip()
}
# 
# 
# # Main text models
unique(VarImp.summ$DataSet) #"allMaaslinBIN" "allMaaslinLM" "allmeta" "dbo8MaaslinBIN" "dbo8MaaslinLM"  "impAllVari" "impVirTaxa" "VirTaxa"
DBO=0
VarImp.summ_DBO <- filter(VarImp.summ, Days_before_onset==DBO)
pvar_meta <- plotVarImp(VarImp.summ_DBO, 'BorutaBacRest', varimp.labels)
pvar_meta

# pvarimp_main <- ggarrange(pvar_allmeta, pvar_nomedication, pvar_noabx, 
#                           pvar_nodiet, pvar_nobase, pvar_noabxbase, nrow=2)


## Compare empirical and null variable importance------------------------------
# columns= c("Days_before_onset","DataSet","Feature","P.val","Mean.Emp","Mean.Null","Diff") 
# varimp.tests = data.frame(matrix(nrow = 0, ncol = length(columns))) 
# colnames(varimp.tests) = columns
# for (DBO in 0:8){
# varimp.tests.DBO <- read.csv(paste("241112_ML_selectTaxa_allbac_VarImpTest_DBO",DBO,"_v1.2.csv", sep=""), header = T)
# varimp.tests.DBO <- filter(varimp.tests.DBO, DataSet == "BacTaxa")
# colnames(varimp.tests.DBO)[1] = "Days_before_onset"
# varimp.tests = rbind(varimp.tests,varimp.tests.DBO)
# }

# write.table(cbind('Days_before_onset','DataSet', 'Feature', 'P.val', 'Mean.Emp', 'Mean.Null', 'Diff', 'P.val.adj', 'Adj_Stars'), 
#             "250415_ML_allmeta_late_VarImpTest2_v1.2.csv", 
#             col.names = FALSE, row.names = FALSE, sep = ",")

for (DBO in 0:8){
  dataset = "BorutaBacRest"
    varimp.tests <- read.csv(paste("250417_ML_ImpVIr_late_rest_VarImpTest_DBO",DBO,"_v1.2.csv", sep=""), header = T)
    varimp.tests <- filter(varimp.tests, DataSet == dataset)
    varimp.tests2 <- subset(varimp.tests, Diff > 0)
    
    # Apply BH (FDR) correction for multiple hypothesis testing and subset to
    # significant comparisons. Add star nomenclature.
    varimp.tests2$P.val.adj <- p.adjust(varimp.tests2$P.val, method='BH')
    varimp.tests2 <- subset(varimp.tests2, P.val.adj < 0.05)
    
    varimp.tests2 <- varimp.tests2 %>%
      mutate(Feature = str_remove(Feature, "_recent")) %>%
      mutate(Feature = str_remove(Feature, "Recent_")) %>%
      mutate(Feature = if_else(
        substring(Feature, nchar(Feature), nchar(Feature)) == "1",
        substring(Feature, 1, nchar(Feature) - 1),
        Feature
      )) %>%
      mutate(Adj_Stars = if_else(P.val.adj <= 0.05 & P.val.adj > 0.01, '*',
                                 if_else(P.val.adj <= 0.01 & P.val.adj > 0.001, '**',
                                         if_else(P.val.adj <= 0.001, '***', 'NS'))))
    
    # write.table(varimp.tests2, 
    #             "250415_ML_allmeta_late_VarImpTest2_v1.2.csv",
    #             append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
    
    # pvar_meta + geom_text(data = varimp.tests2, aes(x = Feature, y = 0.2, label = Adj_Stars), size = 7, vjust = 0.75, color = "white")
    
    VarImp.summ.sub = filter(VarImp.summ, Feature %in% unique(varimp.tests2$Feature) & Days_before_onset==DBO & DataSet == dataset)
    pvar.sub <- plotVarImp(VarImp.summ.sub, dataset, unique(VarImp.summ.sub$Feature))
    pvar.sub
    plot = pvar.sub + geom_text(data = varimp.tests2, aes(x = Feature, y = 0.2, label = Adj_Stars), size = 7, vjust = 0.75, color = "black")
    plot
    ggsave(paste("250417_bar_", dataset,"_ImpVar_late_DBO",DBO ,"_v1.2.pdf", sep=""), pvar.sub , width=20, height=8)
  }
# }

# meta_VarImp_df = read.csv("250415_ML_allmeta_late_VarImpTest2_v1.2.csv",header = 1)
# VarImp_count = as.data.frame(table(meta_VarImp_df$Feature, meta_VarImp_df$DataSet))
# colnames(VarImp_count) = c("Meta_VarImp","DataSet", "Count")
# VarImp_count_desc = arrange(VarImp_count, desc(Count))
# write.table(VarImp_count_desc, "250415_ML_allmeta_late_VarImpTest2SUM_v1.2.csv", col.names = TRUE, row.names = FALSE, sep = ",")
# 
