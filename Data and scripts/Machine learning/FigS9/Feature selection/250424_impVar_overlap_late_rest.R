library(dplyr) 
library(ggvenn)
library(RColorBrewer)

setwd("~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_Boruta")
rest_maaslin <- read.csv('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 3/Fig3A_extension/Fig3A_eachDBO/250416_maaslin2_resistome_sig_features_late.csv', header = T)
# # 241112_maaslin2_phage_host_sig_features.csv
# # 241202_maaslin2_GM_MTX_sig_features_LMLOG.csv
# # 241201_maaslin2_GM_MTX_sig_features_DBO8_LMLOG.csv
# rest_maaslin = filter(rest_maaslin, Days_before_NEC_onset!="0_1")
# rest_maaslin.count = as.matrix(table(rest_maaslin$feature))
# maaslin.var = paste(row.names(rest_maaslin.count)[rest_maaslin.count>1],"vir",sep="_")
maaslin.var = paste(unique(rest_maaslin$feature) ,"rest",sep="_")

varimp.tests2.all <- read.csv('Input/241129_ML_sallTaxa_MAX_VarImpTest2_late_rest_v1.2.csv', header = T) 
varimp.tests2.allrest <- filter(varimp.tests2.all, DataSet=="BacRest" & Days_before_onset>=0) # BacRest or All Features
# varimp.tests2.allrest.count = as.matrix(table(varimp.tests2.allrest$Feature))
# allrest.var = row.names(varimp.tests2.allrest.count)[varimp.tests2.allrest.count>1]
allrest.var = unique(varimp.tests2.allrest$Feature)

rest_boruta <- read.csv('250415_NEC_Boruta_late_rest_sum.csv', header = T)
rest_boruta.feature = filter(rest_boruta, Dataset=="BacRest" & Days_before_onset>=0 & Freq > 25) # BacRest or AllVars
# rest_boruta.feature.count = as.matrix(table(rest_boruta.feature$Feature))
# boruta.var = row.names(rest_boruta.feature.count)[rest_boruta.feature.count>1]
boruta.var = unique(rest_boruta.feature$Feature)

intersect(allrest.var, boruta.var)
paste(intersect(allrest.var, boruta.var), collapse = "', '")

x.rest <- list(MaAsLin=maaslin.var, ML_NEC=allrest.var, Boruta=boruta.var)
p.rest <- ggvenn(x.rest, show_elements = F, label_sep = "\n", fill_color = brewer.pal(name="Set2",n=3),
       fill_alpha = 0.7,  show_outside = "none",show_percentage = F,text_size = 6)
p.rest
intersect(intersect(maaslin.var, allrest.var), boruta.var)
paste(intersect(intersect(maaslin.var, allrest.var), boruta.var), collapse = "', '")
# Reduce(intersect, x.rest)

# dfrE: https://card.mcmaster.ca/ontology/39309




