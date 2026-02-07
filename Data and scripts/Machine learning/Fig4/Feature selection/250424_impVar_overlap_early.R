library(dplyr) 
library(ggvenn)
library(RColorBrewer)

setwd("~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_Boruta")
vir_maaslin <- read.csv('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 3/Fig3A_extension/Fig3A_eachDBO/250416_maaslin2_phage_host_sig_features_early.csv', header = T) 

# vir_maaslin = filter(vir_maaslin, Days_before_NEC_onset!="0_1")
vir_maaslin.count = as.matrix(table(vir_maaslin$feature))
maaslin.var = paste(row.names(vir_maaslin.count)[vir_maaslin.count>=1],"vir",sep="_")
# maaslin.var = paste(unique(vir_maaslin$feature) ,"vir",sep="_")

varimp.tests2.all <- read.csv('Input/240424_ML_allTaxa_MAX_VarImpTest2_early_v1.2.csv', header = T) 
varimp.tests2.allvir <- filter(varimp.tests2.all, DataSet=="VirTaxa" & Days_before_onset>=0) # VirTaxa or All Features
varimp.tests2.allvir.count = as.matrix(table(varimp.tests2.allvir$Feature))
allvir.var = row.names(varimp.tests2.allvir.count)#[varimp.tests2.allvir.count>1]
# allvir.var = unique(varimp.tests2.allvir$Feature)

vir_boruta <- read.csv('250415_NEC_Boruta_early_sum.csv', header = T)
vir_boruta.feature = filter(vir_boruta, Dataset=="VirTaxa" & Days_before_onset>=0 & Freq > 25) # VirTaxa or AllVars
vir_boruta.feature.count = as.matrix(table(vir_boruta.feature$Feature))
boruta.var = row.names(vir_boruta.feature.count)#[vir_boruta.feature.count>1]
# boruta.var = unique(vir_boruta.feature$Feature)

intersect(allvir.var, boruta.var)
paste(intersect(allvir.var, boruta.var), collapse = "', '")

x.vir <- list(MaAsLin=maaslin.var , ML_NEC=allvir.var , Boruta=boruta.var)
p.vir <- ggvenn(x.vir, show_elements = F, label_sep = "\n", fill_color = brewer.pal(name="Set2",n=3),
       fill_alpha = 0.7,  show_outside = "none",show_percentage = F,text_size = 6)
p.vir
intersect(intersect(maaslin.var, allvir.var), boruta.var)
paste(intersect(intersect(maaslin.var, allvir.var), boruta.var), collapse = "', '")
# Reduce(intersect, x.vir)






