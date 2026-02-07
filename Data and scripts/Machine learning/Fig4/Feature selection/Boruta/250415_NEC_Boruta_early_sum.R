
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_Boruta')

write.table(cbind('Days_before_onset','Feature','Freq','Dataset'), 
            "250415_NEC_Boruta_early_sum.csv", 
            col.names = FALSE, row.names = FALSE, sep = ",")

for (DBO in 0:8){
# Load the information of the important variables
ImpVar.df        = read.csv(paste("250415_NEC_boruta_early_DBO", DBO, ".csv", sep = ""), header = T)

write.table(c(ImpVar.df), 
            paste("250415_NEC_Boruta_early_sum.csv", sep = ""),
            append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")

}
