library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library("compositions")
library(reshape)
library(reshape2)
library("xlsx")
library(data.table)
library(dplyr)
library(vegan)

setwd('~/Google Drive/WashU/Phageome/NEC/Data analysis/ShortBRED')
##-----------------------------------------------------------------------------------------------------------------------------------
data_all_nozero <- read_excel("joint_output_nozero.xlsx")
data_nozero<- subset(data_all_nozero, select = c(Sample, Family, Count))


cast.data_nozero = acast(data_nozero, Sample~Family, value.var="Count")
cast.data_nozero[is.na(cast.data_nozero)] <- 0
# write.xlsx2(cast.data_nozero , "240831_NEC_ShortBRED.xlsx", sheetName = "Sheet1", 
#             col.names = TRUE, row.names = FALSE, append = FALSE)

write.table(cast.data_nozero, file = "240831_NEC_ShortBRED.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

NEC_DF <- read_excel("240831_NEC_ShortBRED.xlsx")
NEC_DF = filter(NEC_DF, !grepl("204-01|2156-01|415-01|2132-01|2234-01|234-01|
                               421-01|2218-01|2238-02|49-01|257-01|2022-01|
                               2077-01|2190-01|2228-01", Patient))
NEC_DF = filter(NEC_DF, Patient!="421-01" & Patient!="2077-01")
write.table(NEC_DF, file = "240831_NEC_ShortBRED_V2.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
##-----------------------------------------------------------------------------------------------------------------------------------
# data_all_nozero <- read.delim("joint_output.txt", header = TRUE)
# 
# cast.data_all_nozero  = acast(data_all_nozero , Sample~Family, value.var="Count")
# write.xlsx2(cast.data_all_nozero , "240831_NEC_ShortBRED.xlsx", sheetName = "Sheet1", 
#             col.names = TRUE, row.names = FALSE, append = FALSE)
# write_csv2(as.data.frame(cast.data_all_nozero), "240831_NEC_ShortBRED.csv")

##-----------------------------------------------------------------------------------------------------------------------------------
# Scripts from Anna
#look at overall ARG burden between groups
##load shortbred output
arg.full <- read.delim('joint_output.txt')
arg.wide <- dcast(arg.full, Sample ~ Family, value.var = 'Count')
row.names(arg.wide) <- arg.wide$Sample
arg.wide <- arg.wide[, c(-1)]

##sum counts to get RPKM per sample
arg.wide$RPKM <- rowSums(arg.wide)

##get ARG richness
arg.wide$Richness <- specnumber(arg.wide[, 1:ncol(arg.wide)-1])

#look at associations by ARG class
##group by class
arg.map <- read.delim('ShortBRED_mapping.txt',
                      header = F)
arg.byclass <- merge(arg.full, arg.map, by.x = 'Family', by.y = 'V1') %>% 
  dplyr::rename(Class = V2)
arg.byclass.sum <- arg.byclass %>% 
  dplyr::group_by(Sample, Class) %>% 
  dplyr::summarise(Count = sum(Count))

##convert to wide format for MaAsLin2
arg.byclass.wide <- dcast(arg.byclass.sum, Sample ~ Class, value.var = 'Count')
arg.byclass.wide[is.na(arg.byclass.wide)] <- 0
row.names(arg.byclass.wide) <- arg.byclass.wide$Sample
arg.byclass.wide$Sample <- NULL

arg.byclass.wide = arg.byclass.wide[,colSums(arg.byclass.wide)>1]
arg.byclass.wide$VIRULENCE <- NULL


NEC_arg.df = merge(arg.byclass.wide, dataset_all_resistome, by="row.names")

write.table(NEC_arg.df, file = "241104_NEC_ShortBRED_byclass.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

