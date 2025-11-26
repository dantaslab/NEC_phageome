library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library('caret')

require(doParallel)
cl = makeCluster(6)
registerDoParallel(cl)

df_phage1 = read.csv('BC_phyloseq_hostGenus_between_all_DBO.csv',
                    header = T, row.names = 1)
df_phage1$Source = "Phage Host Genus"

df_bacteria1 = read.csv('BC_phyloseq_bacGenus_between_all_DBO.csv',
                       header = T, row.names = 1)
df_bacteria1$Source = "Bacteria Genus"

df_phage2 = read.csv('BC_phyloseq_vOTU_between_all_DBO.csv',
                    header = T, row.names = 1)
df_phage2$Source = "Phage vOTU"

df_bacteria2 = read.csv('BC_phyloseq_bacSpecies_between_all_DBO.csv',
                       header = T, row.names = 1)
df_bacteria2$Source = "Bacteria Species"


df_between_NECstatus = rbind(df_phage1, df_bacteria1, df_phage2, df_bacteria2)
df_between_NECstatus$NEC_group[df_between_NECstatus$NEC_status_patient1=="case" & df_between_NECstatus$NEC_status_patient2=="case"] = "Case"
df_between_NECstatus$NEC_group[df_between_NECstatus$NEC_status_patient1=="control" & df_between_NECstatus$NEC_status_patient2=="control"] = "Control"
df_between_NECstatus$NEC_group[df_between_NECstatus$NEC_status_patient1!=df_between_NECstatus$NEC_status_patient2] = "Case-Control"

df_between_NECstatus = filter(df_between_NECstatus, DBO_patient1<15)

df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=14] = "13-14"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=12] = "11-12"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=10] = "9-10"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=8] = "7-8"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=6] = "5-6"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=4] = "3-4"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2<=2] = "1-2"
df_between_NECstatus$DBO_group[df_between_NECstatus$DBO_patient2==0]  = "0"

##------------------------------------------------------------------------------------------------------------------------
library(rstatix)
stat.test <- df_between_NECstatus %>%
  group_by(DBO_group, Source) %>%
  t_test(BrayCurtis_dissimilarity ~ NEC_group, ref.group = "Control")
stat.test$statistic_rev = -stat.test$statistic
stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"] = ""

ggheatmap1 <- ggplot(stat.test, aes(DBO_group, Source, fill = statistic_rev))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#006b93", high = "#e32633", mid = "white", 
                       midpoint = 0,
                       # limits = c(-8,17), values = rescale(c(-8,0,17)), 
                       space = "Lab", 
                       name="T-Test\nt-score") + #\n[compare to Case-Control]
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1)) +
  facet_wrap(~group2, scale="fixed") +
  coord_fixed() + scale_x_discrete(limits = c("13-14","11-12","9-10","7-8","5-6","3-4","1-2","0"))
print(ggheatmap1)

ggheatmap_sig1 <- ggheatmap1 + 
  geom_text(aes(DBO_group, Source), label = stat.test$p.adj.signif, color = "black", size = 4) +
  theme(
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(0.6, 0.7),
    legend.direction = "vertical") +
  xlab("Days before NEC onset") +
  ylab("") +
  facet_wrap(~group2, scale="fixed")

ggheatmap_sig1
ggsave("Compare_BC_to_Control_DBO0-14.pdf", ggheatmap_sig1, width=8, height=2.5)

stopCluster(cl)
registerDoSEQ()

