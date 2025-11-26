library("compositions")
library("vegan")
library(readxl)
library(ggplot2)
library("dplyr")
library(xlsx)
library(ggpubr) 
set.seed(123)

SM_DF <- read.csv('NEC_bacteria_genus_RA.csv',
                  header = T,
                  row.names = 1)

metadata<- subset(SM_DF, select = c(NEC_onset2, NEC_status, DOL, Days_before_NEC_onset, Patient))
SM_DF <- subset(SM_DF, select = -c(NEC_onset2, NEC_status, DOL, Days_before_NEC_onset, Patient))

##--------------------------------------------------------------------------------------------------------
## Calculate and compare richness
taxa_rich <- as.data.frame(specnumber(SM_DF))
colnames(taxa_rich)[1] <- "Richness_bacGenus"
taxa_rich_meta <- merge(taxa_rich, metadata, by = "row.names")

taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL>50]  = "50-"
taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL<=50] = "41-50"
taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL<=40] = "31-40"
taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL<=30] = "21-30"
taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL<=20] = "11-20"
taxa_rich_meta$DOL_Group[taxa_rich_meta$DOL<=10] = "1-10"

plot_richness_DOL <- ggplot(taxa_rich_meta, aes(x=DOL_Group, y=Richness_bacGenus, fill=NEC_status)) +
  geom_boxplot() +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25), size=0.5, alpha=1, color="gray") + # , color="gray", size=0.5, alpha=1) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position="none") +
  xlab("Postnatal age [days]") +
  ylab("Bacteria Genus Richness") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~NEC_onset2, scale="fixed", nrow=2) +
  stat_compare_means(
    aes(DOL_Group = NEC_status), label = "p.format", # p.format or p.signif
    method = "wilcox.test", label.y = 33, size = 3)  + ylim(-5,36)
plot_richness_DOL
ggsave("Box_Richness_bacGenus_DOL.pdf", plot_richness_DOL, width=7, height=6)


##--------------------------------------------------------------------------------------------------------
## Calculate and compare SD
taxa_SD <- as.data.frame(diversity(SM_DF, index="shannon"))
colnames(taxa_SD)[1] <- "ShannonDiversity_BacGenus"
taxa_SD_meta <- merge(taxa_SD, metadata, by = "row.names")

taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL>50]  = "50-"
taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL<=50] = "41-50"
taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL<=40] = "31-40"
taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL<=30] = "21-30"
taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL<=20] = "11-20"
taxa_SD_meta$DOL_Group[taxa_SD_meta$DOL<=10] = "1-10"

plot_SD_DOL <- ggplot(taxa_SD_meta, aes(x=DOL_Group, y=ShannonDiversity_BacGenus, fill=NEC_status)) +
  geom_boxplot() +
  scale_fill_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25), size=0.5, alpha=1, color="gray") + # , color="gray", size=0.5, alpha=1) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),legend.position=c(0.85,0.85)) +
  xlab("Postnatal age [days]") +
  ylab("Bacteria Genus H'") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~NEC_onset2, scale="fixed", nrow=2) +
  stat_compare_means(
    aes(DOL_Group = NEC_status), label = "p.format", # p.format or p.signif
    method = "wilcox.test", label.y = max(taxa_SD_meta$ShannonDiversity_BacGenus) * 1.05, size = 3)
plot_SD_DOL
ggsave("Box_SD_bacGenus_DOL.pdf", plot_SD_DOL, width=7, height=6)


