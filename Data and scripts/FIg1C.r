library(ggplot2)
library(readxl)
library("dplyr")
library("xlsx")
set.seed(129)


data = read.csv('Phage_relAbundance_host_genus.csv',
                header = T)

data[,c(2:154)] = ifelse(data[,c(2:154)] != 0, 1, 0)

Sample_count = table(data$DOL)

Clostridium_count = table(data$Clostridium, data$DOL)
Clostridium_prevalence = Clostridium_count[2,]/Sample_count*100

Cutibacterium_count = table(data$Cutibacterium, data$DOL)
Cutibacterium_prevalence = Cutibacterium_count[2,]/Sample_count*100

Enterobacter_count = table(data$Enterobacter, data$DOL)
Enterobacter_prevalence = Enterobacter_count[2,]/Sample_count*100

Enterococcus_count = table(data$Enterococcus, data$DOL)
Enterococcus_prevalence = Enterococcus_count[2,]/Sample_count*100

Escherichia_count = table(data$Escherichia, data$DOL)
Escherichia_prevalence = Escherichia_count[2,]/Sample_count*100

Klebsiella_count = table(data$Klebsiella, data$DOL)
Klebsiella_prevalence = Klebsiella_count[2,]/Sample_count*100

Staphylococcus_count = table(data$Staphylococcus, data$DOL)
Staphylococcus_prevalence = Staphylococcus_count[2,]/Sample_count*100

Streptococcus_count = table(data$Streptococcus, data$DOL)
Streptococcus_prevalence = Streptococcus_count[2,]/Sample_count*100

Veillonella_count = table(data$Veillonella, data$DOL)
Veillonella_prevalence = Veillonella_count[2,]/Sample_count*100

All_prevalence = data.frame(c(Clostridium_prevalence, Cutibacterium_prevalence, Enterobacter_prevalence, 
                              Enterococcus_prevalence, Escherichia_prevalence, Klebsiella_prevalence, 
                              Staphylococcus_prevalence, Streptococcus_prevalence, Veillonella_prevalence))
colnames(All_prevalence) = "Prevalence"
Genus_name = c("Clostridium", "Cutibacterium", "Enterobacter","Enterococcus", "Escherichia", 
               "Klebsiella", "Staphylococcus", "Streptococcus", "Veillonella")

All_prevalence$Host_genus = rep(Genus_name, times=c(74,74,74,74,74,74,74,74,74))
All_prevalence$Date_after_birth = as.numeric(rep(row.names(Sample_count), times=9))


plot_prevalence_genus <- ggplot(All_prevalence, aes(x=Date_after_birth, y=Prevalence, color=Host_genus)) + #geom_point() +
  geom_smooth(method = loess, alpha=0.2) +
  #geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Prevalence [% of cohort]") +
  xlab("Postnatal age [days]") + 
  #scale_color_discrete(name="Host genus") +
  scale_color_brewer(palette = "Spectral") + #Spectral
  theme(legend.position="top", legend.title = element_text(size=0), legend.text = element_text(size = 12, face="italic")) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 60) +
  scale_y_continuous(breaks=c(0,25,50,75,100))

plot_prevalence_genus 
ggsave("Prevalence_Phage_HostGenus.pdf", plot_prevalence_genus, width=6, height=5)




