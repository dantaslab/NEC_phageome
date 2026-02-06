library(ggplot2)
library(readxl)
library("dplyr")
library("xlsx")
set.seed(129)


data = read.csv('Phage_relAbundance_family.csv',
                header = T)

data[,c(2:15)] = ifelse(data[,c(2:15)] != 0, 1, 0)

Sample_count = table(data$DOL)

Autographiviridae_count = table(data$Autographiviridae, data$DOL)
Autographiviridae_prevalence = Autographiviridae_count[2,]/Sample_count*100

Guelinviridae_count = table(data$Guelinviridae, data$DOL)
Guelinviridae_prevalence = Guelinviridae_count[2,]/Sample_count*100

Inoviridae_count = table(data$Inoviridae, data$DOL)
Inoviridae_prevalence = Inoviridae_count[2,]/Sample_count*100

Microviridae_count = table(data$Microviridae, data$DOL)
Microviridae_prevalence = Microviridae_count[2,]/Sample_count*100

Myoviridae_count = table(data$Myoviridae, data$DOL)
Myoviridae_prevalence = Myoviridae_count[2,]/Sample_count*100

Podoviridae_count = table(data$Podoviridae, data$DOL)
Podoviridae_prevalence = Podoviridae_count[2,]/Sample_count*100

Rountreeviridae_count = table(data$Rountreeviridae, data$DOL)
Rountreeviridae_prevalence = Rountreeviridae_count[2,]/Sample_count*100

Siphoviridae_count = table(data$Siphoviridae, data$DOL)
Siphoviridae_prevalence = Siphoviridae_count[2,]/Sample_count*100

All_prevalence = data.frame(c(Autographiviridae_prevalence, Guelinviridae_prevalence, Inoviridae_prevalence, 
                              Microviridae_prevalence, Myoviridae_prevalence, 
                              Podoviridae_prevalence, Rountreeviridae_prevalence, Siphoviridae_prevalence))
colnames(All_prevalence) = "Prevalence"
Genus_name = c("Autographiviridae", "Guelinviridae", "Inoviridae", "Microviridae", 
               "Myoviridae", "Podoviridae", "Rountreeviridae", "Siphoviridae")

All_prevalence$Phage_family = rep(Genus_name, times=c(74,74,74,74,74,74,74,74))
All_prevalence$Date_after_birth = as.numeric(rep(row.names(Sample_count), times=8))


plot_prevalence_genus <- ggplot(All_prevalence, aes(x=Date_after_birth, y=Prevalence, color=Phage_family)) + #geom_point() +
  geom_smooth(method = loess, alpha=0.2) +
  #geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Prevalence [% of cohort]") +
  xlab("Postnatal age [days]") + 
  #scale_color_discrete(name="Host genus") +
  scale_color_brewer(palette = "Accent") +
  theme(legend.position="top", legend.title = element_text(size=0), legend.text = element_text(size = 12, face="italic")) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 60) +
  scale_y_continuous(breaks=c(0,25,50,75,100))

plot_prevalence_genus 
ggsave("Prevalence_Phage_Family.pdf", plot_prevalence_genus, width=6, height=5)




