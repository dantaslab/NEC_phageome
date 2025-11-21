library(ggplot2)
library(readxl)
library("dplyr")
set.seed(129)

data = read_excel("Prevalence_phage_ARG_class.xlsx")
colnames(data) = c("Phage_ARGclass", "Date_after_birth", "Prevalence")
data = filter(data, Phage_ARGclass!="FLUOROQUINOLONE" & Phage_ARGclass!="TRIMETHOPRIM")

plot_prevalence_AMRclass <- ggplot(data, aes(x=Date_after_birth, y=Prevalence, color=Phage_ARGclass)) + #geom_point() +
  geom_smooth(method = loess, alpha=0.2) +
  #geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Prevalence [% of cohort]") +
  xlab("Postnatal age [days]") + 
  #scale_color_discrete(name="Phage_AMRclass") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position="top", legend.title = element_text(size=0), legend.text = element_text(size = 9)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlim(0, 60) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  guides(color = guide_legend(nrow = 6))

plot_prevalence_AMRclass
ggsave("Prevalence_Phage_ARGclass.pdf", plot_prevalence_AMRclass, width=6, height=6)
