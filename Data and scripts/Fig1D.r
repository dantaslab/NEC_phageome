library(ggplot2)
library(readxl)
set.seed(129)

data = read_excel("Prevalence_metabolism.xlsx", col_names = F)
colnames(data) = c("Phage_metabolism", "Date_after_birth", "Prevalence")

# Folding, sorting and degradation => Protein quality control
plot_prevalence_phageAMG <- ggplot(data, aes(x=Date_after_birth, y=Prevalence, color=Phage_metabolism)) + #geom_point() +
  geom_smooth(method = loess, alpha=0.2) +
  #geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Prevalence [% of cohort]") +
  xlab("Postnatal age [days]") + 
  #scale_color_discrete(name="Metabolism") +
  scale_color_brewer(palette = "Paired") +
  theme(legend.position="top", legend.title = element_text(size=0), legend.text = element_text(size = 9)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 60) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  guides(color = guide_legend(nrow = 6))

plot_prevalence_phageAMG
ggsave("Prevalence_Phage_Metabolism.pdf", plot_prevalence_phageAMG, width=6, height=6.4)
