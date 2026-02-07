library(ggplot2)
library(readxl)
library(dplyr)
set.seed(129)


data = read_excel("Data_S2_v1.xlsx")

## Cumulative abx useage
plot_abx_cumu_DOL <- ggplot(data, aes(x=DOL, y=Abx_cumulative, color=NEC_status, linetype=NEC_onset)) + 
  # geom_point(alpha=0.3, size = 0.5) +
  geom_smooth(method = loess, alpha=0.3) +
  #geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Cumulative abx") +
  xlab("Postnatal age [days]") + 
  #scale_color_discrete(name="Host hostGenus") +
  # scale_color_brewer(palette = "Spectral") + 
  scale_color_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  theme(legend.position="none", legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 60) # + scale_y_continuous(breaks=c(0,25,50,75,100)) +
  # facet_wrap(~Host_genus, scale="fixed")
#scale_y_continuous(breaks=c(0,25,50,75,100))
plot_abx_cumu_DOL
ggsave("251204_CumulativeAbx_DOL_smooth.pdf", plot_abx_cumu_DOL, width=5, height=6.2)

plot_abx_cumu_DBO <- ggplot(data, aes(x=Days_before_NEC_onset, y=Abx_cumulative, color=NEC_status, linetype=NEC_onset)) + 
  # geom_point(alpha=0.5, size = 0.5) +
  geom_smooth(method = loess, alpha=0.3) +
  # geom_smooth(alpha=0.2, span = 0.6) +
  ylab("Cumulative abx") +
  xlab("Days before NEC onset") + 
  #scale_color_discrete(name="Host hostGenus") +
  # scale_color_brewer(palette = "Spectral") + 
  scale_color_manual(values=c("case"="#f28a5f","control"="#3291a4"))+
  theme(legend.position=c(0.25,0.85), legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + xlim(0, 60) +
  scale_x_reverse(limits = c(60, 0))
  # facet_wrap(~NEC_onset, scale="fixed")
plot_abx_cumu_DBO
ggsave("251204_CumulativeAbx_DBO_smooth.pdf", plot_abx_cumu_DBO, width=5, height=6.2)

## Overtime cumulative abx -1 day usage
# 1. Group the data by patient and 2. reorder by DOL, then 3. mius the value previously? not day by day frequency?
# Or calculte the slope of the function?
data %>%
  mutate(value_minus_previous = Abx_cumulative - lag(Abx_cumulative))

A=data$Abx_cumulative - lag(data$Abx_cumulative)

## Average abx usage over time


