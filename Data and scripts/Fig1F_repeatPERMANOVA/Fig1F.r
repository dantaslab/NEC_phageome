library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr) 
library(readxl)
library(ggplot2)
library(ggpattern)
library('caret')
library(reshape2)
set.seed(126)

data = read_excel("Data_S4_v1.xlsx", 
                  sheet = "Data_S4")

data$vOTU = data$`-LOG10_BH P_value_vOTU`
plot_bar_RA <-
  ggbarplot(data, x = "Variable", y = "vOTU",
            fill = "Group",               # change fill color by cyl
            color = "black",            # Set bar border colors to white
            legend = "right",
            palette = "Spectral",            # jco journal color palett. see ?ggpar; Spectral
            sort.val = "desc",           # Sort the value in dscending order
            sort.by.groups = TRUE,      # Sort inside each group
            x.text.angle = 90,           # Rotate vertically x axis texts
            alpha = 0.75
  ) +
    geom_hline(yintercept=-log10(0.1), linetype="dashed", color = "black", size = 0.75) +
    xlab("Features") +
    ylab("-log(P value)")

plot_bar_RA 
ggsave("Bar_PERMANOVA_Pvalue_metadata_vOTU.pdf", plot_bar_RA, width=12, height=4)

