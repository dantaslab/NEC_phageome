library(ggplot2)
library(dplyr) 
library(forcats)
library(tidyr)
library(ggpubr) 
library(rstatix)
library(plyr)

## Plot predictive performance metrics------------------------------------------ 

# We'd like to plot the distributions of accuracy, sensitivity, and specificity
# of predictions on the validation cohort for each of the trained models. Also
# plot null distributions of these parameters from the class label shuffling. 

# For ease of faceting in plots, add columns DataSet indicating the base model
# (e.g. allmeta, Beta ..), and Taxa indicating if feature-selected taxa were 
# included as features in the model, to prediction results (all.pred.Val). 

## Load data
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/BW/SingleDB')

# smaller
all.pred.Val.smaller.750 <- read.csv('BW750/250813_ML_allTaxa_MTX_singledb_smaller750_pred_v1.2.csv', header = T)
all.pred.Val.smaller.750$BW = "750"
all.pred.Val.smaller.800 <- read.csv('BW800/250813_ML_allTaxa_MTX_singledb_smaller800_pred_v1.2.csv', header = T)
all.pred.Val.smaller.800$BW = "800"
all.pred.Val.smaller.850 <- read.csv('BW850/250813_ML_allTaxa_MTX_singledb_smaller850_pred_v1.2.csv', header = T)
all.pred.Val.smaller.850$BW = "850"
all.pred.Val.smaller.900 <- read.csv('BW900/250813_ML_allTaxa_MTX_singledb_smaller900_pred_v1.2.csv', header = T)
all.pred.Val.smaller.900$BW = "900"
all.pred.Val.smaller.950 <- read.csv('BW950/250813_ML_allTaxa_MTX_singledb_smaller950_pred_v1.2.csv', header = T)
all.pred.Val.smaller.950$BW = "950"
all.pred.Val.smaller.1000 <- read.csv('BW1000/250813_ML_allTaxa_MTX_singledb_smaller1000_pred_v1.2.csv', header = T)
all.pred.Val.smaller.1000$BW = "1000"

all.pred.Val.smaller = rbind(all.pred.Val.smaller.750,
                             all.pred.Val.smaller.800,
                             all.pred.Val.smaller.850,
                             all.pred.Val.smaller.900,
                             all.pred.Val.smaller.950,
                             all.pred.Val.smaller.1000)

all.pred.Val.smaller$BW_group = "Smaller"

# larger
all.pred.Val.larger.750 <- read.csv('BW750/250813_ML_allTaxa_MTX_singledb_larger750_pred_v1.2.csv', header = T)
all.pred.Val.larger.750$BW = "750"
all.pred.Val.larger.800 <- read.csv('BW800/250813_ML_allTaxa_MTX_singledb_larger800_pred_v1.2.csv', header = T)
all.pred.Val.larger.800$BW = "800"
all.pred.Val.larger.850 <- read.csv('BW850/250813_ML_allTaxa_MTX_singledb_larger850_pred_v1.2.csv', header = T)
all.pred.Val.larger.850$BW = "850"
all.pred.Val.larger.900 <- read.csv('BW900/250813_ML_allTaxa_MTX_singledb_larger900_pred_v1.2.csv', header = T)
all.pred.Val.larger.900$BW = "900"
all.pred.Val.larger.950 <- read.csv('BW950/250813_ML_allTaxa_MTX_singledb_larger950_pred_v1.2.csv', header = T)
all.pred.Val.larger.950$BW = "950"
all.pred.Val.larger.1000 <- read.csv('BW1000/250813_ML_allTaxa_MTX_singledb_larger1000_pred_v1.2.csv', header = T)
all.pred.Val.larger.1000$BW = "1000"

all.pred.Val.larger = rbind(all.pred.Val.larger.750,
                            all.pred.Val.larger.800,
                            all.pred.Val.larger.850,
                            all.pred.Val.larger.900,
                            all.pred.Val.larger.950,
                            all.pred.Val.larger.1000)

all.pred.Val.larger$BW_group = "Larger"

all.pred.Val2 = rbind(all.pred.Val.smaller,all.pred.Val.larger)

# taxincl.order <- c('No Microbiome Data', 'Including Selected Taxa')
# pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity',
#                'F1','Recall','Precision')
pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity')
fillcols <- c('#7f808e', '#ed4b4d', '#1e8df2','#a394d0',"#338327",'#fcc3c3','#00537a') 
linecols <- c('#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e')
model.order <- c('VirTaxa','BacTaxa','BacRest','BacPath','BacIrep','BacMTX')

# Re-order DataSet (model) and Taxa features given provided arguments
all.pred <- all.pred.Val2 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                               levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)

all.pred.long = filter(all.pred.long, DataSet == "VirTaxa" & !is.na(Performance_Measure_Value))
all.pred.long$BW = as.numeric(all.pred.long$BW)



all.pred.long4 = filter(all.pred.long, Days_before_onset>0 & Days_before_onset<=12)
all.pred.long4$Group[all.pred.long4$Days_before_onset<=12] = "DBO 9-12"
all.pred.long4$Group[all.pred.long4$Days_before_onset<=8] = "DBO 5-8"
all.pred.long4$Group[all.pred.long4$Days_before_onset<=4] = "DBO 1-4"

p_rf_pred_line_vir <- ggplot(all.pred.long4, aes(x=BW, y=Performance_Measure_Value,
                                                    fill=BW_group, color=BW_group))+
  # geom_smooth(span=2, alpha=0.25) + # method = loess
  stat_summary(fun=mean, geom='point') +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data="mean_cl_normal", geom = "errorbar", width = 0.3) +
  # geom_line() +
  # geom_point(alpha=0.5, size = 0.2) +
  facet_grid(Group~Performance_Measure,scales="fixed")+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = fillcols)+
  xlab("BW cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.1,0.75),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5)) +
  stat_compare_means(aes(group = as.factor(BW)), label = "p.signif", vjust = 0.5, ref.group = as.factor(2))
  # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long4$BW)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long4$BW)), 2))))),
  #             map_signif_level = T, textsize=8, step_increase = .01)
  # geom_hline(data=mu, aes(yintercept=grp.mean),
  #            colour = "grey60",
  #            linetype = 2) 

p_rf_pred_line_vir
ggsave("250826_lineVir_Pred_allTaxa_MTX_BW_singledb_DBO_v1.2.1.pdf", p_rf_pred_line_vir , width=8, height=7)


all.pred.long5 = filter(all.pred.long, Days_before_onset>0 & Days_before_onset<=14)
all.pred.long5$Group[all.pred.long5$Days_before_onset<=14] = "DBO 8-14"
all.pred.long5$Group[all.pred.long5$Days_before_onset<=7] = "DBO 1-7"

p_rf_pred_line_vir <- ggplot(all.pred.long5, aes(x=BW, y=Performance_Measure_Value,
                                                 fill=BW_group, color=BW_group))+
  # geom_smooth(span=2, alpha=0.25) + # method = loess
  stat_summary(fun=mean, geom='point') +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data="mean_cl_normal", geom = "errorbar", width = 0.3) +
  # geom_line() +
  # geom_point(alpha=0.5, size = 0.2) +
  facet_grid(Group~Performance_Measure,scales="fixed")+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = fillcols)+
  xlab("BW cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.1,0.6),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))
# stat_compare_means(aes(group = BW_group), label = "p.signif", angle = 0, vjust = 1)
# geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$BW)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$BW)), 2))))),
#             map_signif_level = T, textsize=8, step_increase = .01)
# geom_hline(data=mu, aes(yintercept=grp.mean),
#            colour = "grey60",
#            linetype = 2) 

p_rf_pred_line_vir
ggsave("250826_lineVir_Pred_allTaxa_MTX_BW_singledb_DBO_v1.2.2.pdf", p_rf_pred_line_vir , width=8, height=4.8)

