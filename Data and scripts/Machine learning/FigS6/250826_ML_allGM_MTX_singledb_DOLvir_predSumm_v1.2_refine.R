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
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/DOL/SingleDB')

# before
all.pred.Val.before.25 <- read.csv('DOL_group_25/250815_ML_allTaxa_MTX_singledb_DOLbefore25_pred_v1.2.csv', header = T)
all.pred.Val.before.25$DOL = "25"
all.pred.Val.before.30 <- read.csv('DOL_group_30/250815_ML_allTaxa_MTX_singledb_DOLbefore30_pred_v1.2.csv', header = T)
all.pred.Val.before.30$DOL = "30"
all.pred.Val.before.35 <- read.csv('DOL_group_35/250815_ML_allTaxa_MTX_singledb_DOLbefore35_pred_v1.2.csv', header = T)
all.pred.Val.before.35$DOL = "35"
all.pred.Val.before.40 <- read.csv('DOL_group_40/250815_ML_allTaxa_MTX_singledb_DOLbefore40_pred_v1.2.csv', header = T)
all.pred.Val.before.40$DOL = "40"
all.pred.Val.before.50 <- read.csv('DOL_group_50/250815_ML_allTaxa_MTX_singledb_DOLbefore50_pred_v1.2.csv', header = T)
all.pred.Val.before.50$DOL = "50"

all.pred.Val.before = rbind(all.pred.Val.before.25,
                            all.pred.Val.before.30,
                            all.pred.Val.before.35,
                            all.pred.Val.before.40,
                            all.pred.Val.before.50)

all.pred.Val.before$DOL_group = "Before"

# after
all.pred.Val.after.25 <- read.csv('DOL_group_25/250815_ML_allTaxa_MTX_singledb_DOLafter25_pred_v1.2.csv', header = T)
all.pred.Val.after.25$DOL = "25"
all.pred.Val.after.30 <- read.csv('DOL_group_30/250815_ML_allTaxa_MTX_singledb_DOLafter30_pred_v1.2.csv', header = T)
all.pred.Val.after.30$DOL = "30"
all.pred.Val.after.35 <- read.csv('DOL_group_35/250815_ML_allTaxa_MTX_singledb_DOLafter35_pred_v1.2.csv', header = T)
all.pred.Val.after.35$DOL = "35"
all.pred.Val.after.40 <- read.csv('DOL_group_40/250823_ML_allTaxa_MTX_singledb_DOLafter40_pred_v1.2.csv', header = T)
all.pred.Val.after.40$DOL = "40"
all.pred.Val.after.50 <- read.csv('DOL_group_50/250823_ML_allTaxa_MTX_singledb_DOLafter50_pred_v1.2.csv', header = T)
all.pred.Val.after.50$DOL = "50"


all.pred.Val.after = rbind(all.pred.Val.after.25,
                           all.pred.Val.after.30,
                           all.pred.Val.after.35,
                           all.pred.Val.after.40,
                           all.pred.Val.after.50)

all.pred.Val.after$DOL_group = "After"

all.pred.Val2 = rbind(all.pred.Val.before,all.pred.Val.after)

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
all.pred.long$DOL = as.numeric(all.pred.long$DOL)



all.pred.long4 = filter(all.pred.long, Days_before_onset>0 & Days_before_onset<=12)
all.pred.long4$Group[all.pred.long4$Days_before_onset<=12] = "DBO 9-12"
all.pred.long4$Group[all.pred.long4$Days_before_onset<=8] = "DBO 5-8"
all.pred.long4$Group[all.pred.long4$Days_before_onset<=4] = "DBO 1-4"

p_rf_pred_line_vir <- ggplot(all.pred.long4, aes(x=DOL, y=Performance_Measure_Value,
                                                    fill=DOL_group, color=DOL_group))+
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
  xlab("DOL cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.1,0.75),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))
  # stat_compare_means(aes(group = DOL_group), label = "p.signif", angle = 0, vjust = 1)
  # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DOL)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DOL)), 2))))),
  #             map_signif_level = T, textsize=8, step_increase = .01)
  # geom_hline(data=mu, aes(yintercept=grp.mean),
  #            colour = "grey60",
  #            linetype = 2) 

p_rf_pred_line_vir
ggsave("250826_lineVir_Pred_allTaxa_MTX_DOL_singledb_DBO_v1.2.1.pdf", p_rf_pred_line_vir , width=8, height=7)

all.pred.long5 = filter(all.pred.long, Days_before_onset>0 & Days_before_onset<=14)
all.pred.long5$Group[all.pred.long5$Days_before_onset<=14] = "DBO 8-14"
all.pred.long5$Group[all.pred.long5$Days_before_onset<=7] = "DBO 1-7"

p_rf_pred_line_vir <- ggplot(all.pred.long5, aes(x=DOL, y=Performance_Measure_Value,
                                                 fill=DOL_group, color=DOL_group))+
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
  xlab("DOL cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.1,0.63),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))
# stat_compare_means(aes(group = DOL_group), label = "p.signif", angle = 0, vjust = 1)
# geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DOL)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DOL)), 2))))),
#             map_signif_level = T, textsize=8, step_increase = .01)
# geom_hline(data=mu, aes(yintercept=grp.mean),
#            colour = "grey60",
#            linetype = 2) 

p_rf_pred_line_vir
ggsave("250826_lineVir_Pred_allTaxa_MTX_DOL_singledb_DBO_v1.2.2.pdf", p_rf_pred_line_vir , width=8, height=4.8)



