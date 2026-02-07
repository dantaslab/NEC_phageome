library(ggplot2)
library(dplyr) 
library(forcats)
library(tidyr)
library(ggpubr) 
library(rstatix)

## Plot predictive performance metrics------------------------------------------ 

# We'd like to plot the distributions of accuracy, sensitivity, and specificity
# of predictions on the validation cohort for each of the trained models. Also
# plot null distributions of these parameters from the class label shuffling. 

# For ease of faceting in plots, add columns DataSet indicating the base model
# (e.g. allmeta, Beta ..), and Taxa indicating if feature-earlyed taxa were 
# included as features in the model, to prediction results (all.pred.Val). 

## Load data
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/Each_dataset')
all.pred.Val2.early <- read.csv('250424_ML_allTaxa_MTX_singledb_early_pred_v1.2.csv', header = T)
all.pred.Val2.early$NEC_onset = "Early"
# all.pred.Val2.late <- read.csv('241128_ML_allTaxa_MTX_singledb_late_pred_v1.2.csv', header = T)

all.pred.Val2.late = read.csv('250424_ML_allTaxa_MTX_singledb_late_pred_v1.2.csv', header = T)
all.pred.Val2.late$NEC_onset = "Late"

# all.pred.Val2.all <- read.csv('241127_ML_allTaxa_MTX_singledb_pred_v1.2.csv', header = T)
# all.pred.Val2.all$NEC_onset = "All"
all.pred.Val2 = rbind(all.pred.Val2.early, all.pred.Val2.late)
all.pred.Val2 = filter(all.pred.Val2, Days_before_onset>0)

# taxincl.order <- c('No Microbiome Data', 'Including earlyed Taxa')
pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity')
               #'F1','Recall','Precision')
fillcols <- c('#a32b32','#d4edf8') # c('#e2b96e', '#a32b32','#d4edf8') 
linecols <- c('#1f2e2e','#1f2e2e','#1f2e2e')
model.order <- c('VirTaxa', 'BacTaxa', 'BacRest', 'BacPath', 'BacIrep', 'BacMTX')

# Re-order DataSet (model) and Taxa features given provided arguments
all.pred <- all.pred.Val2 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                               levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)


# stat.test <- all.pred.long %>%
#   group_by(DataSet, Performance_Measure) %>%
#   t_test(Performance_Measure_Value ~ NEC_onset, ref.group = "All")
# 
# stat.test <- stat.test %>%
#   add_x_position(x = "DataSet", dodge = 0.8)

# Plotting function call
p_rf_pred_box <- ggboxplot(all.pred.long, x="DataSet", y="Performance_Measure_Value", 
                           fill = "NEC_onset", 
                           facet.by = "Performance_Measure",
                           outlier.size=0.5)+
  stat_summary(fun = mean, geom = "point", shape=4, aes(group = NEC_onset), position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = fillcols)+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.05,0.8),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  # geom_signif(comparisons =list(c("Early", "Late")), map_signif_level = T, textsize=3) +Ç
  # stat_pvalue_manual(stat.test, label = "p.adj.signif",
  #                    vjust = 1,y=1.1, #hide.ns = TRUE,
  #                    bracket.size = 0) +
  stat_compare_means(
    aes(DataSet = NEC_onset), label = "p.signif", # p.format or p.signif
    method = "wilcox.test", label.y = max(all.pred.long$Performance_Measure_Value) * 1.05, size = 3) +
  # ylim(0,NA) +   
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  coord_flip()
p_rf_pred_box
ggsave("250424_box_Pred_onset_allGM_MTX_singledb_DBO>0_v1.2.pdf", p_rf_pred_box , width=10, height=3)

# p_rf_pred_box_DBO <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, 
#                                            fill=NEC_onset, color=NEC_onset))+
#   geom_boxplot(position=position_dodge(width=0.8), outlier.size=0.5)+
#   stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.8))+
#   # facet_wrap(~Performance_Measure, nrow =2)+
#   facet_grid(Days_before_onset~Performance_Measure)+
#   scale_fill_manual(values = fillcols)+
#   scale_color_manual(values = linecols)+
#   theme_classic()+
#   theme(axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         plot.title = element_text(size=12),
#         legend.position = "right",
#         panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
#   # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
#   #             map_signif_level = T, textsize=3, step_increase = .1) +
#   stat_compare_means(method = "wilcox.test", label ="p.signif", label.y = 0.95) + #p.format p.signif
#   scale_x_discrete(limits = rev)+
#   ylim(0,NA)+   
#   coord_flip()
# p_rf_pred_box_DBO
# ggsave("241111_box_Pred_all_vs_early_Taxa_eachDBO_v1.2_class.pdf", p_rf_pred_box_DBO, width=15, height=15)

