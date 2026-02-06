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
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/Each_dataset')
all.pred.Val2 <- read.csv('241127_ML_allTaxa_MTX_singledb_pred_v1.2_BI.csv', header = T)
# all.pred.Val2 = filter(all.pred.Val2, DataSet!="BacIrep" & Days_before_onset>0)
all.pred.Val2 = filter(all.pred.Val2, Days_before_onset>0)

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

mu <- ddply(subset(all.pred.long, DataSet == "VirTaxa"), 
            "Performance_Measure", summarise, grp.mean=mean(Performance_Measure_Value))

# Plotting function call
p_rf_pred_box <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, fill=DataSet))+
  geom_boxplot(position=position_dodge(width=0.5), outlier.size=0.5)+
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5))+
  facet_wrap(~Performance_Measure, nrow =1)+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = linecols)+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5)) +
  # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  # geom_signif(comparisons = list(c('VirTaxa','BacTaxa'),c('VirTaxa','BacRest'),c('VirTaxa','BacIrep'),c('VirTaxa','BacPath'),c('VirTaxa','BacMTX')),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  stat_compare_means(aes(group = DataSet), label = "p.signif", ref.group = "VirTaxa", angle = 0) +
  # scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0,NA),breaks=seq(0, 1, by = 0.25)) +
  # ylim(0,NA) +   
  coord_flip() +
  geom_hline(data=mu, aes(yintercept=grp.mean),
             colour = "grey60",
             linetype = 2) 
p_rf_pred_box
ggsave("241127_box_Pred_allTaxa_MTX_singledb_DBO>0_v1.2_BI.pdf", p_rf_pred_box , width=10, height=3)


# all.pred.long.grouped <- all.pred.long %>%
#   group_by(Days_before_onset, Performance_Measure, DataSet) %>%
#   summarize(mean_Performance_Measure_Value = mean(Performance_Measure_Value),
#             std_Performance_Measure_Value = sd(Performance_Measure_Value))
# 
# stat.test <- all.pred.long %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "BacTaxa") 

t_test_result_combine = c()
for (tp in 0:8) {
  for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
    all.pred.long_tp <- subset(all.pred.long, Days_before_onset == tp & Performance_Measure==measure)
    t_test_result <- all.pred.long_tp %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa", p.adjust.method = "BH")
    t_test_result$Days_before_onset = tp
    t_test_result$Performance_Measure = measure
    t_test_result$mean_Performance_Measure = c(mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BacTaxa']),
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BacRest']),
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BacIrep']),
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BacPath']),
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BacMTX']))
    # t_test_result$mean_Performance_Measure = mean(all.pred.long_tp$Performance_Measure_Value)
    t_test_result_combine = rbind(t_test_result_combine,t_test_result)
    }
}
t_test_result_combine$Performance_Measure <- factor(t_test_result_combine$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
p_rf_pred_line <- ggplot(all.pred.long, aes(x=Days_before_onset, y=Performance_Measure_Value, 
                                       fill=DataSet, color=DataSet))+
  # geom_smooth(method = loess, alpha=0.25) +
  stat_summary(fun=mean, geom='point')+
  stat_summary(fun=mean, geom='line')+
  # geom_line() +
  # geom_point(alpha=1, size = 2) +
  facet_wrap(~Performance_Measure, nrow=1, scales="fixed")+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = fillcols)+
  xlab("Days before NEC onset") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "right",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = Days_before_onset, y = mean_Performance_Measure-0.01, label = p.adj.signif, fill=group2, color=group2), 
             size = 5, vjust = 0) +
  scale_x_reverse(limits = c(8, 0), breaks=c(0:8))
p_rf_pred_line 
ggsave("241127_line_Pred_allTaxa_MTX_singledb_DBO_v1.2_BI.pdf", p_rf_pred_line , width=12, height=3)
