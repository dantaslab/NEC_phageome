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
# (e.g. allmeta, Beta ..), and Taxa indicating if feature-selected taxa were 
# included as features in the model, to prediction results (all.pred.Val). 

## Load data
setwd('~/Google Drive/WashU/Phageome/NEC/Manuscript/Figures_drop/Figure 4_ML/Boruta/AddOn_virome/Select_ImpVar')
# all.pred.Val2.1 = read.csv('241201_ML_ImpGM_MTX_late_rest_pred_v1.2re.csv', header = T)
# all.pred.Val2.1 = filter(all.pred.Val2.1, DataSet == "impAllVari")
# all.pred.Val2.2 = read.csv('250417_ML_ImpVIr_late_rest_Test3_re/250417_ML_ImpVIr_late_rest_pred_v1.2.csv', header = T)
# all.pred.Val2.2 = filter(all.pred.Val2.2, DataSet!="impRestVari8" & DataSet!="impRestVari19" & DataSet!="impRestVari3" & DataSet!="impRestVari2")
# all.pred.Val2   = rbind(all.pred.Val2.1, all.pred.Val2.2)
all.pred.Val2 = read.csv('250417_ML_ImpVIr_late_rest_Test3_re/250417_ML_ImpVIr_late_rest_pred_v1.2.csv', header = T)
all.pred.Val2 = filter(all.pred.Val2, DataSet!="impRestVari8" & DataSet!="impRestVari19" & DataSet!="impRestVari3" & DataSet!="impRestVari2")
all.pred.Val2$NEC_onset = "Late"
unique(all.pred.Val2$DataSet)

# taxincl.order = c('No Microbiome Data', 'Including Selected Taxa')
pred_vars = c('Balanced.Accuracy', 'Sensitivity', 'Specificity')
               #'F1','Recall','Precision')
fillcols = c('#1e8df2', '#E2696F', '#97b1a6', '#f4d3d7', '#e59572', '#ee9c83', '#dd5670', '#851321') 
linecols = c('#1f2e2e', '#1f2e2e', '#1f2e2e', '#1f2e2e', '#1f2e2e', '#1f2e2e', '#1f2e2e', '#1f2e2e')
model.order = c('BacRest', 'impBacRest', 'BorutaBacRest', 'MaaslinBacRest', 'impRestVari9')

# Re-order DataSet (model) and Taxa features given provided arguments
all.pred = all.pred.Val2 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long = gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure = factor(all.pred.long$Performance_Measure,
                                               levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value = as.numeric(all.pred.long$Performance_Measure_Value)


##----------------------------------------------------------------------------------------------------
t_test_result_combine = c()
for (tp in 0:14) {
  for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
    all.pred.long_tp <- subset(all.pred.long, Days_before_onset == tp & Performance_Measure==measure)
    t_test_result <- all.pred.long_tp %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "BacRest", p.adjust.method = "BH")
    t_test_result$Days_before_onset = tp
    t_test_result$Performance_Measure = measure
    t_test_result$mean_Performance_Measure = c(#mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='impAllVari']), 
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='impBacRest']), 
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BorutaBacRest']), 
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='MaaslinBacRest']), 
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='impRestVari9']))
    # t_test_result$mean_Performance_Measure = mean(all.pred.long_tp$Performance_Measure_Value)
    t_test_result_combine = rbind(t_test_result_combine,t_test_result)
  }
}
t_test_result_combine$Performance_Measure <- factor(t_test_result_combine$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity","F1","Recall","Precision"))
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
        legend.position = c(0.25,0.3),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = Days_before_onset, y = mean_Performance_Measure-0.01, label = p.adj.signif, fill=group2, color=group2), 
             size = 5, vjust = 0) +
  scale_x_reverse(limits = c(14, 0), breaks=c(0:14))
p_rf_pred_line 
ggsave("250417_line1_Pred_ImpRest_late_rest_DBO_v1.2.pdf", p_rf_pred_line , width=10, height=3)

all.pred.long.rest = filter(all.pred.long, DataSet == "BorutaBacRest" | DataSet == "BacRest" | DataSet == "impRestVari9")  # impRestVari9 or BacRest
t_test_result_combine.rest = filter(t_test_result_combine, group2 == "BorutaBacRest" | group2 == "impRestVari9")
p_rf_pred_line_rest <- ggplot(all.pred.long.rest, aes(x=Days_before_onset, y=Performance_Measure_Value, 
                                                    fill=DataSet, color=DataSet))+
  # geom_smooth(span=2, alpha=0.25) + # method = loess
  stat_summary(fun=mean, geom='point') +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data="mean_cl_normal", geom = "errorbar", width = 0.3) +
  # geom_line() +
  # geom_point(alpha=0.5, size = 0.2) +
  facet_grid(NEC_onset~Performance_Measure, scales="fixed")+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = fillcols[c(1,3,5)])+
  scale_color_manual(values = fillcols[c(1,3,5)])+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.1,0.25),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  xlab("Days before NEC onset") +
  geom_text(data = t_test_result_combine.rest, aes(x = Days_before_onset, y = mean_Performance_Measure-0.01, label = p.adj.signif, fill=group2, color=group2), 
            size = 5, vjust = 0) +
  scale_x_reverse(limits = c(8, 0), breaks=c(0:8))
p_rf_pred_line_rest
ggsave("250417_line2_Pred_ImpRest_late_rest_DBO0-8_v1.2.1.pdf", p_rf_pred_line_rest , width=10, height=2.85)

# https://rpubs.com/nayefahmad/stat-summary
# mean_se( ) is intended for use with stat_summary. It calculates mean and standard error 
# mean_cl_normal( ) is intended for use with stat_summary. It calculates sample mean and lower and upper Gaussian confidence limits based on the t-distribution

## Generate heatmap to show the differences of dataset compared to virome
library(rstatix)
t_test_result_combine.rest$Days_before_onset = factor(t_test_result_combine.rest$Days_before_onset)
t_test_result_combine.rest$statistic_rev = -t_test_result_combine.rest$statistic
t_test_result_combine.rest$p.adj.signif[t_test_result_combine.rest$p.adj.signif=="ns"]=""
ggheatmap <- ggplot(t_test_result_combine.rest, aes(Days_before_onset, group2, fill = statistic_rev))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#006b93", high = "#e32633", mid = "white",
                       midpoint = 0,
                       # limits = c(-8,17), values = rescale(c(-8,0,17)),
                       space = "Lab",
                       name="T-Test\nt-score") + #\n[compare to Case-Control]
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 0, vjust = 1,
                                   size = 10, hjust = 0.5)) +
  facet_wrap(~Performance_Measure, scale="fixed") +
  coord_fixed() + scale_x_discrete(limits = c("8","7","6","5","4","3","2","1","0"))
  #c("14","13","12","11","10","9","8","7","6","5","4","3","2","1","0")) #c("8","7","6","5","4","3","2","1","0"))
# Print the heatmap
# print(ggheatmap)

ggheatmap_sig <- ggheatmap +
  geom_text(aes(Days_before_onset, group2), label = t_test_result_combine.rest$p.adj.signif, color = "black", size = 4) +
  theme(
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(0.6, 0.7),
    legend.direction = "vertical") +
  xlab("Days before NEC onset") +
  ylab("") +
  facet_wrap(~Performance_Measure, scale="fixed")
# guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
#                              title.position = "top", title.hjust = 0.5))
ggheatmap_sig
ggsave("250417_heatmap_Pred_ImpRest_late_rest_DBO0-8_v1.2.pdf", ggheatmap_sig, width=12, height=3)


##----------------------------------------------------------------------------------------------------
all.pred.Val3 = filter(all.pred.Val2, Days_before_onset>0 & Days_before_onset<=5)
all.pred <- all.pred.Val3 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)
library(plyr)
mu <- ddply(subset(all.pred.long, DataSet == "BorutaBacRest"),  # BorutaBacRest or VirTaxa
            "Performance_Measure", summarise, grp.mean=mean(Performance_Measure_Value))

# Plotting function call
p_rf_pred_box <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, 
                                           fill=DataSet, color=DataSet))+
  geom_boxplot(position=position_dodge(width=0.5), outlier.size=0.5)+
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5))+
  # stat_summary(aes(group=Performance_Measure), fun.y=mean, geom="line", colour="green") +
  facet_grid(NEC_onset~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = linecols)+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  # geom_signif(comparisons = list(c('VirTaxa','impAllVari'),c('VirTaxa','impVirTaxa'),c('VirTaxa','BorutaVirTaxa'),c('BorutaVirTaxa','MaaslinVirTaxa'),c('MaaslinVirTaxa','impAllVari9')),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  stat_compare_means(aes(group = DataSet), label = "p.signif", ref.group = "BorutaBacRest", angle = 0) +
  # scale_x_discrete(limits = rev) +
  ylim(0,NA) +   
  coord_flip() +
  geom_hline(data=mu, aes(yintercept=grp.mean),
             colour = "grey60",
             linetype = 2) 
# geom_hline(yintercept = 1,
#            colour = "grey60",
#            linetype = 2)
p_rf_pred_box
ggsave("250417_box_Pred_ImpRest_late_rest_DBO1-5_v1.2_refine.pdf", p_rf_pred_box , width=10, height=2.58)

##----------------------------------------------------------------------------------------------------
all.pred.Val4 = filter(all.pred.Val2, Days_before_onset>0)
all.pred.Val4$Group[all.pred.Val4$Days_before_onset<=5] = "DBO 1-4"
all.pred.Val4$Group[all.pred.Val4$Days_before_onset>5] = "DBO 5-8"

all.pred <- all.pred.Val4 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)

library(plyr)
mu <- ddply(subset(all.pred.long, DataSet == "BorutaBacRest"), 
            .(Performance_Measure,Group), summarise, grp.mean=mean(Performance_Measure_Value))

# Plotting function call
p_rf_pred_box <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, 
                                           fill=DataSet, color=DataSet))+
  geom_boxplot(position=position_dodge(width=0.5), outlier.size=0.5)+
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5))+
  # stat_summary(aes(group=Performance_Measure), fun.y=mean, geom="line", colour="green") +
  facet_grid(Group~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = linecols)+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  # geom_signif(comparisons = split(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  # geom_signif(comparisons = list(c('VirTaxa','Plus BacTaxa'),c('VirTaxa','Plus BacRest'),c('VirTaxa','Plus BacIrep'),c('VirTaxa','Plus BacPath'),c('VirTaxa','Plus BacMAX'),c('VirTaxa','All Features')),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  stat_compare_means(aes(group = DataSet), label = "p.signif", ref.group = "BorutaBacRest", vjust = 0.3) +
  # scale_x_discrete(limits = rev) +
  ylim(0,NA) +   
  coord_flip() +
  geom_hline(data=mu, aes(yintercept=grp.mean),
             colour = "grey60",
             linetype = 2) 

p_rf_pred_box
ggsave("250417_box_Pred_ImpRest_late_rest_DBO>0_group4_v1.2_refine.pdf", p_rf_pred_box , width=10, height=7)




