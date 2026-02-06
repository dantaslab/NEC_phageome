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


# taxincl.order <- c('No Microbiome Data', 'Including Selected Taxa')
# pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity',
#                'F1','Recall','Precision')
pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity')
fillcols <- c('#7f808e', '#ed4b4d', '#1e8df2','#a394d0',"#338327",'#fcc3c3','#00537a') 
linecols <- c('#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e','#1f2e2e')
model.order <- c('VirTaxa','BacTaxa','BacRest','BacPath','BacIrep','BacMTX')

# Re-order DataSet (model) and Taxa features given provided arguments
all.pred <- all.pred.Val.smaller %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                               levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)
all.pred.long$BW = as.numeric(all.pred.long$BW)

##########################################################################################################################################
##########################################################################################################################################
all.pred.Val1 = filter(all.pred.long, Days_before_onset>0 & Days_before_onset<=12)
all.pred.Val1$Group[all.pred.Val1$Days_before_onset<=12] = "DBO 9-12"
all.pred.Val1$Group[all.pred.Val1$Days_before_onset<=8] = "DBO 5-8"
all.pred.Val1$Group[all.pred.Val1$Days_before_onset<=4] = "DBO 1-4"

# all.pred.Val5$Group[all.pred.Val5$Days_before_onset<=14] = "DBO 8-14"
# all.pred.Val5$Group[all.pred.Val5$Days_before_onset<=7] = "DBO 1-7"

t_test_result_combine = c()
for (group in unique(all.pred.Val1$Group)) {
  for (bw in unique(all.pred.Val1$BW)) {
    for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
      all.pred.long_bw <- subset(all.pred.Val1, Group == group & BW == bw & Performance_Measure==measure)
      t_test_result <- all.pred.long_bw %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa", p.adjust.method = "BH")
      t_test_result$Group = group
      t_test_result$BW = bw
      t_test_result$Performance_Measure = measure
      t_test_result$mean_Performance_Measure = c(mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacTaxa']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacRest']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacPath']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacIrep']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacMTX']))
      # t_test_result$mean_Performance_Measure = mean(all.pred.long_bw$Performance_Measure_Value)
      t_test_result_combine = rbind(t_test_result_combine,t_test_result)
    }
  }
}
t_test_result_combine$Performance_Measure <- factor(t_test_result_combine$Performance_Measure,
                                                    levels = c("Balanced.Accuracy","Sensitivity","Specificity"))

p_rf_pred_line <- ggplot(all.pred.Val1, aes(x=BW, y=Performance_Measure_Value, 
                                            fill=DataSet, color=DataSet))+
  # geom_smooth(method = loess, alpha=0.25) +
  stat_summary(fun=mean, geom='point')+
  stat_summary(fun=mean, geom='line')+
  # geom_line() +
  # geom_point(alpha=1, size = 2) +
  # facet_wrap(~Performance_Measure, nrow=1, scales="fixed")+
  facet_grid(Group~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = fillcols)+
  xlab("BW cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "right",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = BW, y = mean_Performance_Measure-0.01, label = p.adj.signif, fill=group2, color=group2), 
            size = 5, vjust = 0) 
  # scale_x_reverse(limits = c(14, 0), breaks=c(0:14))
p_rf_pred_line 
ggsave("250826_line_Pred_allTaxa_MTX_BWsmaller_singledb_DBO_v1.2.1.pdf", p_rf_pred_line , width=8, height=7)


##---------------------------------------------------------------------------------------------

## Generate heatmap to show the differences of dataset compared to virome
library(rstatix)
t_test_result_combine$BW = factor(t_test_result_combine$BW)
t_test_result_combine$statistic_rev = -t_test_result_combine$statistic
t_test_result_combine$p.adj.signif[t_test_result_combine$p.adj.signif=="ns"]=""
ggheatmap <- ggplot(t_test_result_combine, aes(BW, group2, fill = statistic_rev))+
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
  coord_fixed()
# Print the heatmap
# print(ggheatmap)

ggheatmap_sig <- ggheatmap +
  geom_text(aes(BW, group2), label = t_test_result_combine$p.adj.signif, color = "black", size = 4) +
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
  xlab("BW cutoff") +
  ylab("") +
  facet_grid(Group~Performance_Measure)
  # facet_wrap(~Performance_Measure, scale="fixed")
# guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
#                              title.position = "top", title.hjust = 0.5))
ggheatmap_sig
ggsave("250826_heatmap_Pred_allTaxa_MTX_BWsmaller_singledb_DBO_v1.2.1.pdf", ggheatmap_sig, width=8, height=8)


##########################################################################################################################################
##########################################################################################################################################

all.pred.Val2 = filter(all.pred.long, Days_before_onset>0)
all.pred.Val2$Group[all.pred.Val2$Days_before_onset<=14] = "DBO 8-14"
all.pred.Val2$Group[all.pred.Val2$Days_before_onset<=7] = "DBO 1-7"

t_test_result_combine = c()
for (group in unique(all.pred.Val2$Group)) {
  for (bw in unique(all.pred.Val2$BW)) {
    for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
      all.pred.long_bw <- subset(all.pred.Val2, Group == group & BW == bw & Performance_Measure==measure)
      t_test_result <- all.pred.long_bw %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa", p.adjust.method = "BH")
      t_test_result$Group = group
      t_test_result$BW = bw
      t_test_result$Performance_Measure = measure
      t_test_result$mean_Performance_Measure = c(mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacTaxa']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacRest']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacPath']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacIrep']),
                                                 mean(all.pred.long_bw$Performance_Measure_Value[all.pred.long_bw$DataSet=='BacMTX']))
      # t_test_result$mean_Performance_Measure = mean(all.pred.long_bw$Performance_Measure_Value)
      t_test_result_combine = rbind(t_test_result_combine,t_test_result)
    }
  }
}
t_test_result_combine$Performance_Measure <- factor(t_test_result_combine$Performance_Measure,
                                                    levels = c("Balanced.Accuracy","Sensitivity","Specificity"))

p_rf_pred_line <- ggplot(all.pred.Val2, aes(x=BW, y=Performance_Measure_Value, 
                                            fill=DataSet, color=DataSet))+
  # geom_smooth(method = loess, alpha=0.25) +
  stat_summary(fun=mean, geom='point')+
  stat_summary(fun=mean, geom='line')+
  # geom_line() +
  # geom_point(alpha=1, size = 2) +
  # facet_wrap(~Performance_Measure, nrow=1, scales="fixed")+
  facet_grid(Group~Performance_Measure)+
  scale_fill_manual(values = fillcols)+
  scale_color_manual(values = fillcols)+
  xlab("BW cutoff") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "right",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = BW, y = mean_Performance_Measure-0.01, label = p.adj.signif, fill=group2, color=group2), 
            size = 5, vjust = 0) 
# scale_x_reverse(limits = c(14, 0), breaks=c(0:14))
p_rf_pred_line 
ggsave("250826_line_Pred_allTaxa_MTX_BWsmaller_singledb_DBO_v1.2.2.pdf", p_rf_pred_line , width=8, height=4)


##---------------------------------------------------------------------------------------------

## Generate heatmap to show the differences of dataset compared to virome
library(rstatix)
t_test_result_combine$BW = factor(t_test_result_combine$BW)
t_test_result_combine$statistic_rev = -t_test_result_combine$statistic
t_test_result_combine$p.adj.signif[t_test_result_combine$p.adj.signif=="ns"]=""
ggheatmap <- ggplot(t_test_result_combine, aes(BW, group2, fill = statistic_rev))+
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
  coord_fixed()
# Print the heatmap
# print(ggheatmap)

ggheatmap_sig <- ggheatmap +
  geom_text(aes(BW, group2), label = t_test_result_combine$p.adj.signif, color = "black", size = 4) +
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
  xlab("BW cutoff") +
  ylab("") +
  facet_grid(Group~Performance_Measure)
# facet_wrap(~Performance_Measure, scale="fixed")
# guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
#                              title.position = "top", title.hjust = 0.5))
ggheatmap_sig
ggsave("250826_heatmap_Pred_allTaxa_MTX_BWsmaller_singledb_DBO_v1.2.2.pdf", ggheatmap_sig, width=8, height=5)

