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
setwd('~/Google Drive/WashU/Phageome/NEC/Olm/Results_SRRs/ML/ImpVar_eachDBO')
all.pred.Val2 <- read.csv('250417_ML_Olm_9vir_pred_v1.2_match_Nosplit_VarByDay.csv', header = T)
unique(all.pred.Val2$DataSet)
# all.pred.Val2 = filter(all.pred.Val2, DataSet!="BacIrep" & Days_before_onset>0)
all.pred.Val2$Test = "NEC early cohort training model"

# taxincl.order <- c('No Microbiome Data', 'Including Selected Taxa')
pred_vars <- c('Balanced.Accuracy', 'Sensitivity', 'Specificity')
               #'F1','Recall','Precision')
fillcols <- c('#7f808e', '#97b1a6', '#e59572') 
linecols <- c('#1f2e2e','#1f2e2e','#1f2e2e')
model.order <- c('VirTaxa', 'BorutaVirTaxa', 'impAllVari') # impAllVari is the 9 viral features

# Re-order DataSet (model) and Taxa features given provided arguments
all.pred <- all.pred.Val2 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                               levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)


##----------------------------------------------------------------------------------------------------
t_test_result_combine = c()
for (tp in 0:8) {
  for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
    all.pred.long_tp <- subset(all.pred.long, Days_before_onset == tp & Performance_Measure==measure)
    t_test_result <- all.pred.long_tp %>% t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa", p.adjust.method = "BH")
    t_test_result$Days_before_onset = tp
    t_test_result$Performance_Measure = measure
    t_test_result$mean_Performance_Measure = c(mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='BorutaVirTaxa']), 
                                               mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='impAllVari']))
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
  stat_summary(fun.data="mean_cl_normal", geom = "errorbar", width = 0.3) +
  # geom_line() +
  # geom_point(alpha=1, size = 2) +
  facet_grid(Test~Performance_Measure, scales="fixed")+
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
        legend.position = c(0.75,0.25),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = Days_before_onset, y = mean_Performance_Measure+0.02, label = p.adj.signif, fill=group2, color=group2), 
             size = 5, vjust = 0) +
  scale_x_reverse(limits = c(8, 0), breaks=c(0:8))
p_rf_pred_line 
ggsave("250417_line_Pred_NEC_Olm_9vir_match_Nosplit_DBO_v1.2.pdf", p_rf_pred_line , width=10, height=3)

##----------------------------------------------------------------------------------------------------
all.pred.long = filter(all.pred.long, DataSet!="BorutaVirTaxa")
t_test_result_combine = c()
for (tp in 0:8) {
  for (measure in c('Balanced.Accuracy', 'Sensitivity', 'Specificity')){
    all.pred.long_tp <- subset(all.pred.long, Days_before_onset == tp & Performance_Measure==measure)
    t_test_result <- all.pred.long_tp %>% pairwise_t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa", p.adjust.method = "BH")
    t_test_result$Days_before_onset = tp
    t_test_result$Performance_Measure = measure
    t_test_result$mean_Performance_Measure = c(mean(all.pred.long_tp$Performance_Measure_Value[all.pred.long_tp$DataSet=='impAllVari']))
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
  stat_summary(fun.data="mean_cl_normal", geom = "errorbar", width = 0.3) +
  # geom_line() +
  # geom_point(alpha=1, size = 2) +
  facet_grid(Test~Performance_Measure, scales="fixed")+
  # facet_grid(Data~Performance_Measure)+
  scale_fill_manual(values = c(fillcols[1],fillcols[3]))+
  scale_color_manual(values = c(fillcols[1],fillcols[3]))+
  xlab("NEC early ML model") +
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = c(0.85,0.2),
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5))+
  geom_text(data = t_test_result_combine, aes(x = Days_before_onset, y = mean_Performance_Measure+0.02, label = p.adj.signif, fill=group2, color=group2), 
            size = 5, vjust = 0) +
  scale_x_reverse(limits = c(2, 0), breaks=c(0:2)) + ylim(0.45,0.7)
p_rf_pred_line 
ggsave("250417_line2_Pred_NEC_Olm_9vir_match_Nosplit_DBO0-2_v1.2.1.pdf", p_rf_pred_line , width=4.7, height=2.7)

##----------------------------------------------------------------------------------------------------
# my_sum <- all.pred.long %>%
#   group_by(DataSet, Performance_Measure, Days_before_onset) %>%
#   summarise( 
#     n=n(),
#     mean=mean(Performance_Measure_Value),
#     sd=sd(Performance_Measure_Value)) %>%
#   mutate( se=sd/sqrt(n))  %>%
#   mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

all.pred.long = filter(all.pred.long, Days_before_onset<=3)
stat.test <- all.pred.long %>%
  group_by(Performance_Measure, Days_before_onset) %>%
  t_test(Performance_Measure_Value ~ DataSet, ref.group = "VirTaxa")

stat.test <- stat.test %>%
  add_y_position(fun = "max", step.increase = 0)

# Plotting function call
p_rf_pred_bar <- 
  ggbarplot(all.pred.long, x="DataSet", y="Performance_Measure_Value", color = "black", fill = "DataSet", 
            add = c("mean_ci"), facet = c("Performance_Measure","Days_before_onset")) +
  scale_fill_manual(values = fillcols)+
  # scale_color_manual(values = fillcols)+
  # facet_wrap(~Performance_Measure, nrow=1, scales="fixed")+
  facet_grid(Days_before_onset~Performance_Measure, scales="fixed")+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none",
        panel.border = element_rect(fill=NA, colour="black", linewidth=0.5)) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj.signif",
                     vjust = 0.8, y=0.65, hide.ns = TRUE,
                     remove.bracket = TRUE, color="black") +
  # scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  # ylim(0,NA) +   
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))
  # coord_flip()
p_rf_pred_bar
ggsave("250417_bar_Pred_NEC_Olm_9vir_match_Nosplit_DBO_v1.2.pdf", p_rf_pred_bar , width=8, height=6)


##----------------------------------------------------------------------------------------------------
all.pred.Val3 = filter(all.pred.Val2, Days_before_onset>0)
all.pred <- all.pred.Val3 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)
library(plyr)
mu <- ddply(subset(all.pred.long, DataSet == "VirTaxa"), 
            "Performance_Measure", summarise, grp.mean=mean(Performance_Measure_Value))

# Plotting function call
p_rf_pred_box <- ggplot(all.pred.long, aes(x=DataSet, y=Performance_Measure_Value, 
                                           fill=DataSet, color=DataSet))+
  geom_boxplot(position=position_dodge(width=0.5), outlier.size=0.5)+
  stat_summary(fun=mean, geom='point', shape=4, position=position_dodge(width=0.5))+
  # stat_summary(aes(group=Performance_Measure), fun.y=mean, geom="line", colour="green") +
  facet_wrap(~Performance_Measure, nrow =1)+
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
  # geom_signif(comparisons = Nosplit(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  # geom_signif(comparisons = list(c('VirTaxa','Plus BacTaxa'),c('VirTaxa','Plus BacRest'),c('VirTaxa','Plus BacIrep'),c('VirTaxa','Plus BacPath'),c('VirTaxa','Plus BacMAX'),c('VirTaxa','All Features')),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  stat_compare_means(aes(group = DataSet), label = "p.signif", ref.group = "VirTaxa", angle = 0) +
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
ggsave("250417_box_Pred_NEC_Olm_9vir_match_Nosplit_DBO>0_v1.2.pdf", p_rf_pred_box , width=10, height=2)

##----------------------------------------------------------------------------------------------------
all.pred.Val4 = filter(all.pred.Val2, Days_before_onset>0)
all.pred.Val4$Group[all.pred.Val4$Days_before_onset<=4] = "DBO 1-4"
all.pred.Val4$Group[all.pred.Val4$Days_before_onset>4] = "DBO 5-8"

all.pred <- all.pred.Val4 %>% mutate(DataSet = fct_relevel(Data, model.order))

# Gather performance metrics of interest and fix variable types.
all.pred.long <- gather(all.pred, 'Performance_Measure', 
                        'Performance_Measure_Value', pred_vars)

all.pred.long$Performance_Measure <- factor(all.pred.long$Performance_Measure,
                                            levels = c("Balanced.Accuracy","Sensitivity","Specificity"))
all.pred.long$Performance_Measure_Value <- as.numeric(all.pred.long$Performance_Measure_Value)

library(plyr)
mu <- ddply(subset(all.pred.long, DataSet == "VirTaxa"), 
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
  # geom_signif(comparisons = Nosplit(t(combn(levels(as.factor(all.pred.long$DataSet)), 2)), seq(nrow(t(combn(levels(as.factor(all.pred.long$DataSet)), 2))))),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  # geom_signif(comparisons = list(c('VirTaxa','Plus BacTaxa'),c('VirTaxa','Plus BacRest'),c('VirTaxa','Plus BacIrep'),c('VirTaxa','Plus BacPath'),c('VirTaxa','Plus BacMAX'),c('VirTaxa','All Features')),
  #             map_signif_level = T, textsize=3, step_increase = .1) +
  stat_compare_means(aes(group = DataSet), label = "p.signif", ref.group = "VirTaxa", vjust = 0.3) +
  # scale_x_discrete(limits = rev) +
  ylim(0,NA) +   
  coord_flip() +
  geom_hline(data=mu, aes(yintercept=grp.mean),
             colour = "grey60",
             linetype = 2) 

p_rf_pred_box
ggsave("250417_box_Pred_NEC_Olm_9vir_match_Nosplit_DBO>0_v1.2.pdf", p_rf_pred_box , width=10, height=3)




