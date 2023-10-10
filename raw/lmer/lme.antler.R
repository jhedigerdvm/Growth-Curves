library(dplyr)
library(tidyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(measurements)
library(sjPlot)
library(ggeffects) #necessary to make plots with mixed effects models
library(here)


options(rstudio.help.showDataPreview = FALSE)
getOption("rstudio.help.showDataPreview")
#compare the confidence interval 
confint()
data$age <- as.factor(data$age)
data$birthsite<- as.factor(data$birthsite)

#linear mixed effects model antler size as influenced by birth site and age class as fixed effects, animal ID, cap year, and birth year as random effects 
bcs.lmer<- lmer(bcscm ~ birthsite + age + (1|year_birth) + (1|animal_id)+(1|year_cap) - 1, data = data) #
summary(bcs.lmer)

# Extract the prediction data frame
pred.bcs <- ggpredict(bcs.lmer, terms = c("age", 'bs'))  # this gives overall predictions for the model

#removed 1.5 year olds and 11.5+ from graph due to skewness 
#pred.bcs.2<-pred.bcs[-c(1:3, 31:36),]

# make the base plot and save it in the object "plot_base"
plot_base <- 
  ggplot(data = pred.bcs, mapping = aes(x = x, y = predicted, group = group))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.15,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = 'bold', size = 24, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 20, hjust = 0.5),
        axis.text = element_text(face='bold',size = 16))

(bcs_bar<- plot_base + 
    geom_bar(stat = "identity", width = .5, aes(fill = group), 
             position = position_dodge(width = .6)) + 
    labs(x = "AGE CLASS", y = "ANTLER SIZE (IN)", 
         title = "ANTLER SIZE BY AGE CLASS") +
    scale_fill_manual(labels = c('TREATMENT', 'CONTROL', 'TGT'),
                      values = c("#73B2FF", "#73FFDF","#FF73DF"))+
    geom_errorbar(aes(x = x, ymin= conf.low, ymax=conf.high),
                  width=.2, position=position_dodge(width=.6))+
    coord_cartesian(ylim = c(90, NA))
)
