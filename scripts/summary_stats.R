library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

#Weight of the does harvested in fed pastures from 2006 to 2023
harvest.data<- read.csv('./raw/harvest.weights.does.csv', header = TRUE, na.strings=c("","NA")) 
harvest.data<- harvest.data[,-c(7)]
harvest.data1<-harvest.data[harvest.data$live_weight !=0,] #live weight = dressed weight * 1.3
harvest.data1<-harvest.data1[-c(637),-c(3:4)]

harvest.data1<- harvest.data1 %>%  drop_na(age, live_weight)
write.csv(harvest.data1, './clean/harvest_data.csv') #clean harvest data of does with no NAs

table<-harvest.data1 %>%  group_by(age) %>% summarise(min = min(live_weight),
                                                        max = max(live_weight),
                                                        mean = mean(live_weight),
                                                        sd = sd(live_weight),
                                                        total_count = n())

dmp.does<- read.csv('./raw/dmp_sire_doe_data.csv', header = T, na.strings = c('','NA'))
dmp.does<- dmp.does %>%  drop_na(weight,age) %>% filter(sex=='F')
dmp.does<- dmp.does[-c(233),-c(6:15)]
write.csv(dmp.does, './clean/dmp_doe_data.csv') #live weights of all dmp does used from 2007-2021
table2<-dmp.does %>%  group_by(age) %>% summarise(min = min(weight),
                                                      max = max(weight),
                                                      mean = mean(weight),
                                                      sd = sd(weight),
                                                  total_count = n())
dmp.doe.stats<- table2[c(1,4:13,2,3),]
harvest.doe.stats<- table[c(2,6:14, 3:5),]
write.csv(dmp.doe.stats, 'dmp.doe.stats.csv')
write.csv(harvest.doe.stats, 'harvest.doe.stats.csv')

h.data <- harvest.data1 %>% rename_at('live_weight',~'weight')
h.data<- h.data[,-c(2)]
dmp.data<- dmp.does[,-c(2:3)]
h.data<-h.data[,c(1,3,2)]
dmp.data<-dmp.data %>% mutate(group='dmp')
h.data<- h.data %>%  mutate(group='fed')

doe.weights<-rbind(dmp.data,h.data)
unique(doe.weights$age)
doe.weights$age<-as.factor(doe.weights$age)
doe.weights$group<-as.factor(doe.weights$group)
doe.weights<- filter(doe.weights,age %in% c('1.5','2.5','3.5','4.5','5.5','6.5','7.5','8.5','9.5'))
write.csv(doe.weights, './clean/all_doe_weights.csv')
plot_base <- 
  ggplot(data = doe.weights, mapping = aes(x =factor(age, level=c('1.5','2.5','3.5','4.5','5.5','6.5','7.5','8.5','9.5')), 
                                           y = weight, color = group))+
  #scale_color_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_x_discrete(labels = c('1.5','2.5','3.5','4.5','5.5','6.5','7.5','8.5','9.5'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.1,0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', size = 70, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 52, hjust = 0.5),
        axis.text = element_text(face='bold',size = 48),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg


(doe.plot <- plot_base +
  geom_violin(position = position_dodge(width = 0.6), linewidth = 2, trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75))+
  #geom_jitter(height = 0, width = 0.5)+
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  # facet_grid(mineral~., scales = 'free_y')+
  labs(title='DOE WEIGHT BY AGE CLASS', x = 'AGE CLASS', y='WEIGHT (LBS)')  )

ggsave('doe_weights.png', doe.plot, bg='transparent', width = 20, height = 10)



#weight of sires in DMP and how those weights compare with other pastures
sires<- read.csv('./raw/DMP_sire_data.csv', header = T, na.strings = c('','NA'))
sires<- sires %>%  drop_na(weight,bcs,age) %>% filter(sex=='M')
sires<- sires[,-c(4,7)]
write.csv(sires,'./clean/sires.csv')
sires$bcs<-as.numeric(sires$bcs)
sires$weight<-as.numeric(sires$weight)
sire.weight<- sires %>%  group_by(age) %>% summarise(min = min(weight),
                                                    max = max(weight),
                                                    mean = mean(weight),
                                                    sd = sd(weight),
                                                    total_count=n())
sire.bcs<- sires %>%  group_by(age) %>% summarise(min = min(bcs),
                                                     max = max(bcs),
                                                     mean = mean(bcs),
                                                     sd = sd(bcs),
                                                     total_count=n())
write.csv(sire.bcs,'sire.bcs.csv')
write.csv(sire.weight,'sire.weight.csv')
ey <- data %>% filter(data$birthsite=='ey')
ey <- ey %>%  drop_na(weight,bcsin,age)
ey <- ey[,c(2,4,5,9,11)]
ey.weight<- ey %>%  group_by(age) %>% summarise(min = min(weight),
                                                     max = max(weight),
                                                     mean = mean(weight),
                                                     sd = sd(weight),
                                                     total_count=n())
ey.bcs<- ey %>%  group_by(age) %>% summarise(min = min(bcsin),
                                                  max = max(bcsin),
                                                  mean = mean(bcsin),
                                                  sd = sd(bcsin),
                                                  total_count=n())
ey<- ey[,-c(6:8,10:12)]

write.csv(ey, './clean/ey_nofawns.csv')
write.csv(ey.weight,'ey.weight.csv')
write.csv(ey.bcs,'ey.bcs.csv')

wy <- data %>%filter(data$birthsite=='wy')
wy<- wy %>%  drop_na(weight, bcsin,age)
wy<- wy[,c(2,4,5,9,11)]
wy.weight<- wy %>%  group_by(age) %>% summarise(min = min(weight),
                                                max = max(weight),
                                                mean = mean(weight),
                                                sd = sd(weight),
                                                total_count=n())
wy.bcs<- wy %>%  group_by(age) %>% summarise(min = min(bcsin),
                                             max = max(bcsin),
                                             mean = mean(bcsin),
                                             sd = sd(bcsin),
                                             total_count=n())
write.csv(wy.weight,'wy.weight.csv')
write.csv(wy.bcs,'wy.bcs.csv')

ewy<-rbind(wy,ey)
ewy.weight<- ewy %>%  group_by(age) %>% summarise(min = min(weight),
                                                max = max(weight),
                                                mean = mean(weight),
                                                sd = sd(weight),
                                                total_count=n())
ewy.bcs<- ewy %>%  group_by(age) %>% summarise(min = min(bcsin),
                                             max = max(bcsin),
                                             mean = mean(bcsin),
                                             sd = sd(bcsin),
                                             total_count=n())

write.csv(ewy.weight,'ewy.weight.csv')
write.csv(ewy.bcs, 'ewy.bcs.csv')
sire.data <- sires %>% rename_at('year',~'year_cap')
ewy<- ewy %>%  rename_at('birthsite',~'group')
ewy$group<-as.factor(ewy$group)
sire.data <- sires %>% rename_at('bcs',~'bcsin')
sire.data$age <-as.factor(sire.data$age)
sire.data$age <- recode(sire.data$age, '1.5' = '1', '2.5' = '2', '3.5'='3', '4.5' = '4','5.5' = '5','6.5'='6',
                        '7.5'='7','8.5'='8','9.5'='9','10.5'='10')
sire.data$age<-as.factor(sire.data$age)
sire.data$group<-as.factor(sire.data$group)
sire.data<- sire.data %>% mutate(group='dmp')
ey$age <-as.factor(ey$age)
wy$age<-as.factor(wy$age)
ey.data<- ey %>% mutate(group='ey')
wy.data<- wy %>% mutate(group='wy')
ey.data$age<-as.factor(ey.data$age)
ey.data$group<-as.factor(ey.data$group)
wy.data$age<-as.factor(wy.data$age)
wy.data$group <- as.factor(wy.data$group)
ey.data$year_cap<-as.factor(ey.data$year_cap)
wy.data$year_cap<-as.factor(wy.data$year_cap)
sire.data$year<-as.factor(sire.data$year)
sire.data <- sire.data %>% rename_at('year',~'year_cap')
ewy$year_cap<-as.factor(ewy$year_cap)
ewy$age<-as.factor(ewy$age)

all.data<-rbind(sire.data, ewy)
trim.data <- all.data[all.data$age %in% c('3','4','5','6','7'),]
plot_base <- 
  ggplot(data = trim.data, mapping = aes(x =factor(age, level=c('3','4','5','6','7')), 
                                           y = weight, color = group))+
  #scale_color_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_x_discrete(labels = c('3.5','4.5','5.5','6.5','7.5'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.1,0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', size = 70, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 52, hjust = 0.5),
        axis.text = element_text(face='bold',size = 48),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg


(sire.plot <- plot_base +
   geom_violin(position = position_dodge(width = 0.6), linewidth = 1, trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75))+
    #geom_jitter(height = 0, width = 0.5)+
    # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
    # facet_grid(mineral~., scales = 'free_y')+
    labs(title='BUCK WEIGHT BY AGE CLASS', x = 'AGE CLASS', y='WEIGHT (LBS)')  )

ggsave('sire_weights.png', sire.plot, bg='transparent', width = 20, height = 10)



#number of individuals caught as a 5 6 7 8 9 year old? 
data2 <- data %>% filter(data$age %in% c(5,6,7,8,9)) #contains 284 individuals individuals ages 5 to 9
l<-list()
for (i in unique(data2$animal_id)){
  # i<-unique(data2$animal_id)[1] #debugger coge
  ind<- data2 %>% filter(animal_id==i)
  if (all(c(5,6,7,8,9) %in% ind$age)){ #set up logical vectors where all values are true
    l[[i]]<- ind #if true will be saved to a list
  } else{}
}

data3<-bind_rows(l)


#at what age do most individuals max
max.antler.ind<- data %>%  group_by(animal_id) %>%  slice(which.max(bcscm)) #slice only selects rows containing max bcscm for individuals 
max_antlers_by_age <- hist(max.antler.age$age)

plot_base_lme <- 
  ggplot(data = max.antler.ind, aes(x=age, y=))+
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

lme.plot<- plot_base_lme +
  stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                   labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                              '9.5','10.5','11.5','12.5'))+
  labs(x = "AGE CLASS", y = "ANTLER SIZE (IN)", 
       title = "ANTLER SIZE BY AGE CLASS") 

ggsave('antler_in.png', lme.plot, bg='transparent', width = 10, height = 10)