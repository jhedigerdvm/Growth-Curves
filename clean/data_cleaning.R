#data cleaning

library(ggplot2)
library(tidyverse)
library(here)


#Load dataset with no fawns, just animals with antlers
data <- read.csv('./raw/bucks_nofawns.csv', header=T)
head(data)
data$year <- data$year_cap
#load rainfall data in inches, november to october periods
rain <- read.csv("C:/Users/Joe/Documents/R Projects/Survival/cleaned/rainfall.csv", header = T)

merge <- merge(data, rain, by = c('year', 'birthsite'))
data<- merge


#replace zeros with NA
data[data==0] <- NA

write.csv(data, './clean/morpho_rain.csv', row.names = F)


data<- read.csv('./clean/morpho_rain.csv', header = T)
data1 <- data

hist(data$annual)

#determine top 75%, bottom 25% and median annual rainfall in inches
data %>% summarise(q25 = quantile(annual, .25),
                   median = median(annual),
                   q75 = quantile(annual, .75))

#create a column with categorical rain values, 1 bottom 25%, 2 median, 3 top 75% rainfall
data1$rain.cat <- ifelse(data1$annual <= 19.4, 1, ifelse(data1$annual >19.4 & data1$annual < 29.1, '2',
                                                      ifelse(data1$annual>= 29.1, '3', NA)))

#concatenate birthsite and categorical rain into one column with 9 levels
data1$rain.site<- as.factor(paste(data1$birthsite, data1$rain.cat, sep = "_"))


