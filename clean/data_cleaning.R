#data cleaning

library(ggplot2)
library(tidyverse)
library(here)


#Load dataset with no fawns, just animals with antlers
data <- read.csv('./raw/bucks_nofawns.csv', header=T)  #does not contain 2022 data
head(data)

#load rainfall data in inches, november to october periods
rain <- read.csv("C:/Users/Joe/Documents/3-R Projects/Survival/cleaned/rainfall.csv", header = T)

merge <- merge(data, rain, by = c('year', 'birthsite'))
data<- merge


#replace zeros with NA
data[data==0] <- NA

write.csv(data, './clean/morpho_rain.csv', row.names = F)


data<- read.csv('./clean/morpho_rain.csv', header = T)
data1 <- data

hist(data$annual)

#determine top 75%, bottom 25% and median annual rainfall in inches
data %>% summarise(q25 = quantile(annual.birthyear, .25),
                   median = median(annual.birthyear),
                   q75 = quantile(annual.birthyear, .75))

#now how to create a column with birthyear rain and site, .x is capyear rain, .y is birthyear rain
merge <- merge(data, rain, by.x = c('birthsite','year_birth'), by.y = c('birthsite','year'))

data<- merge[,c(1:12, 17, 19, 24)]
data <- rename(data, 'annual.capyear' = 'annual.x', 'annual.birthyear'='annual.y', 
               'rain.site.cy' = 'rain.site')


#create a column with categorical rain values, 1 bottom 25%, 2 median, 3 top 75% rainfall
data$rain.cat.by <- ifelse(data$annual.birthyear <= 18.3, 1, 
                           ifelse(data$annual.birthyear >18.3 & data$annual.birthyear < 31.8, '2',
                           ifelse(data$annual.birthyear>= 31.8, '3', NA)))

#concatenate birthsite and categorical rain into one column with 9 levels
data$rain.site.by<- as.factor(paste(data$birthsite, data$rain.cat.by, sep = "_"))
data<- data[,-16]
data <- rename(data, 'annual.cy' = 'annual.capyear', 'annual.by'='annual.birthyear')

write.csv(data, './clean/morpho_rain.csv', row.names = F)

#######333
data <- read.csv('./raw/master_data.csv', header = T)
data1 <- data %>%  filter(Age != '0.5') #contains 2022 data
data2 <- data1[,c(1:4, 8,9,15:19, 22,28:44)]
data3 <- data2
names(data3) <- c('animal_id', 'date', 'year_cap', 'year_birth', 'age', 'weight_lb', 'bs', 'points_typ' )
