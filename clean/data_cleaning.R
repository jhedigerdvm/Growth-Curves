#data cleaning

library(ggplot2)
library(tidyverse)
library(dplyr)

library(here)


data <- read.csv('./raw/master_data_2022.csv', header = T) #contains fawn data

#create birthsite column with three groups

data1 <- data %>%
  mutate(
    birthsite = case_when(
      substr(Unique.Buck.Code.., 1, 2) == "90" ~ "ey",
      substr(Unique.Buck.Code.., 1, 2) == "27" & substr(Unique.Buck.Code.., 10, 10) == "0" ~ "wy",
      substr(Unique.Buck.Code.., 1, 2) == "27" & substr(Unique.Buck.Code.., 10, 10) == "1" ~ "dmp",
      TRUE ~ NA_character_  # fallback for other cases
    )
  )

data1 <- data1[,c(1,45, 2:44)] #reorder dataframe bring birthsite to column 2
age_summary <- data1 %>%
  group_by(Capture.Year, birthsite, Age) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Capture.Year, birthsite, Age)

write.csv(age_summary, './output/age_n_byyear.csv', row.names = F)

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
data2 <- data1[,c(1:4, 8,9,15:19, 22,28:42)]
data3 <- data2
names(data3) <- c('animal_id', 'date', 'year_cap', 'year_birth', 'age', 'weight_lb', 'bs', 'points_typ', 
                  'points_total', 'mass', 'bcs_in', 'spike', 'gray_face', 'roman_nose', 'wrinkled_ears', 
                  'loose_face', 'loose_neck', 'loose_chest', 'circ_face_eyes', 'circ_nosebridge', 'eyes_nose',
                  'ear_ear', 'stomach', 'neck', 'tarsal_height', 'tarsal_width', 'stain')
data3[data3 == ""] <- NA

#need to update birth year based upon animal-id
for (i in nrow:data3){
  
  
}



#update with 2022 individuals
master <- read.csv("./raw/master_data_2022.csv", header = T)

master22<- master %>%  filter(Capture.Year== "2022")
master22<- master22[,c(1,3,4,8, 9,19)]
names(master22)<- c('animal_id', 'year_cap','year_birth', 'age', 'weight','bcsin')


#add birth site

for (i in 1:nrow(master22)) {
  x <- substr(master22$animal_id[i], 1, 1)
  if (x == "2") {master22$birthsite[i] <- "wy"} else 
  { master22$birthsite[i] <- "ey"}
}
##make a new column with discernment of DMP and pasture born
for (i in 1:nrow(master22)) {
  x <- substr(master22$animal_id[i], 10, 10)
  if (x == "1") {master22$birthsite[i] <- "dmp"}
}

master22<- master22[master22$age != '0.5',]


merge<- rbind(master22, data)

write.csv(merge, './clean/nofawns22.csv', row.names = F)

data <- read.csv('./clean/nofawns22.csv', header = T)
rain <- read.csv('./clean/rainfall_2022.csv', header = T)

#now how to create a column with birthyear rain and site, .x is capyear rain, .y is birthyear rain
merge <- merge(data, rain, by.x = c('birthsite','year_birth'), by.y = c('birthsite','year'))
merge <- rename(merge, 'annual.by'='annual')
merge <- merge[,-c(8:11)]

#add a column for capture year rain
merge1<- merge(merge, rain, by.x=c('birthsite','year_cap'), by.y=c('birthsite', 'year'))
merge1 <- rename(merge1, 'annual.cy'= 'annual')
merge1 <- merge1[,-c(9:12)] #remove rainvalues other than annual

#rename columns
merge1<- rename(merge1, 'weight.lb'='weight')
merge1<- rename(merge1, 'antler.in'='bcsin')

#create columns for metric system
merge1$antler.cm <- merge1$antler.in*2.54
merge1$weight.kg <- merge1$weight.lb/2.2

#reorder dataframe columns
merge1<- merge1[,c(4,1,3,2,5,6,11,7,10,8,9)]

#round values to two decimals
merge1$weight.kg <- round(merge1$weight.kg, digits = 2)
merge1$antler.cm <- round(merge1$antler.cm, digits = 2)

#save updated file with 2007 to 2022 buck data and rainfall, no fawns
write.csv(merge1, './clean/nofawns22rain.csv', row.names = F)

#read in data if not already
#going to create another column in nofawns22rain.csv to include environmental matching
#this new column will be a 3 class categorical, byrain=cyrain, byrain>cyrain, byrain<cyrain
#we will bin into dry, intermediate and wet
#first lets find the upper and lower quantiles
quantile(data$annual.by, probs = c(0.25, 0.75)) #lower 18.3, upper 31.8
quantile(data$annual.cy, probs = c(0.25, 0.75)) #lower 19.5, upper 29.10

#create new column for binned data with assigned values
data1 <- data %>%
  mutate(by.rain.bin = case_when(
    annual.by <= 18.3 ~ "dry",
    annual.by >= 31.8 ~ "wet",
    TRUE ~ "intermediate"
  ))
data1 <- data1 %>%
  mutate(cy.rain.bin = case_when(
    annual.cy <= 19.5 ~ "dry",
    annual.cy >= 29.10 ~ "wet",
    TRUE ~ "intermediate"
  ))
data2 <- data1 %>%
  mutate(match = case_when(
    cy.rain.bin == by.rain.bin ~ 1,
    by.rain.bin == "wet" & cy.rain.bin %in% c("dry", "intermediate") ~ 2,
    by.rain.bin == "dry" & cy.rain.bin %in% c("intermediate", "wet") ~ 3,
    by.rain.bin == "intermediate" & cy.rain.bin %in% c("dry") ~ 2,
    by.rain.bin == "intermediate" & cy.rain.bin %in% c("wet") ~ 3,
    TRUE ~ NA_real_  # For any unexpected combinations
  ))

data2 %>% count(birthsite,match)
#1 = by and cy rain match, 2 = by rain > cy rain, 3 = by rain < cy rain
data2$match <- as.factor(data2$match)
str(data2)

write.csv(data2, './clean/nofawns22rain.bin.csv', row.names = F)
data<-read.csv('./clean/nofawns22rain.bin.csv', header = T)
