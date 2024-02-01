#linear model looking at interaction of site and rain on antlers 

library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
data <- read.csv('./clean/morpho_rain.csv', header=T)
head(data)

# #create vector with dmp/rain combo for birthyear of individual
# unique(data$rain.site.by)
# rain.bs.by <- as.numeric(as.factor(data$rain.site.by))

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$birthsite))
unique(data$birthsite)
unique(bs) 
bs  #dmp is 1, e yana 2, w yana 3

#create vectors of things of interest
class(data$year_cap)
data$age <- data$year_cap - data$year_birth 
ageclass <- data$age
unique(ageclass)

#create a vector for capture year, to fit random effect later
capyear <- data$year_cap-min(data$year_cap)+1
unique(capyear)

#create vector for birthyear, so I can fit a random effect later
birthyear<-data$year_birth-min(data$year_birth)+1
unique(birthyear)

#create vector with animal_id for the 475 unique individuals
animal_id <- as.numeric(factor(data$animal_id))
length(unique(animal_id))

#create vector with weight
weightkg <- data$weightkg
weightlb <- data$weight

#create vector with antlers
antlercm<- data$bcscm
antlerin <- data$bcsin

#number of observations 
n <- nrow(data)

#create a vector with rainfall data that is scaled and centered
rain.cy <- scale(data$annual.cy)
rain.cy<- as.vector(rain.cy)

age<- c(1,7,10)

nvalues <- 100
rain.sim <- seq(from = min(rain.cy, na.rm = T), to = max(rain.cy, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
rain.sim

data$rain.by <- as.factor(ifelse(data$annual.by <= 18.3, 1, 
                    ifelse(data$annual.by >18.3 & data$annual.by < 31.8, '2',
                        ifelse(data$annual.by>= 31.8, '3', NA))))

rain.by<- as.numeric(as.factor(data$rain.by))

# #cap year and birth year rain
# rain.site.cy <- as.numeric(as.factor(data$rain.site.cy))
# rain.site.by <- as.numeric(as.factor(data$rain.site.by))

#Write linear mixed effects model for antlers ~ age + site + rain + site*rain
set.seed(100)
sink('ant.bs.rain.jags')
cat('
model{

#priors

for(u in 1:12){ #ageclass
  age.beta[u] ~ dnorm(0,0.001)
}

sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

rain.beta ~ dnorm(0,0.001) #rain

for (u in 1:3){     #birth year rain
  rain.by.beta[u] ~ dnorm(0,0.001)
}
  
for (u in 1:3) { #site
  for (j in 1:3){ # birth year rain
    interaction.beta[u,j] ~ dnorm(0, 0.001)
  }
}

# Likelihood
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- rain.beta*rain[i] + bs.beta[bs[i]] + rain.by.beta[rain.by[i]] + age.beta[ageclass[i]] + 
                                                      interaction.beta[bs[i],rain.by[i]]*rain[i]
}


#derived parameter
  for (i in 1:3){ #birthsite, birth rain
    for (j in 1:100){ #rain sim
        for (k in age){ #ageclasses
      bcs[j,i,i,k] <- rain.beta*rain.sim[j] + bs.beta[i] + rain.by.beta[i] + age.beta[k] + interaction.beta[i,i] * rain.sim[j]
    }
    }
  }
 
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, rain.by = rain.by, antlers = antlerin, rain = rain.cy, age=age, ageclass = ageclass, 
                  rain.sim = rain.sim)

#inits function
inits<- function(){list(bs.beta = rnorm(3, 0, 1),  rain.by.beta = rnorm(3,0,1),
                        rain.beta = rnorm(1,0,1), sigma = rlnorm(1), age.beta = rnorm(12,0,1))}
                         #log normal pulls just positive values,interaction.beta = rnorm(9, 0, 1),

#parameters to estimate
parameters <- c('bs.beta', 'rain.beta', 'rain.by.beta', 'age.beta', 'interaction.beta', 'bcs')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

ant.rain.jags<- jagsUI(jags.data, inits, parameters, 'ant.bs.rain.jags', n.thin = nt, n.chains = nc,
                  n.burnin = nb, n.iter = ni)
print(ant.rain.jags)
# MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output, 
        #then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- ant.rain.jags %>% gather_draws(bcs[rain,site,by.rain, age]) 
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 27000 values of rain 1

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual.cy)) + mean(data$annual.cy)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 27000)
  vector1 <- c(vector1,rep_i)

}

gather$rain1 <- vector1

#subset gather to only include dmp
dmp <- gather %>% filter(site == '1')

dmp$by.rain<-as.factor(dmp$by.rain)

#plot with age groups 1, 7, 10
plot.dmp<- dmp %>% filter(age == '7') %>% 
  ggplot(aes(x=rain1, y=.value, color = by.rain, fill = by.rain)) +
  #facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ age + site + rain*site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.facet.jpg', plot, width = 15, height = 10)


#plot with only age group 7
plot<- gather %>% filter(age == '7') %>%
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ age + site + rain*site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.rain.7.jpg', plot, width = 15, height = 10)


# 
# #Write linear mixed effects model for antlers ~ age + site + rain + site*rain + site*age, lower WAIC than excluding site*age
# set.seed(100)
# sink('ant.bs.age.rain.jags')
# cat('
# model{
#   
# #priors
# 
# for(u in 1:12){ #ageclass
#   age.beta[u] ~ dnorm(0,0.001)
# }
# 
# sigma ~ dunif(0,100)
# tau <- 1/(sigma*sigma)
# 
# for (u in 1:3){ #birthsite
#   bs.beta[u] ~ dnorm(0,0.001)
# }
# 
# rain.beta ~ dnorm(0,0.001) #rain
# 
# for (u in 1:3) { #site and rain interaction
#   rain.bs.beta[u] ~ dnorm(0, 0.001)
# }
# 
# for (u in 1:12) { #site and age interaction
#   age.bs.beta[u] ~ dnorm(0,0.001)
# }
# 
# # Likelihood 
# for (i in 1:n){
#  antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
#  mu[i] <- rain.beta*rain[i] + bs.beta[bs[i]] + age.beta[ageclass[i]] + rain.bs.beta[bs[i]]*rain[i] + age.bs.beta[ageclass[i]]*bs[i]
# }
# 
# #derived parameter
#   for (i in 1:3){ #birthsite
#     for (j in 1:1000){ #rain sim
#         for (k in age){ #ageclasses
#       bcs[j,i,k] <- rain.beta*rain.sim[j] + bs.beta[i] + age.beta[k] + rain.bs.beta[i] * rain.sim[j] + age.bs.beta[k] * i
#     }
#     }
#   }
# 
# }   
# ',fill = TRUE)
# sink()
# 
# #bundle data
# jags.data <- list(n=n, bs = bs, antlers = antlerin, rain = rain.cy, age=age, ageclass = ageclass, rain.sim = rain.sim)
# 
# #inits function
# inits<- function(){list(bs.beta = rnorm(3, 0, 1), rain.bs.beta = rnorm(3, 0, 1), age.bs.beta = rnorm(12,0,1),
#                         rain.beta = rnorm(1,0,1), sigma = rlnorm(1), age.beta = rnorm(12,0,1))}
# #log normal pulls just positive values
# 
# #parameters to estimate
# parameters <- c('bs.beta', 'rain.beta', 'age.beta', 'rain.bs.beta', 'age.bs.beta', 'bcs')#
# 
# #MCMC settings
# ni <- 2000
# nb <- 1000
# nt<- 1
# nc <- 3
# 
# ant.bs.age.rain.jags<- jagsUI(jags.data, inits, parameters, 'ant.bs.age.rain.jags', n.thin = nt, n.chains = nc, 
#                        n.burnin = nb, n.iter = ni)
# print(ant.bs.age.rain.jags)
# 
# 
# #gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output, 
# #then index in the order from output so above was bcs[j,i,k], can rename accordingly
# gather<- ant.bs.age.rain.jags %>% gather_draws(bcs[rain, site, age]) 
# gather$site <- as.factor(gather$site)
# gather$age <- as.factor(gather$age)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather$rain == 2)[1] # 27000 values of rain 1
# 
# # unscale and uncenter rain.sim
# rain.sim1 <- (rain.sim * sd(data$annual.cy)) + mean(data$annual.cy)
# 
# #create vector containing simulated rainfall data but in the format to sync up with gather
# vector1 <- numeric(0)
# rain.sim3 <- for (i in rain.sim1) {
#   rep_i <- rep(i, times = 27000)
#   vector1 <- c(vector1,rep_i)
# 
# }
# 
# gather$rain1 <- vector1
# 
# 
# #plot with age groups 1, 7, 10
# plot<- gather %>%
#   ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
#   facet_wrap(vars(age), nrow = 1)+
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(alpha = .2, labels = c('dmp', 'ey', 'wy')) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(labels = c('dmp', 'ey', 'wy'))+ #color of line but no opacification
#   labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain.cy + age + site + rain*site + age*site")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.5,0.3),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# 
# ggsave('./figures/antler.facet1.jpg', plot, width = 15, height = 10)
# 
# 
# #plot with only age group 7
# plot<- gather %>% filter(age == '7') %>%
#   ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
#   facet_wrap(vars(age), nrow = 1)+
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(alpha = .2, labels = c('dmp', 'ey', 'wy')) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(labels = c('dmp', 'ey', 'wy'))+ #color of line but no opacification
#   labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain + age(7) + site + rain*site + age*site")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.1,0.9),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# 
# ggsave('./figures/antler.rain.7a.jpg', plot, width = 15, height = 10)
# 


#Write linear mixed effects model for antlers ~ rain + age + ELC + ELC*age + ELC*rain
#LP version: bs.beta[bs[i]] + rain.by.beta[rain.by[i]] + 
#           rain.beta * rain.cy[i,j] + age.beta[ageclass[i,j]] +
#           interaction.beta[bs[i], rain.by[j]] * rain.cy[i,j]
# dont concatenate because it takes away from the individual contribution of the two variables to explain the variance 
set.seed(100)
sink('rsbycy.jags')
cat('
model{
  
#priors

age.beta[1] <- 0
for(u in 2:12){                        #ageclass
  age.beta[u] ~ dnorm(0,0.001)
}

sigma ~ dunif(0,100)
tau <- 1/(sigma*sigma)

bs.beta[1] <- 0
for (u in 2:3){                         #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

rain.beta ~ dnorm(0,0.001)              #rain

rain.by.beta[1]<-0
for (u in 2:3){                         # birth year rain
  rain.by.beta ~ dnorm(0,0.001)
}

# for (u in 1:9){
#   interaction.beta ~ dnorm(0, 0.001)
# }
# 
# for(i in 1:3){ #site
#   for(j in 1:3){  #birth rain
#     interaction.beta[i,j] ~ dnorm(0,0.001)
#   }
# }

# for (u in 1:12){      #age site interaction
#   age.site.beta ~ dnorm(0, 0.001)
# }

# Likelihood 
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution, CHANGED TO LOG NORMAL, log link function
 mu[i] <- bs.beta[bs[i]] + age.beta[ageclass[i]] + rain.beta*rain.cy[i] + rain.by.beta[rain.by[i]] #+
                  #interaction.beta[bs[i],rain.by[i]]*rain.cy[i] #+ age.site.beta[ageclass[i]]*bs[i]
}
# 
# #derived parameter
#   for (i in 1:3){ #birthsite
#     for (j in 1:9){ #concat by and cy conditions
#         for (k in age){ #ageclasses
#           for (h in 1:1000){
#       bcs[h,j,i,k] <- rain.beta*rain.sim[h] + bs.beta[i] + age.beta[k] + bs.cy.beta[j] + bs.by.beta[j] + rain.bs.beta[j] * j + age.bs.beta[k] * i
#     }
#     }
#   }
# }
}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs=bs, antlers=antlerin, rain.cy = rain.cy, ageclass=ageclass, rain.by = rain.by )#, age=age, rain.sim=rain.simrain.sim = rain.sim, age=age,

#inits function
inits<- function(){list(bs.beta = c(NA,rnorm(2,0,1)), rain.beta = rnorm(1,0,1), age.beta = c(NA,rnorm(11,0,1)), rain.by.beta = c(NA,rnorm(2,0,1)),
                                  interaction.beta = rnorm(9, 0, 1), age.site.beta = rnorm(12,0,1), sigma = rlnorm(1))}

#log normal pulls just positive values

#parameters to estimate
parameters <- c('rain.beta', 'age.beta', 'bs.beta','rain.by.beta', 'interaction.beta', 'age.site.beta' )

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

rsbycy.jags<- jagsUI(jags.data, inits, parameters, 'rsbycy.jags', n.thin = nt, n.chains = nc, 
                              n.burnin = nb, n.iter = ni)
print(rsbycy.jags)

MCMCtrace(rsbycy.jags)

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output, 
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- rsbycy.jags %>% gather_draws(bcs[rain, conditions, site, age]) 
gather$conditions <- as.factor(gather$conditions)
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 27000 values of rain 1

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual.cy)) + mean(data$annual.cy)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 243000)
  vector1 <- c(vector1,rep_i)

}

gather$rain1 <- vector1


#plot with age groups 1, 7, 10
plot<- gather %>%
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('dmp', 'ey', 'wy')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('dmp', 'ey', 'wy'))+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain.cy + age + site + rain*site + age*site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.facet1.jpg', plot, width = 15, height = 10)


#plot with only age group 7
plot<- gather %>% filter(age == '7') %>%
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('dmp', 'ey', 'wy')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('dmp', 'ey', 'wy'))+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain + age(7) + site + rain*site + age*site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.rain.7a.jpg', plot, width = 15, height = 10)

###############################
#model for  birthsite/birthyear rain and current rain with age as a categorical variable
set.seed(100)
sink('rain.bs.by.jags')
cat('
model{
  
#priors

for(u in 1:12){ #ageclass
  age.beta[u] ~ dnorm(0,0.001)
}

sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

for (u in 1:9){ #birthsite
  rain.bs.by.beta[u] ~ dnorm(0,0.001)
}

rain.beta ~ dnorm(0,0.001) #rain

# 
# for (u in 1:3) { #site and rain interaction
#   rain.bs.beta[u] ~ dnorm(0, 0.001)
# }

# Likelihood 
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- rain.beta*rain[i] + rain.bs.by.beta[rain.bs.by[i]] + age.beta[ageclass[i]] #+ rain.bs.beta[bs[i]]*rain[i]
}


#derived parameter
  for (i in 1:9){ #birthsite
    for (j in 1:1000){ #rain sim
        for (k in age){ #ageclasses
      bcs[j,i,k] <- rain.beta*rain.sim[j] + rain.bs.by.beta[i] + age.beta[k] #+ rain.bs.beta[i] * rain.sim[j]
    }
    }
  }
}
  
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, antlers = antlerin, rain = rain.cy, age=age, rain.sim=rain.sim,
                  ageclass=ageclass, rain.bs.by=rain.bs.by) 

#inits function
inits<- function(){list(rain.bs.by.beta = rnorm(9, 0, 1), age.beta = rnorm(12,0,1), 
                        rain.beta = rnorm(1,0,1), sigma = rlnorm(1))}
#log normal pulls just positive values,ageclass.beta = rnorm(12,0,1), 

#parameters to estimate
parameters <- c('rain.bs.by.beta', 'rain.beta', 'age.beta', 'bcs')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

ant.rain.jags<- jagsUI(jags.data, inits, parameters, 'rain.bs.by.jags', n.thin = nt, n.chains = nc, 
                       n.burnin = nb, n.iter = ni)
print(ant.rain.jags)
# MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output, 
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- ant.rain.jags %>% gather_draws(bcs[rain,site,age]) 
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 81000 values of rain 1

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual.cy)) + mean(data$annual.cy)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim2 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 81000)
  vector1 <- c(vector1,rep_i)
  
}

gather$rain1 <- vector1


#plot with age groups 1, 7, 10
plot<- gather %>%
  ggplot(aes(x=rain1, y=.value, color = conditions, fill = conditions)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain + age + site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.facet1.jpg', plot, width = 15, height = 10)


#plot with only age group 7
plot<- gather %>% filter(age == '7') %>% 
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ age + site + rain")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.rain.bs.by.jpg', plot, width = 15, height = 10)




#model for  birthsite/birthyear rain and current rain with interaction of age
set.seed(100)
sink('rain.bs.by.age.jags')
cat('
model{
  
#priors

for(u in 1:12){ #ageclass
  age.beta[u] ~ dnorm(0,0.001)
}

sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

for (u in 1:9){ #birthsite
  rain.bs.by.beta[u] ~ dnorm(0,0.001)
}

rain.beta ~ dnorm(0,0.001) #rain

 
for (u in 1:9) { #site and rain interaction
  rain.bs.by.age.beta[u] ~ dnorm(0, 0.001)
}

# Likelihood 
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- rain.beta*rain[i] + rain.bs.by.beta[rain.bs.by[i]] + age.beta[ageclass[i]] + rain.bs.by.age.beta[rain.bs.by[i]]*ageclass[i]
}

# 
# #derived parameter
#   for (i in 1:9){ #birthsite
#     for (j in 1:1000){ #rain sim
#         for (k in age){ #ageclasses
#       bcs[j,i,k] <- rain.beta*rain.sim[j] + rain.bs.by.beta[i] + age.beta[k] + rain.bs.by.age.beta[i] * ageclass[k]
#     }
#     }
#   }
}
  
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, antlers = antlerin, rain = rain.cy, age=age, rain.sim=rain.sim,
                  ageclass=ageclass, rain.bs.by=rain.bs.by) 

#inits function
inits<- function(){list(rain.bs.by.beta = rnorm(9, 0, 1), age.beta = rnorm(12,0,1), rain.bs.by.age.beta = rnorm(9,0,1),
                        rain.beta = rnorm(1,0,1), sigma = rlnorm(1))}
#log normal pulls just positive values,ageclass.beta = rnorm(12,0,1), 

#parameters to estimate
parameters <- c('rain.bs.by.beta', 'rain.beta', 'age.beta', 'rain.bs.by.age.beta')#, 'bcs'

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

ant.age.rain.jags<- jagsUI(jags.data, inits, parameters, 'rain.bs.by.age.jags', n.thin = nt, n.chains = nc, 
                       n.burnin = nb, n.iter = ni)
print(ant.age.rain.jags)
# MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output, 
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- ant.rain.jags %>% gather_draws(bcs[rain,site,age]) 
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 81000 values of rain 1

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual.cy)) + mean(data$annual.cy)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim2 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 81000)
  vector1 <- c(vector1,rep_i)
  
}

gather$rain1 <- vector1


#plot with age groups 1, 7, 10
plot<- gather %>%
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ rain + age + site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

# ggsave('./figures/antler.facet.jpg', plot, width = 15, height = 10)


#plot with only age group 7
plot<- gather %>% filter(age == '7') %>% 
  ggplot(aes(x=rain1, y=.value, color = site, fill = site)) +
  facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d()+ #color of line but no opacification
  labs(x = "RAINFALL", y = "ANTLER SCORE (IN)", title = "antlers ~ age + site + rain")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antler.rain.bs.by.jpg', plot, width = 15, height = 10)
