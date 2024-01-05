#linear model looking at interaction of site and rain on antlers 

library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
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

#create a vector with rainfall data
rain <- scale(data$annual)
rain<- as.vector(rain)

age<- c(1,7,10)

nvalues <- 1000
rain.sim <- seq(from = min(rain, na.rm = T), to = max(rain, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
rain.sim


#Write linear mixed effects model for interaction between antlers and rain with age as a categorical variable
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

for (u in 1:3) { #site and rain interaction
  rain.bs.beta[u] ~ dnorm(0, 0.001)
}

# Likelihood 
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- rain.beta*rain[i] + bs.beta[bs[i]] + age.beta[ageclass[i]] + rain.bs.beta[bs[i]]*rain[i]
}

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:1000){ #rain sim
        for (k in age){ #ageclasses
      bcs[j,i,k] <- rain.beta*rain.sim[j] + bs.beta[i] + age.beta[k] + rain.bs.beta[i] * rain.sim[j]
    }
    }
  }


# 
# for (i in 1:3){
#   for (j in 1:10){ #age
#     antlers_diff[i,j] <- antlers[1,j] - antlers[i,j]
#   }
#   }
}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, antlers = antlerin, rain = rain, age=age, ageclass = ageclass, rain.sim = rain.sim)#ageclass = ageclass, 

#inits function
inits<- function(){list(bs.beta = rnorm(3, 0, 1), rain.bs.beta = rnorm(3, 0, 1), 
                        rain.beta = rnorm(1,0,1), sigma = rlnorm(1), age.beta = rnorm(12,0,1))}
                         #log normal pulls just positive values,ageclass.beta = rnorm(12,0,1), 

#parameters to estimate
parameters <- c('bs.beta', 'rain.beta', 'age.beta', 'rain.bs.beta', 'bcs')#

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
gather<- ant.rain.jags %>% gather_draws(bcs[rain,site, age]) gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 27000 values of rain 1

# unscale and uncenter rain.sim
rain.sim1 <- (rain.sim * sd(data$annual)) + mean(data$annual)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim1) {
  rep_i <- rep(i, times = 27000)
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


