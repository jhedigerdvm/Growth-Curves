#linear model looking at age birthsite interaction on antlers
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
data <- read.csv('./raw/bucks_nofawns.csv', header=T)
head(data)
# 
# data1<- read.csv('./raw/bucks_nofawns_sires.csv')
# data1[data1 == ""] <- NA                     # Replace blank by NA
# data1$sire<- data1$sire %>% replace_na('no')
# 
# #dataframe containing all ey deer
# ey<- data1 %>%  filter(data1$birthsite=='ey')
# #sire only data
# sire <- data1 %>% filter(data1$sire=='yes')
# 
# ey.sire.data<- rbind(ey,sire)


#replace zeros with NA
data[data==0] <- NA

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


#Write linear mixed effects model for antlers in, this model continues increasing linearly into the older age classes
set.seed(100)
sink('ant.bs.age.jags')
cat('
model{
  
#priors
for(u in 1:12){ #ageclass
  ageclass.beta[u] ~ dunif(0,10)
}

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:3) { #site and age interaction
  age.bs.beta[u] ~ dnorm(0, 0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- ageclass.beta[ageclass[i]] + bs.beta[bs[i]] + age.bs.beta[bs[i]] * ageclass[i]
}

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:12) { #ageclass
      antlers[i,j] <- ageclass.beta[j] + bs.beta[i] + age.bs.beta[i]*j
    }
  }

# 
# for (j in 1:12){ #age
#   antler_difference[j] <- ageclass.beta[6]-ageclass.beta[j]
# }

}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, bs = bs, antler = antlerin)

#inits function
inits<- function(){list(ageclass.beta = runif(0,10), bs.beta = rnorm(3, 0, 1), age.bs.beta = rnorm(3, 0, 1),
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('ageclass.beta','bs.beta', 'age.bs.beta', 'antlers')#'antler_difference'

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'ant.bs.age.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)