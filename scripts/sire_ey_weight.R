#model to determine if DMP sires are lighter than EY deer

library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

data<- read.csv('./clean/dmp_ey_sires.csv',header = T)

#create a vector with 1 for sire, 2 for ey
group<- as.numeric(factor(data$group))
unique(group)

#create vectors of things of interest
ageclass <- as.numeric(factor(data$age))
unique(ageclass)

#create a vector for capture year, to fit random effect later
capyear<- as.numeric(factor(data$year_cap))

#create vector for birthyear, so I can fit a random effect later
birthyear<- as.numeric(factor(data$year_birth))

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


#Write linear mixed effects model for antlers inches,
set.seed(100)
sink('sire_ey.jags')
cat('
model{

#priors
int ~ dunif(0,100)
eps.ageclass[1] <- 0

for(u in 1:2){
  group.beta[u] ~ dnorm(0,tau.group)
}
tau.group <- 1/(sigma*sigma)
sigma ~ dunif(0,10)

for(u in 2:12){
  eps.ageclass[u] ~ dunif(0,10)
}

# Likelihood
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- int + ageclass.beta[ageclass[i]] + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] + eps.id[id[i]] 
 #  +
 #                   
}

#derived parameter
for (j in 1:12){
  antlerin[j] <- int + ageclass.beta[j]
}

for (j in 1:12){ #age
  antler_difference[j] <- antlerin[6]-antlerin[j]
}

}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, 
                  antler = antlerin, capyear = capyear, birthyear = birthyear,id = animal_id, bs = bs)

#inits function
inits<- function(){list(ageclass.beta = c(NA, runif(11,0,100)), int = runif(0,100),
                        sigma = rlnorm(1), eps.capyear = c(NA,rnorm(13,0,1)), sigma.capyear = runif(0,100),
                        eps.birthyear = c(NA, rnorm(13,0,1)), sigma.birthyear = runif(0,100),
                        eps.id = c(NA,rnorm(474,0,1)), sigma.id = runif(0,100),
                        eps.bs = c(NA, rnorm(2,0,1)), sigma.bs = runif(0,100)
) } #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'antlerin', 'antler_difference')

#MCMC settings
ni <- 15000
nb <- 5000
nt<- 10
nc <- 3

lme.ant.jags<- jagsUI(jags.data, inits, parameters, './jags_files/lme.ant.jags', n.thin = nt, n.chains = nc,
                      n.burnin = nb, n.iter = ni)
print(lme.ant.jags)
MCMCtrace(lme.ant.jags)
