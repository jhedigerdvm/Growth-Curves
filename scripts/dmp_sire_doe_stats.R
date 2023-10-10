
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)
library(lme4)
library(lmerTest)

sire.data<-read.csv('./clean/sires.csv')
ey.buck.data<-read.csv('./clean/ey_nofawns.csv')

dmp.ey.sires<- rbind(sire.data,ey.buck.data)
write.csv(dmp.ey.sires, './clean/dmp_ey_sires.csv',row.names = F)

yearcap<-as.numeric(dmp.ey.sires$year_cap)
yearbirth<-as.numeric(dmp.ey.sires$year_birth)
id<-as.factor(dmp.ey.sires$animal_id)
age<-dmp.ey.sires$age
weight<-as.numeric(dmp.ey.sires$weight)
bcs<-dmp.ey.sires$bcs
dmp.ey.sires<-as.matrix(dmp.ey.sires)

#sire,ey,wy, fit model with that and age
ey.sire.data$sire<- as.factor(ey.sire.data$sire) 
ey.sire.data$age <- as.factor(ey.sire.data$age)
ey.sire.data$year_cap <- as.factor(ey.sire.data$year_cap)
ey.sire.data$animal_id <- as.factor(ey.sire.data$animal_id)
lme.weight <- lmer(weight ~ sire + age + (1|year_cap) + (1|animal_id), ey.sire.data)
summary(lme.weight)

#
doe.weights$year<- as.factor(doe.weights$year)
lme.doeweight <- lmer(weight ~ group + age + (1|year), doe.weights)
summary(lme.doeweight)
#Load dataset
data1<- read.csv('./raw/bucks_nofawns_sires.csv')
data1[data1 == ""] <- NA                     # Replace blank by NA
data1$sire<- data1$sire %>% replace_na('no')

#dataframe containing all ey deer
ey<- data1 %>%  filter(data1$birthsite=='ey')
#sire only data
sire <- data1 %>% filter(data1$sire=='yes')

ey.sire.data<- rbind(ey,sire)


#replace zeros with NA
ey.sire.data[ey.sire.data==0] <- NA
ey.sire.data<-ey.sire.data[-c(1013:1016),]

#create a vector with 1 for treatment, 2 for control, 3 for tgt
group<- as.numeric(factor(ey.sire.data$sire))
unique(group) 
group  #no sire 1, yes sire 2

#create vectors of things of interest
class(ey.sire.data$year_cap)
ey.sire.data$age <- ey.sire.data$year_cap - ey.sire.data$year_birth 
ageclass <- ey.sire.data$age
unique(ageclass)

#create a vector for capture year, to fit random effect later
capyear <- ey.sire.data$year_cap-min(ey.sire.data$year_cap)+1
unique(capyear)

# #create vector for birthyear, so I can fit a random effect later
# birthyear<-ey.sire.data$year_birth-min(ey.sire.data$year_birth)+1
# unique(birthyear)

#create vector with animal_id for the 475 unique individuals
animal_id <- as.numeric(factor(ey.sire.data$animal_id))
length(unique(animal_id))

#create vector with weight
weightkg <- ey.sire.data$weightkg
weightlb <- eysire$weight
# 
# #create vector with antlers
# antlercm<- data$bcscm
# antlerin <- data$bcsin

#number of observations 
n <- nrow(ey.sire.data)

# 
# #quick summary stats of data
# dmp<-data %>% filter(data$birthsite=='dmp') 
# ey <- data %>% filter(data$birthsite=='ey') 
# wy <- data %>% filter(data$birthsite=='wy') 
# 
# #assign group two sires
# sire.data<- sire.data
# sire.data$group <- recode(sire.data$group, 'dmp' = 'sire')
# sire.data$animal_id

#Write linear mixed effects model for weight pounds,
set.seed(100)
sink('weightsire.jags')
cat('
model{
  
#priors
for(u in 1:12){
  ageclass.beta[u] ~ dunif(0,100)
}
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each weight is a draw from this distribution
 mu[i] <- ageclass.beta[ageclass[i]] 
}
 
#derived parameter
for (j in 1:12){ #age
  weight_difference[j] <- ageclass.beta[7]-ageclass.beta[j]
}

}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass,
                  weight = weightkg)

#inits function
inits<- function(){list(ageclass.beta = runif(0,100),
                        sigma = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'weight_difference')

#MCMC settings
ni <- 4000
nb <- 1000
nt<- 10
nc <- 3

weight.sire.jags<- jagsUI(jags.data, inits, parameters, 'weightsire.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
MCMCtrace(weight.sire.jags)
print(weight.sire.jags)
write.csv(print(lme.jags$summary), 'weightoutput.csv')
MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results




####################################
#LME bodyweight, birth site, random effects
####################################
set.seed(100)
sink('weightsire.jags')
cat('
model{
  
#priors
  
  ageclass.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate
  # group.beta[1] <- 0     #need to set baseline of zero for >1 categorical covariate
  
  int~dunif(30,50)        #need an intercept to build priors upon

  for (u in 1:2){
    group.beta[u] ~ dunif(0,100) #needed dnorm for bs because of negative values!
  }
  
  for(u in 2:12){
    ageclass.beta[u] ~ dunif(0,150)
  }

  tau <- 1/(sigma*sigma)
  sigma ~ dunif(0,100)


  for(capyear in 1:14){                    #random effect for capyear
    eps.capyear[capyear] ~ dnorm(0, tau.capyear)
  }
  tau.capyear <- 1/(sigma.capyear * sigma.capyear)
  sigma.capyear ~ dunif(0, 100) #random effect SD youre saying that each year-specific error term is coming
                                # from the same distributionm i.e. same mean (0) and standard
                                # deviation (sigma.period). The standard deviation is what connects all the
                                # year-specific random effects to the same distribution.

  # for(birthyear in 1:14){                    #random effect for birthyear
  #     eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
  #   }
  # tau.birthyear <- 1/(sigma.birthyear * sigma.birthyear)
  # sigma.birthyear ~ dunif(0, 150)

# 
#   for(animal_id in 1:475){                       #random effect for animal_id
#     eps.id[animal_id] ~ dnorm(0, tau.id)
#   }
# 
#   tau.id <- 1/(sigma.id * sigma.id)
#   sigma.id ~ dunif(0, 100)

# Likelihood 
for (i in 1:n){
   weight[i] ~ dnorm(mu[i], tau) #each weight is a draw from this distribution
   mu[i] <- int + (ageclass.beta[ageclass[i]])  + group.beta[group[i]] + eps.capyear[capyear[i]] #+ eps.birthyear[birthyear[i]] + eps.id[animal_id[i]]
  }
# 
#derived parameter
# for (i in 1:2){ #group
#   for (j in 1:12){ #age
#     derived_weight[i,j]<- int+ group.beta[i] + ageclass.beta[j]
# }
# }
# 
group.diff <- group.beta[1] - group.beta[2]
  
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n,  weight = weightkg, group=group, ageclass = ageclass, 
                  capyear=capyear)#birthyear = birthyear, animal_id = animal_id

#inits function
inits<- function(){list(int = runif(1,30,50), sigma = rlnorm(1), group.beta = runif(2,0,1),
                        eps.capyear = rnorm(14, 0, 1),sigma.capyear = rlnorm(1), 
                        ageclass.beta = c(NA,runif(11,0,150))
                        
                        
)} #log normal pulls just positive values eps.birthyear = rnorm(14,0,1),sigma.birthyear = rlnorm(1)eps.id = runif(475, 0, 100), sigma.id = rlnorm(1)

#parameters to estimate
parameters <- c('int','group.beta','ageclass.beta', 
                'group.diff')#'derived_weight', 

#MCMC settings
ni <- 20000
nb <- 15000
nt<- 10
nc <- 3

weight.sire.jags<- jagsUI(jags.data, inits, parameters, 'weightsire.jags', n.thin = nt, n.chains = nc, 
                     n.burnin = nb, n.iter = ni)
MCMCtrace(weight.sire.jags)
print(weight.sire.jags)
