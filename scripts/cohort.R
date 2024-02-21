#cohort effects
#linear model looking at the effect of age and site interaction on antlers
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
data <- read.csv('./raw/bucks_nofawns.csv', header=T)
head(data)

#replace zeros with NA
data[data==0] <- NA

#filter for just 6 year olds
data <- data %>% filter(age=="6")

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


#model evaluating birth year effect on antler scores in 6 year old ageclass
set.seed(100)
sink('ant.cohort.jags')
cat('
model{
  
#priors
# 
# for (u in 1:3){ #birthsite
#   bs.beta[u] ~ dnorm(0,0.001)
# }
# 
# for (u in 1:12) { #site and age interaction
#   age.bs.beta[u] ~ dnorm(0, 0.001)
# }

for (u in 1:9) { #birth year effect
  birthyear.beta[u] ~ dnorm(0,0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:9) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }


# for (u in 1:475){ #random effect for animal id
#   eps.id[u] ~ dnorm(0,sigma)
# }


# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- birthyear.beta[birthyear[i]] + eps.capyear[capyear[i]] 
 
}

#derived parameter
  for (i in 1:9){ #birthsite
      antlers[i] <- birthyear.beta[i]
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
jags.data <- list(n=n,  antler = antlerin, birthyear = birthyear, capyear = capyear) #bs = bs,capyear = capyear, 

#inits function
inits<- function(){list(birthyear.beta = rnorm(9, 0,1), eps.capyear = rnorm(9,0,1), 
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('birthyear.beta', 'eps.capyear','antlers')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'ant.cohort.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(antlers[birthyear]) #this creates a dataframe in long format with indexing
gather$birthyear <- as.factor(gather$birthyear)


plot2<- gather %>% 
  ggplot(aes(x=birthyear, y=.value)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels=c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015'))+
  # scale_fill_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  # scale_color_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  # 
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "YEAR", y = "ANTLER SCORE (IN)", title = "ANTLERS BY YEAR (6 Y.O.)")+
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
plot2
ggsave('./figures/ant.cohort.jpg', plot2, width = 12, height = 6)



#model evaluating birth year effect on antler scores in 6 year old ageclass
set.seed(100)
sink('weight.cohort.jags')
cat('
model{
  
#priors
# 
# for (u in 1:3){ #birthsite
#   bs.beta[u] ~ dnorm(0,0.001)
# }
# 


for (u in 1:9) { #birth year effect
  birthyear.beta[u] ~ dnorm(0,0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:9) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }



# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- birthyear.beta[birthyear[i]] + eps.capyear[capyear[i]] 
 
}

#derived parameter
  for (i in 1:9){ #birthsite
      bodymass[i] <- birthyear.beta[i]
    }

}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n,  weight = weightlb, birthyear = birthyear, capyear = capyear) #bs = bs,capyear = capyear, 

#inits function
inits<- function(){list(birthyear.beta = rnorm(9, 0,1), eps.capyear = rnorm(9,0,1), 
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('birthyear.beta', 'eps.capyear','bodymass')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'weight.cohort.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(bodymass[birthyear]) #this creates a dataframe in long format with indexing
gather$birthyear <- as.factor(gather$birthyear)


plot2<- gather %>% 
  ggplot(aes(x=birthyear, y=.value)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels=c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015'))+
  # scale_fill_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  # scale_color_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  # 
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "YEAR", y = "WEIGHT (LB)", title = "WEIGHT BY YEAR (6 Y.O.)")+
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
plot2
ggsave('./figures/weight.cohort.jpg', plot2, width = 12, height = 6)

################################33
#Birthsite effect on cohorts
set.seed(100)
sink('ant.cohort.jags')
cat('
model{
  
#priors

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

# for (u in 1:12) { #site and age interaction
#   age.bs.beta[u] ~ dnorm(0, 0.001)
# }

for (u in 1:9) { #birth year effect
  birthyear.beta[u] ~ dnorm(0,0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:9) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }


# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- birthyear.beta[birthyear[i]] + bs.beta[bs[i]] + eps.capyear[capyear[i]] 
 
}

#derived parameter
  for (i in 1:9){ #birth year
    for (j in 1:3){ #birth site
      antlers[i,j] <- birthyear.beta[i] + bs.beta[j]
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
jags.data <- list(n=n,  antler = antlerin, bs = bs, birthyear = birthyear, capyear = capyear) #bs = bs,capyear = capyear, 

#inits function
inits<- function(){list(birthyear.beta = rnorm(9, 0,1), eps.capyear = rnorm(9,0,1), 
                        bs.beta = rnorm(3,0,1), sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('birthyear.beta', 'bs.beta', 'eps.capyear','antlers')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'ant.cohort.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(antlers[birthyear, birthsite]) #this creates a dataframe in long format with indexing
gather$birthyear <- as.factor(gather$birthyear)
gather$birthsite <- as.factor(gather$birthsite)


plot2<- gather %>% 
  ggplot(aes(x=birthyear, y=.value, group = birthsite, color = birthsite)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels=c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015'))+
  scale_fill_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +

  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "YEAR", y = "ANTLER SCORE (IN)", title = "ANTLERS BY YEAR (6 Y.O.)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/ant.site.cohort.jpg', plot2, width = 12, height = 6)


################################33
#interaction between site and cohort
set.seed(100)
sink('interact.cohort.jags')
cat('
model{
  
#priors

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:9) { #site and cohort interaction
  bs.by.beta[u] ~ dnorm(0, 0.001)
}

for (u in 1:9) { #birth year effect
  birthyear.beta[u] ~ dnorm(0,0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:9) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }


# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- birthyear.beta[birthyear[i]] + bs.beta[bs[i]] + bs.by.beta[bs[i]]*birthyear[i]+
                    eps.capyear[capyear[i]] 
 
}

#derived parameter
  for (i in 1:9){ #birth year
    for (j in 1:2){ #birth site
      antlers[i,j] <- birthyear.beta[i] + bs.beta[j] + bs.by.beta[j]*i
    }
}


}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n,  antler = antlerin, bs = bs, birthyear = birthyear, capyear = capyear) #bs = bs,capyear = capyear, 

#inits function
inits<- function(){list(birthyear.beta = rnorm(9, 0,1), eps.capyear = rnorm(9,0,1), 
                        bs.beta = rnorm(3,0,1), bs.by.beta = rnorm(9, 0, 1), 
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('birthyear.beta', 'bs.beta','bs.by.beta', 'eps.capyear','antlers')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'interact.cohort.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(antlers[birthyear, birthsite]) #this creates a dataframe in long format with indexing
gather$birthyear <- as.factor(gather$birthyear)
gather$birthsite <- as.factor(gather$birthsite)


plot2<- gather %>% 
  ggplot(aes(x=birthyear, y=.value, group = birthsite, color = birthsite)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels=c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015'))+
  scale_fill_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL")) +
  scale_color_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL")) +
  
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "YEAR", y = "ANTLER SCORE (IN)", title = "ANTLERS BY YEAR (6 Y.O.)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/ant.interaction.cohort.jpg', plot2, width = 12, height = 6)

################################33
#interaction between site and cohort weight
set.seed(100)
sink('weight.interact.cohort.jags')
cat('
model{
  
#priors

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:9) { #site and cohort interaction
  bs.by.beta[u] ~ dnorm(0, 0.001)
}

for (u in 1:9) { #birth year effect
  birthyear.beta[u] ~ dnorm(0,0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:9) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }


# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- birthyear.beta[birthyear[i]] + bs.beta[bs[i]] + bs.by.beta[bs[i]]*birthyear[i]+
                    eps.capyear[capyear[i]] 
 
}

#derived parameter
  for (i in 1:9){ #birth year
    for (j in 1:2){ #birth site
      mass[i,j] <- birthyear.beta[i] + bs.beta[j] + bs.by.beta[j]*i
    }
}


}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n,  weight = weightlb, bs = bs, birthyear = birthyear, capyear = capyear) #bs = bs,capyear = capyear, 

#inits function
inits<- function(){list(birthyear.beta = rnorm(9, 0,1), eps.capyear = rnorm(9,0,1), 
                        bs.beta = rnorm(3,0,1), bs.by.beta = rnorm(9, 0, 1), 
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('birthyear.beta', 'bs.beta','bs.by.beta', 'eps.capyear','mass')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'weight.interact.cohort.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(mass[birthyear, birthsite]) #this creates a dataframe in long format with indexing
gather$birthyear <- as.factor(gather$birthyear)
gather$birthsite <- as.factor(gather$birthsite)


plot2<- gather %>% 
  ggplot(aes(x=birthyear, y=.value, group = birthsite, color = birthsite)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels=c('2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015'))+
  scale_fill_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL")) +
  scale_color_discrete(name = "BIRTH SITE", labels = c("DMP", "CONTROL")) +
  
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "YEAR", y = "WEIGHT (LB)", title = "WEIGHT BY YEAR (6 Y.O.)")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/WEIGHT.interaction.cohort.jpg', plot2, width = 12, height = 6)
