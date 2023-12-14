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

# 
# #quick summary stats of data
# dmp<-data %>% filter(data$birthsite=='dmp') 
# ey <- data %>% filter(data$birthsite=='ey') 
# wy <- data %>% filter(data$birthsite=='wy') 

# #assign group two sires
# sire.data<- sire.data
# sire.data$group <- recode(sire.data$group, 'dmp' = 'sire')
# sire.data$animal_id

#Write linear mixed effects model for antlers in, this model continues increasing linearly into the older age classes
set.seed(100)
sink('lme.jags')
cat('
model{
  
#priors
for(u in 1:12){
  ageclass.beta[u] ~ dunif(0,1000)
}
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- ageclass.beta[ageclass[i]] 
}
 
#derived parameter
for (j in 1:12){ #age
  antler_difference[j] <- ageclass.beta[6]-ageclass.beta[j]
}

}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass,
                  antler = antlerin)

#inits function
inits<- function(){list(ageclass.beta = runif(0,1000),
                        sigma = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'antler_difference')

#MCMC settings
ni <- 2000
nb <- 500
nt<- 10
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'lme.jags', n.thin = nt, n.chains = nc, 
                             n.burnin = nb, n.iter = ni)
print(lme.jags)
write.csv(print(lme.jags$summary), 'antleroutput.csv')
MCMCtrace(lme.jags)
MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results


#manually calculate F value to see what percent of the posterior overlaps zero
#first select column of interest, then calculate the proportion greater than zero
#overlap of zero, what is it? 
posterior.lme<- tidy_draws(lme.jags)
diff6<-posterior.lme$`antler_difference[6]`
six<-length(which(diff6>0))/length(diff6)
diff8<-posterior.lme$`antler_difference[8]`
eight<-length(which(diff8>0))/length(diff8)
diff9<-posterior.lme$`antler_difference[9]`
nine<-length(which(diff9>0))/length(diff9)
diff12<-posterior.lme$`antler_difference[12]`
twelve<-length(which(diff12>0))/length(diff12)


#at what age do most individuals max
max.antler.ind<- data3 %>%  group_by(animal_id) %>%  slice(which.max(bcscm)) #slice only selects rows containing max bcscm for individuals 
max.antler.age<-max.antler.ind[niners,] #create df with the individuals that lived nine years and the age they maxed antlers 
max_antlers_by_age <- hist(max.antler.ind$age)
table(max.antler.ind$age, useNA='always')  

#get output into a tibble before you can manipulate in ggplot
posterior.lme<- tidy_draws(lme.jags)
posterior.lme <- posterior.lme[,c(4:15)]

#pivot longer puts them in a tibble format
lme.long <- posterior.lme %>% pivot_longer(everything())

#create age class column and assign 1-12 to appropriate ageclass
lme.long$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(lme.long)/12))

#remove 'names' column 
lme.long <- lme.long[,-c(1)]

#rename values column to antlers in cm
lme.long <- lme.long %>% rename_at('value', ~'antlers')

#ggplot
plot_base_lme <- 
  ggplot(data = lme.long, aes(x=ageclass, y=antlers))+
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


################################################

#Write linear mixed effects model for weight pounds,
set.seed(100)
sink('lme.jags')
cat('
model{
  
#priors
for(u in 1:12){
  ageclass.beta[u] ~ dunif(0,1000)
}
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
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
                  weight = weightlb)

#inits function
inits<- function(){list(ageclass.beta = runif(0,1000),
                        sigma = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'weight_difference')

#MCMC settings
ni <- 4000
nb <- 1000
nt<- 10
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'lme.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
write.csv(print(lme.jags$summary), 'weightoutput.csv')
MCMCtrace(lme.jags)
MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results


#get output into a tibble before you can manipulate in ggplot
posterior.lme<- tidy_draws(lme.jags)
posterior.lme <- posterior.lme[,c(4:15)]

#pivot longer puts them in a tibble format
lme.long <- posterior.lme %>% pivot_longer(everything())

#create age class column and assign 1-12 to appropriate ageclass
lme.long$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(lme.long)/12))

#remove 'names' column 
lme.long <- lme.long[,-c(1)]

#rename values column to antlers in cm
lme.long <- lme.long %>% rename_at('value', ~'weight')


plot_base_lme <- 
  ggplot(data = lme.long, aes(x=ageclass, y=weight))+
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
  labs(x = "AGE CLASS", y = "WEIGHT (LB)", 
       title = "ANTLER SIZE BY AGE CLASS") 

ggsave('weight_lb.png', lme.plot, bg='transparent', width = 10, height = 10)




####################################
#LME antlers, birth site, random effects
####################################
set.seed(100)
sink('antlers.jags')
cat('
model{
  
#priors
  
  ageclass.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate
  bs.beta[1] <- 0     #need to set baseline of zero for >1 categorical covariate
  
  int~dunif(0,1000)        #need an intercept to build priors upon

  for (u in 2:3){
    bs.beta[u] ~ dnorm(0,0.001) #needed dnorm for bs because of negative values! 
  }
  
  for(u in 2:12){
    ageclass.beta[u] ~ dunif(0,1000)
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

  for(birthyear in 1:14){                    #random effect for birthyear
      eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
    }
  tau.birthyear <- 1/(sigma.birthyear * sigma.birthyear)
  sigma.birthyear ~ dunif(0, 100)
  

  for(animal_id in 1:475){                       #random effect for animal_id
    eps.id[animal_id] ~ dnorm(0, tau.id)
  }

  tau.id <- 1/(sigma.id * sigma.id)
  sigma.id ~ dunif(0, 100)

# Likelihood 
for (i in 1:n){
   antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
   mu[i] <- int + bs.beta[bs[i]] + (ageclass.beta[ageclass[i]]) + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] + eps.id[animal_id[i]]
  }

#derived parameter
for (i in 1:3){ #birthsite
  for (j in 1:12){ #age
    derived_antlers[i,j]<- int+ bs.beta[i] + ageclass.beta[j]
}
}

for(i in 1:2){
   site.diff[i] <- bs.beta[3] - bs.beta[i]
  }
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, bs = bs, antler = antlercm, 
                  capyear = capyear, birthyear = birthyear, animal_id = animal_id)

#inits function
inits<- function(){list(int = runif(1,0,1000), bs.beta = c(NA,rnorm(2,0,1)),
                        ageclass.beta = c(NA,runif(11,0,1000)), sigma = rlnorm(1),
                        eps.capyear = rnorm(14, 0, 1),sigma.capyear = rlnorm(1), 
                        eps.birthyear = rnorm(14,0,1),sigma.birthyear = rlnorm(1),
                        eps.id = rnorm(475, 0, 1), sigma.id = rlnorm(1)
                        
                        )} #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'bs.beta', 'derived_antlers', 'site.diff')

#MCMC settings
ni <- 15000
nb <- 10000
nt<- 10
nc <- 3

antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
MCMCtrace(antlers.jags)
print(antlers.jags)
write.csv(print(antlers.jags$summary), 'antlers_site_re.csv')
MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results

#get output into a tibble before you can manipulate in ggplot
posterior.antlers<- tidy_draws(antlers.jags)
posterior.antlers <- posterior.antlers[,-c(1:18)]
posterior.antlers <- posterior.antlers[,-c(34:39)]

#pivot longer puts them in a tibble format
antlers.long <- posterior.antlers %>% pivot_longer(everything())

#create age class column and assign 1-12 to appropriate ageclass
antlers.long$ageclass <- as.factor(rep(c('1','1','1','2','2','2','3','3','3','4','4','4',
                                     '5', '5', '5', '6','6', '6', '7', '7', '7', '8', '8','8',
                                     '9', '9', '9', '10', '10', '10', '11', '11', '11'), nrow(antlers.long)/33))

antlers.long$bs <- as.factor(rep(c('dmp','ey','wy'), nrow(antlers.long)/3))

#remove 'names' column 
antlers.long <- antlers.long[,-c(1)]

#rename values column to antlers in cm
antlers.long <- antlers.long %>% rename_at('value', ~'antlers')

plot_base_lme <- 
  ggplot(data = antlers.long, aes(x=ageclass, y=antlers, group = bs))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.15,0.90),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.title = element_text(face = 'bold', size = 40, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

antlers.plot<- plot_base_lme +
  stat_pointinterval(aes(color = bs), alpha = 1, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                     values=c("royalblue2", "green", "darkorchid"))+
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'), 
                   labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                              '9.5','10.5','11.5'))+
  labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
       title = "ANTLER SIZE BY AGE CLASS") 

ggsave('antlers_site_re10x10.png', antlers.plot, bg='transparent', width = 12, height = 10)






####################################
#LME bodyweight, birth site, random effects
####################################
set.seed(100)
sink('weight.jags')
cat('
model{
  
#priors
  
  ageclass.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate
  bs.beta[1] <- 0     #need to set baseline of zero for >1 categorical covariate
  
  int~dunif(30,50)        #need an intercept to build priors upon

  for (u in 2:3){
    bs.beta[u] ~ dnorm(0,0.001) #needed dnorm for bs because of negative values!
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

  for(birthyear in 1:14){                    #random effect for birthyear
      eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
    }
  tau.birthyear <- 1/(sigma.birthyear * sigma.birthyear)
  sigma.birthyear ~ dunif(0, 150)


  for(animal_id in 1:475){                       #random effect for animal_id
    eps.id[animal_id] ~ dnorm(0, tau.id)
  }

  tau.id <- 1/(sigma.id * sigma.id)
  sigma.id ~ dunif(0, 100)

# Likelihood 
for (i in 1:n){
   weight[i] ~ dnorm(mu[i], tau) #each weight is a draw from this distribution
   mu[i] <- int + (ageclass.beta[ageclass[i]])  + bs.beta[bs[i]] + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] + eps.id[animal_id[i]]
  }

#derived parameter
for (i in 1:3){ #birthsite
  for (j in 1:12){ #age
    derived_weight[i,j]<- int+ bs.beta[i] + ageclass.beta[j]
}
}

for(i in 1:2){
   site.diff[i] <- bs.beta[3] - bs.beta[i]
  }
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n,  weight = weightkg, bs = bs, ageclass = ageclass, 
                  capyear=capyear, birthyear = birthyear, animal_id = animal_id)

#inits function
inits<- function(){list(int = runif(1,30,50), sigma = rlnorm(1), bs.beta = c(NA,rnorm(2,0,1)),
                        eps.capyear = rnorm(14, 0, 1),sigma.capyear = rlnorm(1), 
                        eps.birthyear = rnorm(14,0,1),sigma.birthyear = rlnorm(1),
                        ageclass.beta = c(NA,runif(11,0,150)),
                        eps.id = runif(475, 0, 100), sigma.id = rlnorm(1)
                        
)} #log normal pulls just positive values

#parameters to estimate
parameters <- c('int','bs.beta','ageclass.beta', 
                  'derived_weight', 'site.diff')

#MCMC settings
ni <- 30000
nb <- 15000
nt<- 10
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
MCMCtrace(weight.jags)
print(weight.jags)
write.csv(print(weight.jags$summary), 'weight_site_re.csv')
MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results

#get output into a tibble before you can manipulate in ggplot
posterior.weight<- tidy_draws(weight.jags)
posterior.weight <- posterior.weight[,-c(1:19)]
posterior.weight <- posterior.weight[,-c(34:39)]

#pivot longer puts them in a tibble format
weight.long <- posterior.weight %>% pivot_longer(everything())

#create age class column and assign 1-12 to appropriate ageclass
weight.long$ageclass <- as.factor(rep(c('1','1','1','2','2','2','3','3','3','4','4','4',
                                     '5', '5', '5', '6','6', '6', '7', '7', '7', '8', '8','8',
                                     '9', '9', '9', '10', '10', '10', '11', '11', '11'), nrow(weight.long)/33))

weight.long$bs <- as.factor(rep(c('dmp','ey','wy'), nrow(weight.long)/3))

#remove 'names' column 
weight.long <- weight.long[,-c(1)]

#rename values column to antlers in cm
weight.long <- weight.long %>% rename_at('value', ~'weight')
plot_base_lme <- 
  ggplot(data = weight.long, aes(x=ageclass, y=weight, group = bs))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.15,0.90),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.title = element_text(face = 'bold', size = 40, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

weight.plot<- plot_base_lme +
  stat_pointinterval(aes(color = bs), alpha = 1, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                     values=c("royalblue2", "green", "darkorchid"))+
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'), 
                   labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                              '9.5','10.5','11.5'))+
  labs(x = "AGE CLASS", y = "BODY WEIGHT (KG)", 
       title = "WEIGHT BY AGE CLASS") 

ggsave('weight.png', weight.plot, bg='transparent', width = 12, height = 10)


####################################### 
#Write generalized linear mixed effects model for antlers cm with poisson distr
#######################################
sink('glm.jags')
cat('
model{
  
#priors
 for (u in 1:3){                             
    site.beta[u] ~ dnorm(0, 0.01)               #priors for birthsite
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

age.beta ~dnorm(0,0.001) 
 
# Likelihood 
for (i in 1:n){
 antler[i] ~ dpois(lambda[i])
 lambda[i] <- exp(site.beta[bs[i]] + (age.beta * ageclass[i])) 
 }


#derived parameter, creating a 3x12 matrix
for(i in 1:3){ #site
for (j in 1:12){ #age
  derived_antler[i,j] <- exp(site.beta[i] + (age.beta * j))
}
}
}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, bs = bs,
                  antler = antlercm)

#inits function
inits<- function(){list(site.beta = rnorm(3,0,1),
                        age.beta = rnorm(1,0,1),
                        sigma = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('site.beta', 'age.beta', 'derived_antler')

#MCMC settings
ni <- 2000
nb <- 500
nt<- 10
nc <- 3

glm.jags<- jagsUI(jags.data, inits, parameters, 'glm.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
#write.csv(out2$summary, file = 'antleroutput.csv')
print(glm.jags)
MCMCtrace(antlercm.jags)
MCMCplot(antlercm.jags, horiz = FALSE, params = 'derived_antler' ) #quick plot to look at results


#write quadratic model for antlers in cm

sink('quadratic_site.jags')
cat('
model{
  
#priors

#age.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate
#site.beta[1] <- 0     #need to set baseline of zero for >1 categorical covariate

 for (u in 1:3){                             
    site.beta[u] ~ dnorm(0, 0.01)               #priors for birthsite
}

#int~dnorm(0,0.001)        #need an intercept to build priors upon
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

age.beta ~dnorm(0,0.001) # need agebeta and agebeta2 for quadratic effect
age.beta2~dnorm(0,0.001)


for(capyear in 1:14){                    #random effect for capyear 
  eps.capyear[capyear] ~ dnorm(0, tau.capyear)
}
tau.capyear <- 1/(sigma.capyear * sigma.capyear)
sigma.capyear ~ dunif(0, 100) #random effect SD youre saying that each year-specific error term is coming
                              # from the same distributionm i.e. same mean (0) and standard
                              # deviation (sigma.period). The standard deviation is what connects all the
                              # year-specific random effects to the same distribution.
                              
for(birthyear in 1:14){                    #random effect for birthyear
    eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
  }
  tau.birthyear <- 1/(sigma.birthyear * sigma.birthyear)
  sigma.birthyear ~ dunif(0, 100)

for(animal_id in 1:475){                       #random effect for animal_id
  eps.id[animal_id] ~ dnorm(0, tau.id)
}

tau.id <- 1/(sigma.id * sigma.id)
sigma.id ~ dunif(0, 100)
 
# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau)
 mu[i] <- site.beta[bs[i]] + (age.beta * ageclass[i]) + (age.beta2 * ageclass[i]*ageclass[i]) + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] + eps.id[animal_id[i]]
}


#derived parameter, creating a 3x12 matrix
for(i in 1:3){ #site
for (j in 1:12){ #age
  derived_antler[i,j] <- site.beta[i] + age.beta*j + age.beta2*j*j #im still confused why we multiple by j and not [j]
}
}
}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, bs = bs, capyear=capyear, birthyear=birthyear, animal_id = animal_id,
                antler = antlercm)

#inits function
inits<- function(){list(eps.id = rnorm(475, 0, 1) , eps.capyear = rnorm(14, 0, 1), eps.birthyear = rnorm(14,0,1),
                      site.beta = rnorm(3,0,1),
                      age.beta = rnorm(1,0,1), age.beta2 = rnorm(1,0,1), 
                      sigma = rlnorm(1),sigma.capyear = rlnorm(1), sigma.id = rlnorm(1), 
                      sigma.birthyear = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('site.beta', 'age.beta', 'age.beta2', 'derived_antler')

#MCMC settings
ni <- 5000
nb <- 500
nt<- 10
nc <- 3

quadratic_site.jags<- jagsUI(jags.data, inits, parameters, 'quadratic_site.jags', n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni)
#write.csv(out2$summary, file = 'antleroutput.csv')
print(quadratic_site.jags)
MCMCtrace(antlercm.jags)
MCMCplot(antlercm.jags, horiz = FALSE, params = 'derived_antler' ) #quick plot to look at results

#get output into a tibble before you can manipulate in ggplot
posterior<- tidy_draws(out2)
antlers <- posterior[,-c(1:8,45)]
#pivot longer puts them in a tibble format
antlers_long <- antlers %>% pivot_longer(everything())
#create site column
antlers_long$site<- as.factor(rep(c(1:3), nrow(antlers_long)/3))
#create age class column
antlers_long$ageclass <- as.factor(rep(c('1','1','1','2','2','2','3','3','3','4','4','4',
                             '5', '5', '5', '6','6', '6', '7', '7', '7', '8', '8','8',
                             '9', '9', '9', '10', '10', '10', '11', '11', '11', 
                             '12', '12', '12'), nrow(antlers_long)/36))

#remove 'names' column that contained derived_antlers
antlers_long <- antlers_long[,-c(1)]

#rename values column to antlers in cm
antlers_long <- antlers_long %>% rename_at('value', ~'antlers')

plot_base <- 
ggplot(data = antlers_long, mapping = aes(x = ageclass, y = antlers, color = site))+
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

A<- plot_base +
stat_pointinterval(aes(color = site), alpha = 0.7, .width = c(0.5, 0.95), 
                   position = 'dodge') +
scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                   values=c("royalblue2", "yellow3", "darkorchid"))+
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
A


#shortcut to piping is control  shift m









#write model for weight in kg, quadratic equation

sink('weightkg.jags')
cat('
model{
  
#priors

 for (u in 1:3){                             
    site.beta[u] ~ dnorm(0, 0.01)               #priors for birthsite
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

age.beta ~dnorm(0,0.001)
age.beta2~dnorm(0,0.001)


for(capyear in 1:14){                    #random effect for capyear
  eps.capyear[capyear] ~ dnorm(0, tau.capyear)
}
tau.capyear <- 1/(sigma.capyear * sigma.capyear)
sigma.capyear ~ dunif(0, 100) #random effect SD youre saying that each year-specific error term is coming
                              # from the same distributionm i.e. same mean (0) and standard
                              # deviation (sigma.period). The standard deviation is what connects all the
                              # year-specific random effects to the same distribution.
                              
for(birthyear in 1:14){                    #random effect for birthyear
    eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
  }
  tau.birthyear <- 1/(sigma.birthyear * sigma.birthyear)
  sigma.birthyear ~ dunif(0, 100)

for(animal_id in 1:475){                       #random effect for animal_id
  eps.id[animal_id] ~ dnorm(0, tau.id)
}

tau.id <- 1/(sigma.id * sigma.id)
sigma.id ~ dunif(0, 100)
 
# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau)
 mu[i] <- site.beta[bs[i]] + (age.beta * ageclass[i]) + (age.beta2 * ageclass[i]*ageclass[i]) + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] + eps.id[animal_id[i]]
}


#derived parameter, creating a 3x12 matrix
for(i in 1:3){ #site
for (j in 1:12){ #age
  derived_weight[i,j] <- site.beta[i] + age.beta*j + age.beta2*j*j
}
}
}   
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, 
                ageclass = ageclass, 
                bs = bs, 
                capyear=capyear, 
                birthyear=birthyear, 
                animal_id = animal_id,
                weight = weightkg)

#inits function
inits<- function(){list(eps.id = rnorm(475, 0, 1) , 
                      eps.capyear = rnorm(14, 0, 1), 
                      eps.birthyear = rnorm(14,0,1),
                      site.beta = rnorm(3,0,1),
                      age.beta = rnorm(1,0,1), 
                      age.beta2 = rnorm(1,0,1), 
                      sigma = rlnorm(1),
                      sigma.capyear = rlnorm(1), 
                      sigma.id = rlnorm(1), 
                      sigma.birthyear = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('site.beta', 'age.beta', 'age.beta2', 'derived_weight')

#MCMC settings
ni <- 5000
nb <- 500
nt<- 10
nc <- 3

weightkg.jags<- jagsUI(jags.data, inits, parameters, 'weightkg.jags', n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni)
#write.csv(out2$summary, file = 'antleroutput.csv')
print(weightkg.jags)
MCMCtrace(weightkg.jags)
MCMCplot(weightkg.jags, horiz = FALSE, params = 'derived_weight' )

#get output into a tibble before you can manipulate in ggplot
posterior<- tidy_draws(out3)
weight <- posterior[,-c(1:8,45)]
#pivot longer puts them in a tibble format
weight_long <- weight %>% pivot_longer(everything())
head(weight_long)
#create site column
weight_long$site<- as.factor(rep(c(1:3), nrow(weight_long)/3))
#create age class column
weight_long$ageclass <- as.factor(rep(c('1','1','1','2','2','2','3','3','3','4','4','4',
                                       '5', '5', '5', '6','6', '6', '7', '7', '7', '8', '8','8',
                                       '9', '9', '9', '10', '10', '10', '11', '11', '11', 
                                       '12', '12', '12'), nrow(weight_long)/36))

head(weight_long)
#remove 'names' column that contained derived_antlers
weight_long <- weight_long[,-c(1)]

head(weight_long)
#rename values column to antlers in cm
weight_long <- weight_long %>% rename_at('value', ~'weight')

plot_base <- 
ggplot(data = weight_long, mapping = aes(x = ageclass, y = weight, color = site))+
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

weight.plot<- plot_base +
stat_pointinterval(aes(color = site), alpha = 0.7, .width = c(0.5, 0.95), 
                   position = 'dodge') +
scale_color_manual(name="BIRTH SITE", labels=c("TREATMENT", "CONTROL", "TGT"),
                   values=c("royalblue2", "yellow3", "darkorchid"))+
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "WEIGHT (KG)", 
     title = "WEIGHTS BY AGE CLASS") 
weight.plot









#write model for antler von bertalanffy growth curves, hyperparameter model

sink('vb.jags')
cat('
model{
  
#priors

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)
logLinf ~ dnorm(0,1/1000)#T(-10,10)
logK ~ dnorm(0,1/1000)#T(-5,5)
t0 ~ dnorm(0,1/1000)                  #theoretical age when antlers are zero cm


#Exponentiate parameters
Linf = exp(logLinf) #theoretical maximum mean antler size achieved, log scale to restrict to positive values
K = exp(logK)       #brody growth coefficient , log scale to restrict to positive values 

#Likelihood 
for (i in 1:n){
   antler[i] ~ dnorm(mu[i], tau)
   mu[i] <- Linf*(1-exp(-K*(ageclass[i]-t0))) 
}

#derived parameter, creating a 3x12 matrix
for (j in 1:12){ #age
  derived_antlers[j] <-  Linf*(1-exp(-K*(j-t0))) 
}
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, antler = antlercm)

#inits function
inits<- function(){list(sigma = rlnorm(1),
                      logLinf = rlnorm(1),
                        logK =  rlnorm(1),
                        t0 = rnorm(1,1,1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('Linf', 'K', 't0', 'sigma', 'derived_antlers')

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 10
nc <- 3

out5<- jagsUI(jags.data, inits, parameters, 'vb.jags', n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni)
output<-print(out5$summary)
write.csv(output, 'vb.output.csv')
MCMCtrace(out5)
MCMCplot(out5, horiz = FALSE, params = 'derived_antlers' )

#get output into a tibble before you can manipulate in ggplot
posterior<- tidy_draws(out5)
antlers <- posterior[,c(8:19)]
#pivot longer puts them in a tibble format
antlers_long <- antlers %>% pivot_longer(everything())
#create site column
#antlers_long$site<- as.factor(rep(c(1:3), nrow(antlers_long)/3))
#create age class column
antlers_long$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(antlers_long)/12))

#remove 'names' column that contained derived_antlers
antlers_long <- antlers_long[,-c(1)]

#rename values column to antlers in cm
antlers_long <- antlers_long %>% rename_at('value', ~'antlers')

#make new dataframe and assign group to raw data
antlers_raw<- data[,c(4,10)]
antlers_raw$group<-as.factor(rep(c('raw'),nrow(antlers_raw)))
antlers_raw <- antlers_raw %>% rename_at('age', ~'ageclass')
antlers_raw <- antlers_raw %>% rename_at('bcscm', ~'antlers')


#assign group to posterior data
antlers_long$group<-as.factor(rep(c('VB'),nrow(antlers_long)))

#merge raw data with posterior data
antlers_merge<- rbind(antlers_long,antlers_raw)

plot_base <- 
ggplot(data = antlers_merge, aes(x=ageclass, y=antlers, color = group))+
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

A<- plot_base +
stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
geom_smooth(se=FALSE) +
# stat_pointinterval(data = data, aes(x = age, y = bcscm),
#                    alpha = 0.7, .width = c(0.5, 0.95)) +
# geom_smooth(data=data, aes(x=age, y=bcscm), fill="red",
#             colour="red", linewidth=1, se=FALSE) +
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
A
ggsave('vb.png', A, bg='transparent', width = 10, height = 10)


#now lets work to graph the raw data, first calculate summary stats
data_no_na<- data %>% drop_na(bcscm)
bcscm_sum<- data_no_na %>% group_by(age) %>% summarise_at(vars(bcscm), list(mean=mean, sd=sd)) %>% as.data.frame()

plot_base <- 
ggplot(data = bcscm_sum, mapping = aes(x = age, y = mean))+
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


B<- plot_base +
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.3)+
geom_point()+
#geom_smooth()+
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
B
ggsave('mean_sd_raw.jpg', B, width = 10, height = 10)





##################
#Gompertz curve
##################

sink('gompertz.jags')
cat('
model{
  
#priors, original analysis by Doll and Jacquemin 2018 used informative priors that were specified by taking the mean and SD of each param

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)
logLinf ~ dnorm(0,1/1000)               #asymptote
b ~ dnorm(0,tau)                    #inflection point
t0 ~ dnorm(0,1/1000)                  #theoretical age when antlers are zero cm


#Exponentiate parameters
Linf = exp(logLinf) #theoretical maximum mean antler size achieved, log scale to restrict to positive values

#Likelihood 
for (i in 1:n){
  antler[i] ~ dnorm(mu[i], tau)
   mu[i] <-  Linf*exp(-exp(-b*ageclass[i]-t0)) 
}

#derived parameter
for (j in 1:12){
  derived_antlers[j] <- Linf*exp(-exp(-b*j-t0))
}
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, antler = antlercm)

#inits function
inits<- function(){list(sigma = rlnorm(1),
                      logLinf = rlnorm(1),
                      b=rlnorm(1),
                      t0 = rnorm(1,1,1))} #log normal pulls just positive values, why dont we include tau here? 

#parameters to estimate
parameters <- c('Linf', 'b', 't0', 'derived_antlers')

#MCMC settings
ni <- 3000
nb <- 1000
nt<- 10
nc <- 3

out6<- jagsUI(jags.data, inits, parameters, 'gompertz.jags', n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni)
print(out6)
MCMCtrace(out6)
MCMCplot(out6, horiz = FALSE, params = 'derived_antlers' )

#get output into a tibble before you can manipulate in ggplot
posterior1<- tidy_draws(out6)
antlers1 <- posterior1[,c(7:18)]

#pivot longer puts them in a tibble format
antlers_long1 <- antlers1 %>% pivot_longer(everything())

#create age class column
antlers_long1$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(antlers_long1)/12))

#remove 'names' column that contained derived_antlers
antlers_long1 <- antlers_long1[,-c(1)]

#rename values column to antlers in cm
antlers_long1 <- antlers_long1 %>% rename_at('value', ~'antlers')

#make new dataframe and assign group to raw data
# antlers_raw<- data[,c(4,10)]
# antlers_raw$group<-as.factor(rep(c('raw'),nrow(antlers_raw)))
# antlers_raw <- antlers_raw %>% rename_at('age', ~'ageclass')
# antlers_raw <- antlers_raw %>% rename_at('bcscm', ~'antlers')


#assign group to posterior data
antlers_long1$group<-as.factor(rep(c('gompertz'),nrow(antlers_long1)))

#merge raw data with posterior data
antlers_merge1<- rbind(antlers_long1,antlers_raw)

plot_base1 <- 
ggplot(data = antlers_merge1, aes(x=ageclass, y=antlers, color = group))+
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

gomp.plot<- plot_base1 +
stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
geom_smooth(se=FALSE) +
# stat_pointinterval(data = data, aes(x = age, y = bcscm),
#                    alpha = 0.7, .width = c(0.5, 0.95)) +
# geom_smooth(data=data, aes(x=age, y=bcscm), fill="red",
#             colour="red", linewidth=1, se=FALSE) +
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
gomp.plot
ggsave('gompertz.png', gomp.plot, bg='transparent', width = 10, height = 10)


############################
#logistic growth curve
###########################


#write model for antler logistic growth curves, hyperparameter model

sink('log.jags')
cat('
model{
  
#priors, original analysis by Doll and Jacquemin 2018 used informative priors that were specified by taking the mean and SD of each param

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)
logLinf ~ dnorm(0,1/1000)
logK ~ dnorm(0,1/1000)
t0 ~ dnorm(0,1/1000)                  #theoretical age when antlers are zero cm


#Exponentiate parameters
Linf = exp(logLinf) #theoretical maximum mean antler size achieved, log scale to restrict to positive values
K = exp(logK)       #brody growth coefficient , log scale to restrict to positive values 

#Likelihood 
for (i in 1:n){
   antler[i] ~ dnorm(mu[i], tau)
   mu[i] <- Linf/(1+exp(-K*(ageclass[i]-t0))) 
}

#derived parameter, creating a 3x12 matrix
for (j in 1:12){ #age
  derived_antlers[j] <-  Linf/(1+exp(-K*(j-t0))) 
}
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, antler = antlercm)

#inits function
inits<- function(){list(sigma = rlnorm(1),
                      logLinf = rlnorm(1),
                      logK =  rlnorm(1),
                      t0 = rnorm(1,1,1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('Linf', 'K', 't0', 'sigma', 'derived_antlers')

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 10
nc <- 3

out7<- jagsUI(jags.data, inits, parameters, 'log.jags', n.thin = nt, n.chains = nc, 
            n.burnin = nb, n.iter = ni)

print(out7)
MCMCtrace(out7)
MCMCplot(out7, horiz = FALSE, params = 'derived_antlers' )


#get output into a tibble before you can manipulate in ggplot
posterior2<- tidy_draws(out7)
antlers2 <- posterior2[,c(8:19)]

#pivot longer puts them in a tibble format
antlers_long2 <- antlers2 %>% pivot_longer(everything())

#create age class column
antlers_long2$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(antlers_long2)/12))

#remove 'names' column that contained derived_antlers
antlers_long2 <- antlers_long2[,-c(1)]

#rename values column to antlers in cm
antlers_long2 <- antlers_long2 %>% rename_at('value', ~'antlers')

#make new dataframe and assign group to raw data
# antlers_raw<- data[,c(4,10)]
# antlers_raw$group<-as.factor(rep(c('raw'),nrow(antlers_raw)))
# antlers_raw <- antlers_raw %>% rename_at('age', ~'ageclass')
# antlers_raw <- antlers_raw %>% rename_at('bcscm', ~'antlers')


#assign group to posterior data
antlers_long2$group<-as.factor(rep(c('logistic'),nrow(antlers_long2)))

#merge raw data with posterior data
antlers_merge2<- rbind(antlers_long2,antlers_raw)

plot_base2 <- 
ggplot(data = antlers_merge2, aes(x=ageclass, y=antlers, color = group))+
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

log.plot<- plot_base2 +
stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
geom_smooth(se=FALSE) +
# stat_pointinterval(data = data, aes(x = age, y = bcscm),
#                    alpha = 0.7, .width = c(0.5, 0.95)) +
# geom_smooth(data=data, aes(x=age, y=bcscm), fill="red",
#             colour="red", linewidth=1, se=FALSE) +
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
log.plot
ggsave('logistic.png', log.plot, bg='transparent', width = 10, height = 10)



###########################
# Quadratic hyperparameter
##########################

#write model for antlers in cm

sink('quadratic.jags')
cat('
model{
  
#priors

# age.beta[1] <- 0      #need to set baseline of zero for >1 categorical covariate
#int~dnorm(0,0.001)        #need an intercept to build priors upon
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

age.beta ~dnorm(0,0.001) # need agebeta and agebeta2 for quadratic effect
age.beta2~dnorm(0,0.001)

# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau)
 mu[i] <- (age.beta * ageclass[i]) + (age.beta2 * ageclass[i]*ageclass[i]) 
}


#derived parameter, creating a 3x12 matrix
for (j in 1:12){ #age
  derived_antler[j] <- age.beta*j + age.beta2*j*j 
}
}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, antler = antlercm)

#inits function
inits<- function(){list(age.beta = rnorm(1,0,1), age.beta2 = rnorm(1,0,1), 
                      sigma = rlnorm(1))} #log normal pulls just positive values

#parameters to estimate
parameters <- c('age.beta', 'age.beta2', 'derived_antler')

#MCMC settings
ni <- 3000
nb <- 1000
nt<- 10
nc <- 3

quadratic.jags<- jagsUI(jags.data, inits, parameters, 'quadratic.jags', n.thin = nt, n.chains = nc, 
                      n.burnin = nb, n.iter = ni)
#write.csv(out2$summary, file = 'antleroutput.csv')
print(quadratic.jags)
MCMCtrace(quadratic.jags)
MCMCplot(quadratic.jags, horiz = FALSE, params = 'derived_antler' ) #quick plot to look at results

#get output into a tibble before you can manipulate in ggplot
posterior3<- tidy_draws(quadratic.jags)
#first select column of interest, then calculate the proportion greater than zero
#overlap of zero, what is it? 
#at what age do individuals max
data$age<-as.integer(data$age)
data_no_na<- data %>% drop_na(bcscm)
bcscm_sum<- data_no_na %>% group_by(age) %>% summarise_at(vars(bcscm), list(mean=mean, sd=sd)) %>% as.data.frame()
data2<-data %>% group_by(animal_id) %>%  summarise(oldest=max(age))
max_ant<- which(data2$oldest>=9)
slice<- data %>%  group_by(animal_id) %>%  slice(which.max(bcscm))
data2<-slice[max_ant,]
hist(data2$age)

antlers3 <- posterior3[,c(6:17)]
#pivot longer puts them in a tibble format
antlers_long3 <- antlers3 %>% pivot_longer(everything())
#create age class column
antlers_long3$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(antlers_long3)/12))

#remove 'names' column that contained derived_antlers
antlers_long3 <- antlers_long3[,-c(1)]

#rename values column to antlers in cm
antlers_long3 <- antlers_long3 %>% rename_at('value', ~'antlers')

#assign group to posterior data
antlers_long3$group<-as.factor(rep(c('quadratic'),nrow(antlers_long3)))

#merge raw data with posterior data
antlers_merge3<- rbind(antlers_long3,antlers_raw)

plot_base3 <- 
ggplot(data = antlers_merge3, aes(x=ageclass, y=antlers, color = group))+
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

quad.plot<- plot_base3 +
stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
geom_smooth(se=FALSE) +
# stat_pointinterval(data = data, aes(x = age, y = bcscm),
#                    alpha = 0.7, .width = c(0.5, 0.95)) +
# geom_smooth(data=data, aes(x=age, y=bcscm), fill="red",
#             colour="red", linewidth=1, se=FALSE) +
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
quad.plot
ggsave('quad.png', quad.plot, bg='transparent', width = 10, height = 10)


###############################
#merge all growth curves and raw data into one plot
all.curves<-rbind(antlers_long,antlers_long1,antlers_long2,antlers_long3, antlers_raw)

plot_base3 <- 
ggplot(data = all.curves, aes(x=ageclass, y=antlers, color = group))+
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

all.curves.plot<- plot_base3 +
stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
geom_smooth(se=FALSE) +
# stat_pointinterval(data = data, aes(x = age, y = bcscm),
#                    alpha = 0.7, .width = c(0.5, 0.95)) +
# geom_smooth(data=data, aes(x=age, y=bcscm), fill="red",
#             colour="red", linewidth=1, se=FALSE) +
scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'), 
                 labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                            '9.5','10.5','11.5','12.5'))+
labs(x = "AGE CLASS", y = "ANTLER SIZE (CM)", 
     title = "ANTLER SIZE BY AGE CLASS") 
all.curves.plot
ggsave('allcurves.png', all.curves.plot, bg='transparent', width = 10, height = 10)

############################
#WAIC for models
# compute wAIC for quadratic
library(R2jags)
quadratic.waic <- jags.samples(quadratic.jags$model, 
                         c("WAIC","deviance"), 
                         type = "mean", 
                         n.iter = 5000,
                         n.burnin = 1000,
                         n.thin = 1)
quadratic.waic$p_waic <- quadratic.waic$WAIC
quadratic.waic$waic <- quadratic.waic$deviance + quadratic.waic$p_waic
tmp <- sapply(quadratic.waic, sum)
waic.mquad <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

# compute wAIC for logistic model
log.waic <- jags.samples(out7$model, 
                         c("WAIC","deviance"), 
                         type = "mean", 
                         n.iter = 5000,
                         n.burnin = 1000,
                         n.thin = 1)
log.waic$p_waic <- log.waic$WAIC
log.waic$waic <- log.waic$deviance + log.waic$p_waic
tmp <- sapply(log.waic, sum)
waic.mlog <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


# compute wAIC for gompertz model
gomp.waic <- jags.samples(out6$model, 
                       c("WAIC","deviance"), 
                       type = "mean", 
                       n.iter = 5000,
                       n.burnin = 1000,
                       n.thin = 1)
gomp.waic$p_waic <- gomp.waic$WAIC
gomp.waic$waic <- gomp.waic$deviance + gomp.waic$p_waic
tmp <- sapply(gomp.waic, sum)
waic.mgomp <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


# compute wAIC for vb model
vb.waic <- jags.samples(out5$model, 
                        c("WAIC","deviance"), 
                        type = "mean", 
                        n.iter = 5000,
                        n.burnin = 1000,
                        n.thin = 1)
vb.waic$p_waic <- vb.waic$WAIC
vb.waic$waic <- vb.waic$deviance + vb.waic$p_waic
tmp <- sapply(vb.waic, sum)
waic.mvb <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)

# The difference in wAIC tells us that the covariate has some effect
# wAIC of m1 the model with covariate << wAIC of m0 intercept only
waic.df<-data.frame(von_Bertalanffy = waic.mvb, Gompertz = waic.mgomp, Logistic = waic.mlog, Quadratic = waic.mquad)
write.csv(waic.df,'waic.csv')


##############
# Summarstats
##############

data$age<-as.integer(data$age)
data_no_na<- data %>% drop_na(bcscm)
bcscm_sum<- data_no_na %>% group_by(age) %>% summarise_at(vars(bcscm), list(mean=mean, sd=sd)) %>% as.data.frame()
data2<-data %>% group_by(animal_id) %>%  summarise(oldest=max(age))
max_ant<- which(data2$oldest>=9)
slice<- data %>%  group_by(animal_id) %>%  slice(which.max(bcscm))
data2<-slice[max_ant,]
hist(data2$age)
