#usethis::use_git() #connect rproject to git
#senescent life long curves

library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

data <- read.csv('./raw/bucks_nofawns.csv', header=T)


#replace zeros with NA
data[data==0] <- NA

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$birthsite))
unique(data$birthsite)
unique(bs) 
bs  #dmp is 1, e yana 2, w yana 3

#create vectors of things of interest
data$age <- data$year_cap - data$year_birth 
ageclass <- data$age
ageclass

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
# 
# #Write linear mixed effects model for weight pounds,
# set.seed(100)
# sink('lme.jags')
# cat('
# model{
# 
# #priors
# int ~ dunif(0,100)
# ageclass.beta[1] <- 0
# 
# for(u in 2:12){
#   ageclass.beta[u] ~ dunif(0,100)
# }
# tau <- 1/(sigma*sigma)
# sigma ~ dunif(0,100)
# 
# for(capyear in 1:14){                    #random effect for capyear
#     eps.capyear[capyear] ~ dnorm(0, tau.capyear)
#   }
#   tau.capyear <- 1/(sigma.capyear * sigma.capyear)
#   sigma.capyear ~ dunif(0, 100)
# 
# for(birthyear in 1:14){
#   eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
# }
#   tau.birthyear <- 1/(sigma.birthyear*sigma.birthyear)
#   sigma.birthyear ~ dunif(0,100)
# 
# for(id in 1:475){
#   eps.id[id] ~ dnorm(0, tau.id)
# }
#   tau.id <- 1/(sigma.id*sigma.id)
#   sigma.id ~ dunif(0,100)
# 
# # Likelihood
# for (i in 1:n){
#  weight[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
#  mu[i] <- int + ageclass.beta[ageclass[i]] + eps.capyear[capyear[i]] + eps.birthyear[birthyear[i]] +
#                     eps.id[id[i]]
# }
# 
# #derived parameter
# for (j in 1:12){
#   weightlb[j] <- int + ageclass.beta[j]
# }
# 
# for (j in 1:12){ #age
#   weight_difference[j] <- ageclass.beta[7]-ageclass.beta[j]
# }
# 
# }
# ',fill = TRUE)
# sink()
# 
# #bundle data
# jags.data <- list(n=n, ageclass = ageclass, birthyear = birthyear,
#                   weight = weightlb, capyear = capyear, id = animal_id)
# 
# #inits function
# inits<- function(){list(ageclass.beta = c(NA, runif(11,0,100)), int = runif(0,100),
#                         sigma = rlnorm(1), eps.capyear = rnorm(14,0,1), sigma.capyear = runif(0,100),
#                         eps.birthyear = rnorm(14,0,1), sigma.birthyear = runif(0,100),
#                         eps.id = rnorm(475,0,1), sigma.id = runif(0,100)
#                         ) } #log normal pulls just positive values
# 
# #parameters to estimate
# parameters <- c('ageclass.beta', 'weightlb', 'weight_difference')
# 
# #MCMC settings
# ni <- 7000
# nb <- 3000
# nt<- 10
# nc <- 3
# 
# lme.jags<- jagsUI(jags.data, inits, parameters, 'lme.jags', n.thin = nt, n.chains = nc,
#                   n.burnin = nb, n.iter = ni)
# print(lme.jags)
# MCMCtrace(lme.jags)
# # MCMCplot(lme.jags, horiz = FALSE, params = 'ageclass.beta' ) #quick plot to look at results
# # write.csv(print(lme.jags$summary), 'weightoutput.csv')
# 
# 
# #get output into a tibble before you can manipulate in ggplot
# posterior.lme<- tidy_draws(lme.jags)
# posterior.lme <- posterior.lme[,c(16:27)]
# 
# #pivot longer puts them in a tibble format
# lme.long <- posterior.lme %>% pivot_longer(everything())
# 
# #create age class column and assign 1-12 to appropriate ageclass
# lme.long$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11','12'), nrow(lme.long)/12))
# 
# #remove 'names' column
# lme.long <- lme.long[,-c(1)]
# 
# #rename values column to antlers in cm
# lme.long <- lme.long %>% rename_at('value', ~'weight')
# 
# 
# plot_base_lme <-
#   ggplot(data = lme.long, aes(x=ageclass, y=weight))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.15,0.9),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 16),
#         plot.title = element_text(face = 'bold', size = 24, hjust = 0.5 ),
#         axis.title = element_text(face = 'bold',size = 20, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 16))
# 
# lme.plot<- plot_base_lme +
#   stat_pointinterval(alpha = 0.7, .width = c(0.5, 0.95)) +
#   scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'),
#                    labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
#                               '9.5','10.5','11.5','12.5'))+
#   labs(x = "AGE CLASS", y = "WEIGHT (LB)",
#        title = "WEIGHT BY AGE CLASS")
# 
# ggsave('weight_lb.png', lme.plot, bg='transparent', width = 10, height = 5)

######################################


#Write linear mixed effects model for antlers inches,
set.seed(100)
sink('lme.ant.jags')
cat('
model{

#priors
int ~ dunif(0,100)
ageclass.beta[1] <- 0
eps.capyear[1]<-0
eps.birthyear[1]<-0
eps.id[1]<-0
eps.bs[1] <- 0

for(u in 2:12){
  ageclass.beta[u] ~ dunif(0,100)
}
tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)

for(capyear in 2:14){                    #random effect for capyear
    eps.capyear[capyear] ~ dnorm(0, tau.capyear)
  }
  tau.capyear <- 1/(sigma.capyear * sigma.capyear)
  sigma.capyear ~ dunif(0, 100)

for(birthyear in 2:14){
  eps.birthyear[birthyear] ~ dnorm(0, tau.birthyear)
}
  tau.birthyear <- 1/(sigma.birthyear*sigma.birthyear)
  sigma.birthyear ~ dunif(0,100)

for(id in 2:475){
  eps.id[id] ~ dnorm(0, tau.id)
}
  tau.id <- 1/(sigma.id*sigma.id)
  sigma.id ~ dunif(0,100)

for(bs in 2:3){
  eps.bs[bs] ~ dnorm(0, tau.bs)
}
  tau.bs <- 1/(sigma.bs*sigma.bs)
  sigma.bs ~ dunif(0,100)

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

summary_ants<- lme.ant.jags$summary
write.csv(summary_ants, 'antlersum.csv', row.names = T)


#get output into a tibble before you can manipulate in ggplot
posterior.lme<- tidy_draws(lme.ant.jags)
posterior.lme <- posterior.lme[,c(16:26)]

#pivot longer puts them in a tibble format
lme.long <- posterior.lme %>% pivot_longer(everything())

#create age class column and assign 1-12 to appropriate ageclass
lme.long$ageclass <- as.numeric(rep(c('1','2','3','4','5','6','7','8','9','10','11'), nrow(lme.long)/11))

#remove 'names' column
lme.long <- lme.long[,-c(1)]

#rename values column to antlers in cm
lme.long <- lme.long %>% rename_at('value', ~'antler')


plot_base_lme <-
  ggplot(data = lme.long, aes(x=ageclass, y=antler))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.15,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5 ),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 24),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) #transparent legend panel)


lme.plot<- plot_base_lme +
  # stat_halfeye()+
  stat_pointinterval(alpha = .9, .width = c(0.5, 0.95)) +
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'),
                   labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                              '9.5','10.5','11.5'))+
  labs(x = "AGE CLASS", y = "ANTLER SCORE (IN)",
       title = "ANTLER SCORE BY AGE CLASS")

ggsave('./figures/antler.in.png', lme.plot, bg='transparent', width = 12, height = 9)

count<- data %>% group_by(age) %>%  summarise(n=n())
unique(data$animal_id)

