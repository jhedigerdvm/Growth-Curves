#linear model looking at the effect of age and site interaction on antlers
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
# data <- read.csv('./raw/bucks_nofawns.csv', header=T)
data<- read.csv('./clean/nofawns22.csv', header=T)
head(data)

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


#Write linear mixed effects model for antlers in, 
set.seed(100)
sink('ant.bs.age.jags')
cat('
model{
  
#priors

for(u in 1:12){ #ageclass
  ageclass.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:12) { #site and age interaction
  age.bs.beta[u] ~ dnorm(0, 0.001)
}
# 
# for (u in 1:14) { #birth year effect
#   birthyear.beta[u] ~ dnorm(0,0.001)
# }

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:15) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }

for (u in 1:494){ #random effect for animal id
  eps.id[u] ~ dnorm(0,sigma)
}

for (u in 1:15) { #random effect for birth year
  eps.birthyear[u] ~ dnorm(0, sigma)
}

# Likelihood 
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- ageclass.beta[ageclass[i]] + bs.beta[bs[i]] + age.bs.beta[ageclass[i]] * bs[i] + 
                                                        eps.capyear[capyear[i]] + 
                                                        eps.id[animal_id[i]] +
                                                        eps.birthyear[birthyear[i]]
}

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:10) { #ageclass
      antlers[i,j] <- ageclass.beta[j] + bs.beta[i] + age.bs.beta[j]*i
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
jags.data <- list(n=n, ageclass = ageclass, bs = bs, antler = antlerin, capyear = capyear, 
                  animal_id=animal_id, birthyear = birthyear)

#inits function
inits<- function(){list(ageclass.beta = rnorm(12,0,1), bs.beta = rnorm(3, 0, 1), age.bs.beta = rnorm(12, 0, 1),
                        eps.birthyear = rnorm(15, 0,1),eps.capyear = rnorm(15, 0,1), eps.id = rnorm(494,0,1),
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('ageclass.beta','bs.beta','eps.birthyear', 'age.bs.beta', 'antlers', 'antlers_diff')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

lme.jags<- jagsUI(jags.data, inits, parameters, 'ant.bs.age.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(lme.jags)
MCMCtrace(lme.jags)
# write.csv(lme.jags$summary, 'antlers.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- lme.jags %>% gather_draws(antlers[site, age]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)



plot2<- gather %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#0072B2", "black", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "black", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_x_discrete(labels = c("1.5", "2.5", "3.5", '4.5', '5.5', '6.5', '7.5', '8.5', '9.5', '10.5'))+
# scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "AGE", y = "ANTLER SCORE (IN)", title = "ANTLER SCORE BY AGE AND SITE")+
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
ggsave('./figures/antlers.age.site1.jpg', plot2, width = 10, height = 7)


plot2<- gather %>% filter(age==c("5","6","7","8", "9")) %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  # scale_x_discrete(labels = c("1.5", "2.5", "3.5", '4.5', '5.5', '6.5', '7.5', '8.5', '9.5', '10.5'))+
  
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "AGE", y = "ANTLER SCORE (IN)", title = "ANTLER SCORE BY AGE AND SITE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.5,0.15),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/antlers.age.site56789.jpg', plot2, width = 10, height = 7)


plot2<- gather %>% filter(age==c("5","6","7","8", "9")) %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.4)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  
  # scale_fill_viridis_d(begin = 0, end = 0.6, option = 'F', alpha = .2, 
  #                      labels = "DMP", "CONTROL", "TGT") + #this allowed me to opacify the ribbon but not the line
  # scale_color_viridis_d(option = 'F', begin = 0, end = 0.6,
  #                       labels = "DMP", "CONTROL", "TGT")+ #color of line but no opacification
  labs(x = "AGE", y = "ANTLER SCORE (IN)", title = "ANTLER SCORE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.7,0.12),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/antlers.TIGHT.jpg', plot2, width = 5, height = 7)

#Model predicting bodyweight in pounds as a function of site, age, and the interaction bt site and age
set.seed(100)
sink('weight.jags')
cat('
model{
  
#priors

for(u in 1:12){ #ageclass
  ageclass.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

for (u in 1:12) { #site and age interaction
  age.bs.beta[u] ~ dnorm(0, 0.001)
}
# 
# for (u in 1:15) { #birth year effect
#   birthyear.beta[u] ~ dnorm(0,0.001)
# }

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)


for (u in 1:15) {   #random effect for capture year
  eps.capyear[u] ~ dnorm(0, sigma)
  }

for (u in 1:494){ #random effect for animal id
  eps.id[u] ~ dnorm(0,sigma)
}

for (u in 1:15) { #random effect for birth year
  eps.birthyear[u] ~ dnorm(0, sigma)
}

# Likelihood 
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each weight is a draw from this distribution
 mu[i] <- ageclass.beta[ageclass[i]] + bs.beta[bs[i]] + age.bs.beta[ageclass[i]] * bs[i] + 
                                                        eps.capyear[capyear[i]] + 
                                                        eps.id[animal_id[i]] +
                                                        eps.birthyear[birthyear[i]]
}

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:10) { #ageclass
      bodymass[i,j] <- ageclass.beta[j] + bs.beta[i] + age.bs.beta[j]*i
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
jags.data <- list(n=n, ageclass = ageclass, bs = bs, weight = weightlb, capyear = capyear, 
                  animal_id=animal_id, birthyear = birthyear)

#inits function
inits<- function(){list(ageclass.beta = rnorm(12,0,1), bs.beta = rnorm(3, 0, 1), age.bs.beta = rnorm(12, 0, 1),
                        eps.birthyear = rnorm(15, 0,1),eps.capyear = rnorm(15, 0,1), eps.id = rnorm(494,0,1),
                        sigma = rlnorm(1))} #log normal pulls just positive values,

#parameters to estimate
parameters <- c('ageclass.beta','bs.beta', 'age.bs.beta', 'bodymass')#

#MCMC settings
ni <- 2000
nb <- 1000
nt<- 1
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc, 
                  n.burnin = nb, n.iter = ni)
print(weight.jags)
# MCMCtrace(lme.jags)
write.csv(weight.jags$summary, './output/weight.bs.age.intrxn.csv')


#create a tibble of the posterior draws
gather<- weight.jags %>% gather_draws(bodymass[site, age]) #this creates a dataframe in long format with indexing
gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)


plot3<- gather %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#0072B2", "black", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "black", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_x_discrete(labels = c("1.5", "2.5", "3.5", '4.5', '5.5', '6.5', '7.5', '8.5', '9.5', '10.5'))+
  labs(x = "AGE", y = "WEIGHT (LBS)", title = "WEIGHT BY AGE AND SITE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.6,0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot3
ggsave('./figures/weight.age.site1.jpg', plot3, width = 10, height = 7)



plot2<- gather %>% filter(age==c("5","6","7","8", "9")) %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  labs(x = "AGE", y = "WEIGHT (LBS)", title = "WEIGHT BY AGE AND SITE")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.45,0.12),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot2
ggsave('./figures/weight.age.site56789.jpg', plot2, width = 10, height = 8)


#squeeze 5-9
plot2<- gather %>% filter(age==c("5","6","7","8", "9")) %>% 
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(alpha = .5, .width = c(0.5, 0.95), 
                     position = position_dodge(width = 0.3)) +
  scale_fill_manual(values = c("#0072B2", "grey", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  scale_color_manual(values = c("#0072B2", "grey", "#CC79A7"),name = "BIRTH SITE", labels = c("DMP", "CONTROL", "TGT")) +
  labs(x = "AGE", y = "WEIGHT (LBS)", title = "WEIGHT")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('./figures/weight.tight.jpg', plot2, width = 5, height = 7)



#how many of each ageclass by site
count<- data %>% group_by(birthsite, age) %>% summarise(count=n())
