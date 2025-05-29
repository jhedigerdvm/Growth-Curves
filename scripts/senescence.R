#usethis::use_git() #connect rproject to git
#senescent life long curves

library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)

# data <- read.csv('./raw/bucks_nofawns.csv', header=T)

data<- read.csv('./clean/nofawns22rain.csv', header = T)
data1<- data
data1 <- data1 %>% filter(year_birth >= '2011')
data <- data1
# ey <- data %>% filter(birthsite == "ey")
# unique(ey$animal_id)
# data <- ey
head(data)
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
weightkg <- data$weight.kg
weightlb <- data$weight.lb

#create vector with antlers
antlercm<- data$antler.cm
antlerin <- data$antler.in

#number of observations 
n <- nrow(data)


######################################


#Write linear mixed effects model for antlers inches,
set.seed(100)
sink('lme.ant.jags')
cat('
model{

#priors
int ~ dunif(0,100)
ageclass.beta[1] <- 0
# 
# eps.capyear[1]<-0
# 
# # eps.birthyear[1]<-0
# eps.id[1]<-0

bs.beta[1] <- 0
age.site.beta[1] <- 0

for(u in 2:11){ # age class
  ageclass.beta[u] ~ dunif(0,100)
}

for(u in 2:3){ #birth site
  bs.beta[u] ~ dnorm(0, 0.001)
}

for (u in 2:11) {#age site interaction
  age.site.beta[u] ~ dnorm(0, 0.001)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,100)
# 
# for(u in 2:494){ #random effect for id
#   eps.id[u] ~ dnorm(0, tau.id)
# }
#   tau.id <- 1/(sigma.id * sigma.id)
#   sigma.id ~ dunif(0, 100)
#   
# 
# for(u in 2:15){                    #random effect for capyear
#     eps.capyear[u] ~ dnorm(0, tau.capyear)
#   }
#   tau.capyear <- 1/(sigma.capyear * sigma.capyear)
#   sigma.capyear ~ dunif(0, 100)
  

# Likelihood
for (i in 1:n){
 antler[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- int + ageclass.beta[ageclass[i]] 
                  + bs.beta[bs[i]] 
                  + age.site.beta[ageclass[i]]*bs[i] 
                  # + eps.id[id[i]] 
                  # + eps.capyear[capyear[i]]
         
}

#derived parameter
for (j in 1:11){ #age
  for (k in 1:3){ #site
  antlerin[j,k] <- int + ageclass.beta[j] + bs.beta[k] + age.site.beta[j]*k
  }
}
# 
# for (j in 7:12){ #age
#   for (k in 1:3){
#   antler_diff[j, k] <- antlerin[j,k]-antlerin[8,k]
# }
# }

}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, ageclass = ageclass, bs = bs, id = animal_id,
                  antler = antlerin, capyear = capyear)#, birthyear = birthyear,id = animal_id)#, bs = bs

#inits function
inits<- function(){list(ageclass.beta = c(NA, runif(10,0,100)), int = runif(0,100),
                        sigma = rlnorm(1), bs.beta = c(NA, rnorm(2,0, 0.001)),
                        age.site.beta = c(NA, rnorm(10, 0, 0.001))#,
                      #  eps.id = c(NA, rnorm(493, 0, 0.001), eps.capyear = c(NA, rnorm(14,0,1)))
                        #, eps.capyear = c(NA,rnorm(13,0,1)), sigma.capyear = runif(0,100),
                        # eps.birthyear = c(NA, rnorm(13,0,1)), sigma.birthyear = runif(0,100),
                        # eps.id = c(NA,rnorm(474,0,1)), sigma.id = runif(0,100)#,
                       # eps.bs = c(NA, rnorm(2,0,1)), sigma.bs = runif(0,100)
) } #log normal pulls just positive values

#parameters to estimate
parameters <- c('ageclass.beta', 'bs.beta', 'age.site.beta', 'antlerin', 'eps.capyear', 'antler_diff')

#MCMC settings
ni <- 3000
nb <- 1000
nt<- 10
nc <- 3

lme.ant.jags<- jagsUI(jags.data, inits, parameters, 'lme.ant.jags', n.thin = nt, n.chains = nc,
                  n.burnin = nb, n.iter = ni, parallel = TRUE)
print(lme.ant.jags)
MCMCtrace(lme.ant.jags)
# 
# summary_ants<- lme.ant.jags$summary
# write.csv(summary_ants, 'antlersum.csv', row.names = T)
# 

#get output into a tibble before you can manipulate in ggplot
gather<- lme.ant.jags %>% gather_draws(antlerin[age,site])

gather$age <- as.factor(gather$age)
gather$site <- as.factor(gather$site)



# 
# plot_base_lme <-
#   ggplot(data = gather, aes(x=age, y=.value, color = site, fill = site))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.15,0.9),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 16),
#         plot.title = element_text(face = 'bold', size = 32, hjust = 0.5 ),
#         axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 24),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#         legend.background = element_rect(fill='transparent'), #transparent legend bg
#         legend.box.background = element_rect(fill='transparent')) #transparent legend panel)
# 
# 
# lme.plot<- plot_base_lme +
#   # stat_halfeye()+
#   stat_pointinterval(alpha = .9, .width = c(0.5, 0.95)) +
#   scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11'),
#                    labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
#                               '9.5','10.5','11.5'))+
#   labs(x = "AGE CLASS", y = "ANTLER SCORE (IN)",
#        title = "ANTLER SCORE BY AGE CLASS")


plot<- gather  %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  # facet_wrap(vars(age), nrow = 1)+
  stat_pointinterval(alpha=0.9, .width=c(0.5,0.95), position = position_dodge(width = 0.5))+
  scale_fill_viridis_d(alpha = .2, labels = c('DMP', 'EY', 'WY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('DMP', 'EY', 'WY'))+ #color of line but no opacification
  labs(x = "Age", y = "Antler (in)", title = "Antlers by Age and Site")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        # legend.position = c(0.15,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 21, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)



ggsave('./figures/antler.in2.png', plot, bg='transparent', width = 12, height = 9)

summary <- data %>%
  count(birthsite, age)
unique(data$animal_id)

