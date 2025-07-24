#final analysis
#linear model looking at the effect of age and site interaction on antlers
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
# data.allyears<- read.csv('./clean/nofawns22rain.csv', header=T)
# head(data)
# 
# data.subset <- data %>% filter(year_birth >= '2011')
# hist(data.allyears$annual.by)
# hist(data.subset$annual.by)
# data<- data.allyears

data <- read.csv('./clean/nofawns22rain.bin.csv', header = T)

# data.allo <- data %>% filter(year_birth >= '2011', !is.na(antler.cm), !is.na(weight.kg)) #only use when calculating allometric scaling

#filter for data to include years where all three sites are on line
data <- data %>% filter(year_birth >= '2011')

#replace zeros with NA
data[data==0] <- NA

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$birthsite))
bs

#create vectors of things of interest
data$age <- data$year_cap - data$year_birth 
ageclass <- data$age
count <- data %>% count(birthsite,age)

ageclass.sc <- as.vector(scale(ageclass)) #scale and center ageclass to treat as continuous cov.

#create a vector for capture year
capyear <- data$year_cap-min(data$year_cap)+1

#create vector for birthyear, so I can fit a random effect later
birthyear<-data$year_birth-min(data$year_birth)+1

#create vector with animal_id for the 494 unique individuals
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
n

#rainfall data
data$annual.cy.sc <- scale(data$annual.cy)
annualcy <- as.vector(data$annual.cy.sc)  #total annual rainfall for year of capture
data$annual.by.sc <- scale(data$annual.by)
annualby<-as.vector(data$annual.by.sc)  #total annual rainfall for birthyear of individual

#create simulated rainfall data
nvalues <- 100
annual.cy.sim <- seq(from = min(annualcy, na.rm = T), to = max(annualcy, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
annual.cy.sim #capture year rainfall sim

annual.by.sim <-seq(from = min(annualby, na.rm = T), to = max(annualby, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data
annual.by.sim #birth year rainfall sim

#create vector of environmental matching
match <- as.numeric(data$match)

#create simulated age data
nvalues <- 100
ageclass.sim <- seq(from = min(ageclass.sc, na.rm = T), to = max(ageclass.sc, na.rm = T), length.out = nvalues) #obtained to and from values from max and min of annual rainfall in data

#create log of antlers and weight
data$antlers.log <- log(data$antler.cm)
data$weight.log <- log(data$weight.kg)

antlers.log <- data$antlers.log
weight.log <- data$weight.log


# #function for weight matrix
# weight.init <- weight.log
# weight.init[is.na(weight.init)]<-mean(weight.log, na.rm = T) #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
# weight.init[!is.na(weight.log)]<-NA
# 
# 
# #need to find where NA values are in weight matrix, will use this information to build priors
# indices <- which(is.na(weight.log), arr.ind=T)
# indices <- indices %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
# NA_indices_weight <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
# for(i in 1:nrow(indices)){
#   NA_indices_weight[indices[[i,1]],indices[[i,3]]] <- indices[[i,2]]
# }
# weight<-as.matrix(weight)
# 
# #how many occasions does each individual have of an NA weight
# occasions_weight <- rowSums(is.na(weight))
# 
# ##now do the same thing we did with weight with antlers
# antlers<- pivot_wider(data, names_from = 'year', values_from = 'bcsin', id_cols = 'animal_id' )
# antlers<- as.matrix(antlers[,-1])
# 
# antlers <- scale(antlers) # scale and center
# 
# 
# #function for weight matrix
# antlers.init <- antlers
# antlers.init[is.na(antlers.init)]<-0 #applying mean weight to initial values for NA observations, because its scaled and centered, we can just use zero? 
# antlers.init[!is.na(antlers)]<-NA
# for (i in 1:dim(antlers.init)){ #cant have mean weight for years before animal was first captured
#   antlers.init[i,1:f[i]]<-NA
# }
# 
# #need to find where NA values are in weight matrix, will use this information to build priors
# antlers <- as.data.frame(antlers)
# indices_antlers <- as.data.frame(which(is.na(antlers), arr.ind=T))
# indices_antlers <- indices_antlers %>% arrange(row) %>%  group_by(row) %>%  mutate(n=1:n()) %>% ungroup()
# NA_indices_antlers <- matrix(NA, nrow=nrow(ch), ncol=ncol(ch))
# for(i in 1:nrow(indices_antlers)){
#   NA_indices_antlers[indices_antlers[[i,1]],indices_antlers[[i,3]]] <- indices_antlers[[i,2]]
# }
# antlers<-as.matrix(antlers)
# 
# #how many occasions does each individual have of an NA weight
# occasions_antlers <- rowSums(is.na(antlers))


# ---- Model1: weight ~ site + age (continuous)+ site x age ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.jags')
cat('
model{

#priors

beta1 ~ dnorm(0, 0.001) #age as continuous
          

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
  beta3[u] ~ dnorm(0, 0.001)  # site x age interaction
}



# 
# for (u in 1:11){ #random effect for capture year
#   eps1[u] ~ dnorm(0,sigma)
# }

#
# for (u in 1:494){ #random effect for animal-id
#   eps2[u] ~ dnorm(0, sigma)
# }


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]  #ageclass continuous
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[bs[i]]*ageclass[i]  #interaction betw by age and site
              
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #ageclass
      response[j,i] <- beta1*ageclass.sim[j] + beta2[i] + beta3[i]*ageclass.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = weightlb,
                  ageclass = ageclass.sc,
                  ageclass.sim = ageclass.sim,
                  capyear=capyear,
                  animal_id = animal_id,
                  match=match
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #interaction
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3','response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(weight.jags)
MCMCtrace(weight.jags)
write.csv(weight.jags$summary, './output/weight.agexsite.cont.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$age == 2)[1] 

#unscale and uncenter weight
age.sim.usc <- (ageclass.sim * sd(data$age, na.rm = T)) + mean(data$age, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
age.sim.usc1 <- for (i in age.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$ageclass <- vector

plot<- gather %>%
  ggplot(aes(x=ageclass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Body Mass (lbs)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/weight.agexsite.cont.jpg', plot, width = 15, height = 10)


# ---- Model2: antlers ~ site + age (continuous)+ site x age ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('antlers.jags')
cat('
model{

#priors

beta1 ~ dnorm(0, 0.001) #age as continuous
          

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
  beta3[u] ~ dnorm(0, 0.001)  # site x age interaction
}



# 
# for (u in 1:11){ #random effect for capture year
#   eps1[u] ~ dnorm(0,sigma)
# }

#
# for (u in 1:494){ #random effect for animal-id
#   eps2[u] ~ dnorm(0, sigma)
# }


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]  #ageclass continuous
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[bs[i]]*ageclass[i]  #interaction betw by age and site
              
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #ageclass
      response[j,i] <- beta1*ageclass.sim[j] + beta2[i] + beta3[i]*ageclass.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = antlerin,
                  ageclass = ageclass.sc,
                  ageclass.sim = ageclass.sim,
                  capyear=capyear,
                  animal_id = animal_id,
                  match=match
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #interaction
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3','response')  #'beta3',

#MCMC settings
ni <- 5000
nb <- 2000
nt<- 10
nc <- 3

antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(antlers.jags)
MCMCtrace(antlers.jags)
write.csv(antlers.jags$summary, './output/antlers.agexsite.cont.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$age == 2)[1] 

#unscale and uncenter weight
age.sim.usc <- (ageclass.sim * sd(data$age, na.rm = T)) + mean(data$age, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
age.sim.usc1 <- for (i in age.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$ageclass <- vector

plot<- gather %>%
  ggplot(aes(x=ageclass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/antlers.agexsite.cont.jpg', plot, width = 15, height = 10)



# ---- Model3: allometric scaling  ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('allo.jags')
cat('
model{

#priors
beta1 ~ dnorm(0, 0.001) #log weight as continuous
beta4 ~ dnorm(0, 0.001)   # age continuous

for (u in 1:3){
  beta2[u] ~ dnorm(0, 0.001)  #site effect
  beta3[u] ~ dnorm(0, 0.001)  #site and log weight interaction
}

#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*morpho[i]   #log weight continuous
              + beta2[bs[i]]  # site effect
              + beta3[bs[i]]*morpho[i]
              + beta4*age[i]
          
 }





}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, y = antlers.log,
                  morpho = weight.log,
                  bs = bs,
                  age = ageclass.sc
                
                  
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #log weight
                        beta2 = rnorm(3,0,1), # site effect
                        beta3 = rnorm(3,0,1), #site x log weight
                        beta4 = rnorm(1,0,1)
                        )}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4')  #,

#MCMC settings
ni <- 5000
nb <- 2000
nt<- 10
nc <- 3

allo.jags<- jagsUI(jags.data, inits, parameters, 'allo.jags', n.thin = nt, n.chains = nc,
                      n.burnin = nb, n.iter = ni, parallel = T)
print(allo.jags)
MCMCtrace(allo.jags)
write.csv(allo.jags$summary, './output/allometrytest.csv')




# ---- Model4: weight ~ site + age (categorical)+ age x site ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.jags')
cat('
model{

#priors

for(u in 1:11){
  beta1[u] ~ dnorm(0, 0.001) #age as categorical
  beta3[u] ~ dnorm(0, 0.001)  # age x site interaction

}  

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
}



# 
# for (u in 1:11){ #random effect for capture year
#   eps1[u] ~ dnorm(0,sigma)
# }

#
# for (u in 1:494){ #random effect for animal-id
#   eps2[u] ~ dnorm(0, sigma)
# }


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1[ageclass[i]]  #ageclass categorical
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[ageclass[i]]*bs[i]  #interaction betw by age and site
              
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:11){ #ageclass
      response[j,i] <- beta1[j] + beta2[i] + beta3[j]*i

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = weightlb,
                  ageclass = ageclass,
                  capyear=capyear,
                  animal_id = animal_id
)

#inits function
inits<- function(){list(beta1 = rnorm(11,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(11,0,1), #interaction
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3','response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(weight.jags)
MCMCtrace(weight.jags)
write.csv(weight.jags$summary, './output/weight.agexsite.cat.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)


plot<- gather %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = 'dodge')+
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Body Mass (lbs)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/weight.agexsite.cont.jpg', plot, width = 15, height = 10)



# ---- Model5: weight ~ site + age (categorical)+ site x age ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.jags')
cat('
model{

#priors

for(u in 1:11){
  beta1[u] ~ dnorm(0, 0.001) #age as categorical

}  

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
  beta3[u] ~ dnorm(0, 0.001)  # site x age interaction
}

beta4 ~ dnorm(0, 0.001) #rainfall



for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1[ageclass[i]]  #ageclass categorical
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[bs[i]]*ageclass[i]  #interaction betw by age and site
              + beta4*rain[i]
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:11){ #ageclass
      response[j,i] <- beta1[j] + beta2[i] + beta3[i]*j

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = weightlb,
                  ageclass = ageclass,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualby
)

#inits function
inits<- function(){list(beta1 = rnorm(11,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #interaction
                        beta4 = rnorm(1,0,1), # birthyear rainfall
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(weight.jags)
MCMCtrace(weight.jags)
write.csv(weight.jags$summary, './output/weight.sitexage.cat.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)


plot<- gather %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = 'dodge')+
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Body Mass (lbs)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/weight.sitexage.cont.jpg', plot, width = 15, height = 10)



# ---- Model6: antlers ~ site + age (categorical)+ site x age ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('antlers.jags')
cat('
model{

#priors

for(u in 1:11){
  beta1[u] ~ dnorm(0, 0.001) #age as categorical

}  

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
  beta3[u] ~ dnorm(0, 0.001)  # site x age interaction
}

beta4 ~ dnorm(0, 0.001) #rainfall



for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1[ageclass[i]]  #ageclass categorical
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[bs[i]]*ageclass[i]  #interaction betw by age and site
              + beta4*rain[i]
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:11){ #ageclass
      response[j,i] <- beta1[j] + beta2[i] + beta3[i]*j

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = antlerin,
                  ageclass = ageclass,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualby
)

#inits function
inits<- function(){list(beta1 = rnorm(11,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #interaction
                        beta4 = rnorm(1,0,1), # birthyear rainfall
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(antlers.jags)
MCMCtrace(antlers.jags)
write.csv(antlers.jags$summary, './output/antlers.sitexage.cat.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)


plot<- gather %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = 'dodge')+
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/antlers.sitexage.cat.jpg', plot, width = 15, height = 10)



# ---- Model 7: weight ~ site + birth rain + age (continuous)+ site x rain ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.rain.jags')
cat('
model{

#priors

  beta1 ~ dnorm(0, 0.001) #age as continuous
  beta2 ~ dnorm(0, 0.001) #rain 

 

for (u in 1:3){ 
  beta3[u] ~ dnorm(0, 0.001) #birthsite
  beta4[u] ~ dnorm(0, 0.001)  # site x rain interaction
}


for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]   #ageclass continous
              + beta2*rain[i]   #rainfall continuous
              + beta3[bs[i]]    # site
              + beta4[bs[i]]*rain[i]  #interaction betw by site and rain
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #rainfall
      response[j,i] <- beta2*rain.sim[j] + beta3[i] + beta4[i]*rain.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = weightlb,
                  ageclass = ageclass.sc,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualby,
                  rain.sim = annual.by.sim
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(1, 0, 1), #rain
                        beta3 = rnorm(3,0,1), #site
                        beta4 = rnorm(3,0,1), # interaction
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = runif(1,0,10))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

weight.rain.jags<- jagsUI(jags.data, inits, parameters, 'weight.rain.jags', n.thin = nt, n.chains = nc,
                      n.burnin = nb, n.iter = ni, parallel = T)
print(weight.rain.jags)
MCMCtrace(weight.rain.jags)
write.csv(weight.rain.jags$summary, './output/weight.byrainxsite.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.rain.jags %>% gather_draws(response[rain,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (annual.by.sim * sd(data$annual.by, na.rm = T)) + mean(data$annual.by, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$byrainfall <- vector

plot<- gather %>%
  ggplot(aes(x=byrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Birthyear Rainfall (in)", y = "Weight (lbs)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/byrainxweight.jpg', plot, width = 15, height = 10)


# ---- Model 8: weight ~ site + cap year rain + age (continuous)+ site x rain ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.rain.jags')
cat('
model{

#priors

  beta1 ~ dnorm(0, 0.001) #age as continuous
  beta2 ~ dnorm(0, 0.001) #rain 

 

for (u in 1:3){ 
  beta3[u] ~ dnorm(0, 0.001) #birthsite
  beta4[u] ~ dnorm(0, 0.001)  # site x rain interaction
}


for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]   #ageclass continous
              + beta2*rain[i]   #rainfall continuous
              + beta3[bs[i]]    # site
              + beta4[bs[i]]*rain[i]  #interaction betw by site and rain
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #rainfall
      response[j,i] <- beta2*rain.sim[j] + beta3[i] + beta4[i]*rain.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = weightlb,
                  ageclass = ageclass.sc,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualcy,
                  rain.sim = annual.cy.sim
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(1, 0, 1), #rain
                        beta3 = rnorm(3,0,1), #site
                        beta4 = rnorm(3,0,1), # interaction
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = runif(1,0,10))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

weight.rain.jags<- jagsUI(jags.data, inits, parameters, 'weight.rain.jags', n.thin = nt, n.chains = nc,
                          n.burnin = nb, n.iter = ni, parallel = T)
print(weight.rain.jags)
MCMCtrace(weight.rain.jags)
write.csv(weight.rain.jags$summary, './output/weight.cyrainxsite.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.rain.jags %>% gather_draws(response[rain,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (annual.cy.sim * sd(data$annual.cy, na.rm = T)) + mean(data$annual.cy, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$cyrainfall <- vector

plot<- gather %>%
  ggplot(aes(x=cyrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Capture Year Rainfall (in)", y = "Weight (lbs)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.1),          # x, y inside the plot area
        legend.justification = c("right", "bottom"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/cyrainxweight.jpg', plot, width = 15, height = 10)


# ---- Model9: antlers ~ site + birth rain + age (continuous)+ site x rain ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('antlers.rain.jags')
cat('
model{

#priors

  beta1 ~ dnorm(0, 0.001) #age as continuous
  beta2 ~ dnorm(0, 0.001) #rain 

 

for (u in 1:3){ 
  beta3[u] ~ dnorm(0, 0.001) #birthsite
  beta4[u] ~ dnorm(0, 0.001)  # site x rain interaction
}


for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]   #ageclass continous
              + beta2*rain[i]   #rainfall continuous
              + beta3[bs[i]]    # site
              + beta4[bs[i]]*rain[i]  #interaction betw by site and rain
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #rainfall
      response[j,i] <- beta2*rain.sim[j] + beta3[i] + beta4[i]*rain.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = antlerin,
                  ageclass = ageclass.sc,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualby,
                  rain.sim = annual.by.sim
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(1, 0, 1), #rain
                        beta3 = rnorm(3,0,1), #site
                        beta4 = rnorm(3,0,1), # interaction
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = runif(1,0,10))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

antlers.rain.jags<- jagsUI(jags.data, inits, parameters, 'antlers.rain.jags', n.thin = nt, n.chains = nc,
                          n.burnin = nb, n.iter = ni, parallel = T)
print(antlers.rain.jags)
MCMCtrace(antlers.rain.jags)
write.csv(antlers.rain.jags$summary, './output/antlers.byrainxsite.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.rain.jags %>% gather_draws(response[rain,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (annual.by.sim * sd(data$annual.by, na.rm = T)) + mean(data$annual.by, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$byrainfall <- vector

plot<- gather %>%
  ggplot(aes(x=byrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Birthyear Rainfall (in)", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/byrainxantlers.jpg', plot, width = 15, height = 10)


# ---- Model 10: antlers ~ site + cap year rain + age (continuous)+ site x rain ----

#basic model looking at the effect of site on morpho
set.seed(100)
sink('antlers.rain.jags')
cat('
model{

#priors

  beta1 ~ dnorm(0, 0.001) #age as continuous
  beta2 ~ dnorm(0, 0.001) #rain 

 

for (u in 1:3){ 
  beta3[u] ~ dnorm(0, 0.001) #birthsite
  beta4[u] ~ dnorm(0, 0.001)  # site x rain interaction
}


for (u in 1:11){ #random effect for capture year
  eps1[u] ~ dnorm(0,sigma)
}


for (u in 1:377){ #random effect for animal-id
  eps2[u] ~ dnorm(0, sigma)
}


#hyperparameters
  sigma ~ dunif(0,100)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]   #ageclass continous
              + beta2*rain[i]   #rainfall continuous
              + beta3[bs[i]]    # site
              + beta4[bs[i]]*rain[i]  #interaction betw by site and rain
              + eps1[capyear[i]]
              + eps2[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #rainfall
      response[j,i] <- beta2*rain.sim[j] + beta3[i] + beta4[i]*rain.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = antlerin,
                  ageclass = ageclass.sc,
                  capyear=capyear,
                  animal_id = animal_id,
                  rain = annualcy,
                  rain.sim = annual.cy.sim
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(1, 0, 1), #rain
                        beta3 = rnorm(3,0,1), #site
                        beta4 = rnorm(3,0,1), # interaction
                        eps1 = rnorm(11, 0, 1), #capture year
                        eps2 = rnorm(377,0,1), #animal id
                        sigma = runif(1,0,10))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4','eps1','eps2', 'response')  #'beta3',

#MCMC settings
ni <- 10000
nb <- 5000
nt<- 10
nc <- 3

antlers.rain.jags<- jagsUI(jags.data, inits, parameters, 'antlers.rain.jags', n.thin = nt, n.chains = nc,
                          n.burnin = nb, n.iter = ni, parallel = T)
print(antlers.rain.jags)
MCMCtrace(antlers.rain.jags)
write.csv(antlers.rain.jags$summary, './output/antlers.cyrainxsite.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.rain.jags %>% gather_draws(response[rain,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 4500 values of antler 1

#unscale and uncenter rain.sim
rain.sim.usc <- (annual.cy.sim * sd(data$annual.cy, na.rm = T)) + mean(data$annual.cy, na.rm = T)

#create vector containing simulated antler data but in the format to sync up with gather
vector <- numeric(0)
rain.sim.usc1 <- for (i in rain.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$cyrainfall <- vector

plot<- gather %>%
  ggplot(aes(x=cyrainfall, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Capture Year Rainfall (in)", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/cyrainxantlers.jpg', plot, width = 15, height = 10)











# Old models #############################  
#basic model looking at the effect of site on morpho
set.seed(100)
sink('antlers.jags')
cat('
model{

#priors

beta1 ~ dnorm(0, 0.001) #age as continuous
          

for (u in 1:3){ 
  beta2[u] ~ dnorm(0, 0.001) #birthsite
  beta3[u] ~ dnorm(0, 0.001)  # site x age interaction
}



# 
# for (u in 1:11){ #random effect for capture year
#   eps1[u] ~ dnorm(0,sigma)
# }

#
# for (u in 1:494){ #random effect for animal-id
#   eps2[u] ~ dnorm(0, sigma)
# }


#hyperparameters
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 y[i] ~ dnorm(mu[i], tau) #each morpho is a draw from this distribution
 mu[i] <-   beta1*ageclass[i]  #ageclass continuous
              + beta2[bs[i]]      #birthsite categorical 
              + beta3[bs[i]]*ageclass[i]  #interaction betw by age and site
              
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #ageclass
      response[j,i] <- beta1*ageclass.sim[j] + beta2[i] + beta3[i]*ageclass.sim[j]

      }
    }


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, y = antlerin,
                  ageclass = ageclass.sc,
                  ageclass.sim = ageclass.sim,
                  capyear=capyear,
                  animal_id = animal_id,
                  match=match
)

#inits function
inits<- function(){list(beta1 = rnorm(1,0,1), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #interaction
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3','response')  #'beta3',

#MCMC settings
ni <- 5000
nb <- 2000
nt<- 10
nc <- 3

antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc,
                      n.burnin = nb, n.iter = ni, parallel = T)
print(antlers.jags)
MCMCtrace(antlers.jags)
write.csv(antlers.jags$summary, './output/antlers.agexsite.cont.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$age == 2)[1] 

#unscale and uncenter weight
age.sim.usc <- (ageclass.sim * sd(data$age, na.rm = T)) + mean(data$age, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
age.sim.usc1 <- for (i in age.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$ageclass <- vector

plot<- gather %>%
  ggplot(aes(x=ageclass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/antlers.agexsite.cont.jpg', plot, width = 15, height = 10)







#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- antlers.jags %>% gather_draws(response[age,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$age == 2)[1] 

#unscale and uncenter weight
age.sim.usc <- (ageclass.sim * sd(data$age, na.rm = T)) + mean(data$age, na.rm = T)

#create vector containing simulated morpho data but in the format to sync up with gather
vector <- numeric(0)
age.sim.usc1 <- for (i in age.sim.usc) {
  rep_i <- rep(i, times = first_idx-1) #change times to match the number of first_idx
  vector <- c(vector,rep_i)
  
}

gather$ageclass <- vector

plot<- gather %>%
  ggplot(aes(x=ageclass, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(option = 'turbo', alpha = .2, labels = c("DMP", "CONTROL", "TGT") ) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(option = 'turbo', labels = c("DMP", "CONTROL", "TGT"))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers (in)", title = "")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "inside",
        legend.position.inside = c(0.1,0.9),          # x, y inside the plot area
        legend.justification = c("left", "top"),        # anchor point of the legend box        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 32, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 28, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
plot
ggsave('./figures/antlers.agexsite.cont.jpg', plot, width = 15, height = 10)



plot<- gather  %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = 'dodge')+
  # stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('Match', 'BY>CY', 'BY<CY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d( labels = c('Match', 'BY>CY', 'BY<CY'))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers", title = "Antlers by Age and Env Match")+
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

plot
# ggsave('.//figures/weight.rain.site.intrx.jpg', plot, width = 10, height = 8)






#basic model looking at the effect of site on morpho
set.seed(100)
sink('weight.jags')
cat('
model{

#priors

for (u in 1:12){ #age as categorical
  beta1[u] ~ dnorm( 0 , 0.001 ) 
}

for (u in 1:3){ #birthsite
  beta2[u] ~ dnorm(0,0.001)
}

for (u in 1:3){ #environmental matching
  beta3[u] ~ dnorm(0,0.0001)
}

for (u in 1:3){ #interaction between match and site
  beta4[u] ~ dnorm(0,0.0001)
}

# 
# for (u in 1:11){ #random effect for capture year
#   eps1[u] ~ dnorm(0,sigma)
# }

#
# for (u in 1:494){ #random effect for animal-id
#   eps2[u] ~ dnorm(0, sigma)
# }
# 
# for (u in 1:3){ #rain birthsite interaction
#   beta3[u] ~ dnorm( 0 , 0.001 )
# }

#hyperparameters
sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 morpho[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- beta1[ageclass[i]]   #ageclass categorical
              + beta2[bs[i]]  #birthsite categorical 
              + beta3[match[i]] #environmental match
              + beta4[bs[i]]*match[i]  #interaction betw by match and site
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:11){ #ageclass
      for (k in 1:3){ #environmental match
      response[j,i,k] <- beta1[j] + beta2[i] + beta3[k]+ beta4[i]*k

      }
    }}


  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, morpho = antlerin,
                  ageclass = ageclass,
                  capyear=capyear,
                  animal_id = animal_id,
                  match=match
                  )

#inits function
inits<- function(){list(beta1 = rnorm(12,0,10), #age
                        beta2 = rnorm(3, 0, 1), #site
                        beta3 = rnorm(3,0,1), #env match
                        beta4 = rnorm(3,0,1), 
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('beta1', 'beta2', 'beta3', 'beta4', 'response')  #'beta3',

#MCMC settings
ni <- 5000
nb <- 2000
nt<- 1
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni, parallel = T)
print(weight.jags)
MCMCtrace(weight.jags)

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(response[age,site,match])

gather$site <- as.factor(gather$site)
gather$age <- as.factor(gather$age)
gather$match<- as.factor(gather$match)

plot<- gather  %>%
  ggplot(aes(x=age, y=.value, color = site, fill = site)) +
  stat_pointinterval(position = 'dodge')+
  facet_wrap('match')+
  # stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('Match', 'BY>CY', 'BY<CY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d( labels = c('Match', 'BY>CY', 'BY<CY'))+ #color of line but no opacification
  labs(x = "Age", y = "Antlers", title = "Antlers by Age and Env Match")+
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

plot
# ggsave('.//figures/weight.rain.site.intrx.jpg', plot, width = 10, height = 8)











#weight ~ int + age + site + by.rain + 1|capyear
set.seed(100)
sink('weight.jags')
cat('
model{

#priors
int ~ dunif(0,100)
beta1[1] <- 0 #age
beta2[1] <- 0 #site
beta4[1] <- 0 #site * rain

for (u in 2:12){ #ageclass categorical
  beta1[u] ~ dnorm(0, 0.001)
}

# beta1 ~ dunif( 0 , 10 ) #age continuous

for (u in 2:3){ #birthsite
  beta2[u] ~ dnorm(0,0.001)
}

beta3 ~ dnorm(0, 0.001) #birth year rain, continuous


for (u in 2:3){ #rain birthsite interaction
  beta4[u] ~ dnorm( 0 , 0.001 )
}

#hyperparameters
sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 morpho[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- int + beta1[ageclass[i]]   #ageclass categorical
              + beta2[bs[i]]  #birthsite categorical 
              + beta3*rain[i] #continuous birth year rain
              + beta4[bs[i]]*rain[i]  #interaction betw by rain and site
              # + eps1[capyear[i]]
              # + eps.id[animal_id[i]]
 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #birth rain sim
      for (k in 7){ #ageclass 7
      response[j,i] <- int + beta1[k] + 
                              beta2[i] + 
                              beta3*rain.sim[j] +
                              beta4[i]*rain.sim[j]
                                    

      }
    }
}

  # for (i in 1:3){ #birthsite
  #   site_diff[i] <- beta2[i] - beta2[1]
  # }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, rain = annualcy, morpho = weightlb,
                                                      ageclass = ageclass,
                                                      capyear=capyear,
                                                      animal_id = animal_id,
                                                      rain.sim = annual.cy.sim
                  )

#inits function
inits<- function(){list(int = runif(1,0,100), 
                        beta1 = c(NA, rnorm(11,0,1)), #age
                        beta2 = c(NA, rnorm(2, 0, 1)), #site
                        beta3 = rnorm(1,0,1), #rain
                        beta4 = c(NA, rnorm(2,0,1)), #site x rain
                        sigma = rlnorm(1))}

#parameters to estimate
parameters <- c('int','beta1', 'beta2', 'beta3', 'beta4', 'response')  #, 'rain.bs.beta'

#MCMC settings
ni <- 5000
nb <- 2000
nt<- 1
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                       n.burnin = nb, n.iter = ni, parallel = T)
print(weight.jags)
MCMCtrace(weight.jags)

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(response[rain,site])

gather$site <- as.factor(gather$site)

#find first row for 2nd rain value
first_idx <- which(gather$rain == 2)[1] # 27000 values of rain 1

# unscale and uncenter rain.sim, for BIRTH YEAR
# rain.sim <- (annual.by.sim * sd(data$annual.by)) + mean(data$annual.by)
rain.sim <- (annual.cy.sim * sd(data$annual.cy)) + mean(data$annual.cy)


#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in rain.sim) {
  rep_i <- rep(i, times = 27000)
  vector1 <- c(vector1,rep_i)

}

gather$rain.sim <- vector1



plot<- gather  %>%
  ggplot(aes(x=rain.sim, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('DMP', 'EY', 'WY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('DMP', 'EY', 'WY'))+ #color of line but no opacification
  labs(x = "Cap Year ANNUAL RAINFALL", y = "Weight", title = "Weight by SITE*CY Rain")+
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

plot
# ggsave('.//figures/weight.rain.site.intrx.jpg', plot, width = 10, height = 8)

# # 




# # ########################################################################################
# #top performing model from survival data
# #antlers ~ int + age + site + by.rain + 1|capyear
# set.seed(100)
# sink('antlers.jags')
# cat('
# model{
# 
# #priors
# 
# age.beta ~ dunif( 0 , 10 )
# 
# 
# for (u in 1:3){ #birthsite
#   bs.beta[u] ~ dnorm(0,0.001)
# }
# 
# rain.by.beta ~ dnorm(0, 0.001) #birth year rain, continuous
# 
# for (u in 1:15){ #random effect for capture year
#   eps.capyear[u] ~ dnorm(0,sigma)
# }
# 
# for (u in 1:494){ #random effect for animal-id
#   eps.id[u] ~ dnorm(0, sigma)
# }
# 
# # for (u in 1:3 ){ #birth year rain and site intrx
# #   rain.bs.beta[u] ~ dnorm ( 0 , 0.01 )
# # }
# 
# #hyperparameters
# sigma ~ dunif(0,10)
# tau <- 1/(sigma*sigma)
# 
# # Likelihood
# for (i in 1:n){
#  antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
#  mu[i] <- age.beta*ageclass[i] + bs.beta[bs[i]]
#                                 + rain.by.beta*annualby[i]
#                                 + eps.capyear[capyear[i]]
#                                 + eps.id[animal_id[i]]
#                                # +rain.bs.beta[bs[i]]*annualby[i]
# 
#  }
# 
# #derived parameter
#   for (i in 1:3){ #birthsite
#     for (j in 1:100){ #birth rain sim
#       bcs[j,i] <- bs.beta[i] + rain.by.beta*annual.by.sim[j] 
#                              # + rain.bs.beta[i]*annual.by.sim[j]
# }
#   }
# 
#   
#     for (i in 1:3){ #birthsite
#     site_diff[i] <- bs.beta[i] - bs.beta[1]
#     }
# 
# 
# }
# ',fill = TRUE)
# sink()
# 
# #bundle data
# jags.data <- list(n=n, bs = bs, annualby = annualby, antlers = antlerin,
#                   ageclass = ageclass,
#                   capyear=capyear,
#                   animal_id = animal_id,
#                   annual.by.sim = annual.by.sim
# )
# 
# #inits function
# inits<- function(){list(bs.beta = rnorm(3, 0, 1), rain.by.beta = rnorm(1,0,1),
#                         sigma = rlnorm(1))}
# #log normal pulls just positive values,interaction.beta = rnorm(9, 0, 1),age.beta = rbinom(2226,1,1
# 
# #parameters to estimate
# parameters <- c('bs.beta', 'rain.by.beta', 'age.beta', 'site_diff', 'bcs')#'age.beta', , 'bcs'
# 
# #MCMC settings
# ni <- 1000
# nb <- 500
# nt<- 1
# nc <- 3
# 
# antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc,
#                      n.burnin = nb, n.iter = ni)
# print(antlers.jags)
# write.csv(antlers.jags$summary, './/output/antlers.rain.site.csv')
# #gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
# #then index in the order from output so above was bcs[j,i,k], can rename accordingly
# gather.bcs<- antlers.jags %>% gather_draws(bcs[by.rain,site])
# 
# #filter for min, max, and median birthyear rain values
# # gather<-gather %>% filter(by.rain %in% c(1,50,100))
# # gather$by.rain <- as.factor(gather$by.rain)
# gather.bcs$site <- as.factor(gather.bcs$site)
# 
# #find first row for 2nd rain value
# first_idx <- which(gather.bcs$by.rain == 2)[1] # 4500 values of rain 1
# 
# # unscale and uncenter rain.sim
# annual.by.sim.1 <- (annual.by.sim * sd(data$annual.by)) + mean(data$annual.by)
# # rain.by.sim1 <- (rain.by.sim * sd(data$annual.by)) + mean(data$annual.by)
# 
# #create vector containing simulated rainfall data but in the format to sync up with gather
# vector1 <- numeric(0)
# rain.sim3 <- for (i in annual.by.sim.1) {
#   rep_i <- rep(i, times = 4500) #adjust based upon first_idx
#   vector1 <- c(vector1,rep_i)
# 
# }
# 
# gather.bcs$annual.by.sim <- vector1
# 
# 
# 
# plot1<- gather.bcs  %>%
#   ggplot(aes(x=annual.by.sim, y=.value, color = site, fill = site)) +
#   stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
#   scale_fill_viridis_d(alpha = .2, labels = c('DMP', 'EY', 'WY')) + #this allowed me to opacify the ribbon but not the line
#   scale_color_viridis_d(labels = c('DMP', 'EY', 'WY'))+ #color of line but no opacification
#   labs(x = "BIRTHYEAR ANNUAL RAINFALL", y = "ANTLER SCORE (IN)", title = "ANTLERS BY SITE AND BIRTH RAIN")+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         legend.position = c(0.15,0.9),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 28),
#         plot.title = element_text(face = 'bold', size = 21, hjust = 0.5),
#         axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
#         axis.text = element_text(face='bold',size = 28),
#         # axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)
# 
# ggsave('.//figures/antlers.rain.bs.intrx.jpg', plot1, width = 10, height = 8)

