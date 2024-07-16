#final analysis
#linear model looking at the effect of age and site interaction on antlers
library(jagsUI)
library(ggplot2)
library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(here)


#Load dataset with no fawns, just animals with antlers
# data <- read.csv('./raw/bucks_nofawns.csv', header=T)
data<- read.csv('./clean/nofawns22rain.csv', header=T)
head(data)

#replace zeros with NA
data[data==0] <- NA

#create a vector with 1 for treatment, 2 for control, 3 for tgt
bs<- as.numeric(factor(data$birthsite))

#create vectors of things of interest
data$age <- data$year_cap - data$year_birth 
ageclass <- data$age
ageclass <- as.vector(scale(ageclass)) #scale and center ageclass to treat as continuous cov.

#create a vector for capture year, to fit random effect later
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


#top performing model from survival data
#weight ~ int + age + site + by.rain + 1|capyear
set.seed(100)
sink('weight.jags')
cat('
model{

#priors

age.beta ~ dunif( 0 , 10 )

for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

rain.by.beta ~ dnorm(0, 0.001) #birth year rain, continuous

for (u in 1:15){ #random effect for capture year
  eps.capyear[u] ~ dnorm(0,sigma)
}

for (u in 1:494){ #random effect for animal-id
  eps.id[u] ~ dnorm(0, sigma)
}
# 
# for (u in 1:3){ #rain birthsite interaction
#   rain.bs.beta[u] ~ dnorm( 0 , 0.001 )
# }

#hyperparameters
sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 weight[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- age.beta*ageclass[i] + bs.beta[bs[i]] 
                                + rain.by.beta*annualby[i] 
                                + eps.capyear[capyear[i]]
                                + eps.id[animal_id[i]]
                                # + rain.bs.beta[bs[i]]*annualby[i]
 
 }                                   

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #birth rain sim
      bodymass[j,i] <- bs.beta[i] + rain.by.beta*annual.by.sim[j] 
                                    # + rain.bs.beta[i]*annual.by.sim[j]

      }
    }
  

  for (i in 1:3){ #birthsite
    site_diff[i] <- bs.beta[i] - bs.beta[1]
  }
  

}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, annualby = annualby, weight = weightlb, 
                                                      ageclass = ageclass, 
                                                      capyear=capyear,
                                                      animal_id = animal_id,
                                                      annual.by.sim = annual.by.sim
                  )

#inits function
inits<- function(){list(bs.beta = rnorm(3, 0, 1), age.beta = runif(1,0,10), rain.by.beta = rnorm(1,0,1),
                        sigma = rlnorm(1))}
#log normal pulls just positive values,interaction.beta = rnorm(9, 0, 1),age.beta = rbinom(2226,1,1

#parameters to estimate
parameters <- c('bs.beta', 'rain.by.beta', 'age.beta',  'bodymass')  #, 'rain.bs.beta'

#MCMC settings
ni <- 1000
nb <- 500
nt<- 1
nc <- 3

weight.jags<- jagsUI(jags.data, inits, parameters, 'weight.jags', n.thin = nt, n.chains = nc,
                       n.burnin = nb, n.iter = ni)
print(weight.jags)
write.csv(weight.jags$summary, './/output/weight.rain.site.csv')

#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather<- weight.jags %>% gather_draws(bodymass[by.rain,site])
#
#filter for min, max, and median birthyear rain values
# gather<-gather %>% filter(by.rain %in% c(1,50,100))
# gather$by.rain <- as.factor(gather$by.rain)
gather$site <- as.factor(gather$site)
# gather$age <- as.factor((gather$age))

#find first row for 2nd rain value
first_idx <- which(gather$by.rain == 2)[1] # 4500 values of rain 1

# unscale and uncenter rain.sim
annual.by.sim.1 <- (annual.by.sim * sd(data$annual.by)) + mean(data$annual.by)
# rain.by.sim1 <- (rain.by.sim * sd(data$annual.by)) + mean(data$annual.by)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in annual.by.sim.1) {
  rep_i <- rep(i, times = 4500)
  vector1 <- c(vector1,rep_i)

}

gather$annual.by.sim <- vector1



plot<- gather  %>%
  ggplot(aes(x=annual.by.sim, y=.value, color = site, fill = site)) +
  # facet_wrap(vars(age), nrow = 1)+
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('DMP', 'EY', 'WY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('DMP', 'EY', 'WY'))+ #color of line but no opacification
  labs(x = "BIRTHYEAR ANNUAL RAINFALL", y = "WEIGHT (LBS)", title = "WEIGHT BY SITE AND BIRTH RAIN")+
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

ggsave('.//figures/weight.rain.site.intrx.jpg', plot, width = 10, height = 8)
# # 
# ########################################################################################
#top performing model from survival data
#antlers ~ int + age + site + by.rain + 1|capyear
set.seed(100)
sink('antlers.jags')
cat('
model{

#priors

age.beta ~ dunif( 0 , 10 )


for (u in 1:3){ #birthsite
  bs.beta[u] ~ dnorm(0,0.001)
}

rain.by.beta ~ dnorm(0, 0.001) #birth year rain, continuous

for (u in 1:15){ #random effect for capture year
  eps.capyear[u] ~ dnorm(0,sigma)
}

for (u in 1:494){ #random effect for animal-id
  eps.id[u] ~ dnorm(0, sigma)
}

# for (u in 1:3 ){ #birth year rain and site intrx
#   rain.bs.beta[u] ~ dnorm ( 0 , 0.01 )
# }

#hyperparameters
sigma ~ dunif(0,10)
tau <- 1/(sigma*sigma)

# Likelihood
for (i in 1:n){
 antlers[i] ~ dnorm(mu[i], tau) #each antler is a draw from this distribution
 mu[i] <- age.beta*ageclass[i] + bs.beta[bs[i]]
                                + rain.by.beta*annualby[i]
                                + eps.capyear[capyear[i]]
                                + eps.id[animal_id[i]]
                               # +rain.bs.beta[bs[i]]*annualby[i]

 }

#derived parameter
  for (i in 1:3){ #birthsite
    for (j in 1:100){ #birth rain sim
      bcs[j,i] <- bs.beta[i] + rain.by.beta*annual.by.sim[j] 
                             # + rain.bs.beta[i]*annual.by.sim[j]
}
  }

  
    for (i in 1:3){ #birthsite
    site_diff[i] <- bs.beta[i] - bs.beta[1]
    }


}
',fill = TRUE)
sink()

#bundle data
jags.data <- list(n=n, bs = bs, annualby = annualby, antlers = antlerin,
                  ageclass = ageclass,
                  capyear=capyear,
                  animal_id = animal_id,
                  annual.by.sim = annual.by.sim
)

#inits function
inits<- function(){list(bs.beta = rnorm(3, 0, 1), rain.by.beta = rnorm(1,0,1),
                        sigma = rlnorm(1))}
#log normal pulls just positive values,interaction.beta = rnorm(9, 0, 1),age.beta = rbinom(2226,1,1

#parameters to estimate
parameters <- c('bs.beta', 'rain.by.beta', 'age.beta', 'site_diff', 'bcs')#'age.beta', , 'bcs'

#MCMC settings
ni <- 1000
nb <- 500
nt<- 1
nc <- 3

antlers.jags<- jagsUI(jags.data, inits, parameters, 'antlers.jags', n.thin = nt, n.chains = nc,
                     n.burnin = nb, n.iter = ni)
print(antlers.jags)
write.csv(antlers.jags$summary, './/output/antlers.rain.site.csv')
#gatherdraws creates a dataframe in long format, need to subset by the variable of interest in jags output,
#then index in the order from output so above was bcs[j,i,k], can rename accordingly
gather.bcs<- antlers.jags %>% gather_draws(bcs[by.rain,site])

#filter for min, max, and median birthyear rain values
# gather<-gather %>% filter(by.rain %in% c(1,50,100))
# gather$by.rain <- as.factor(gather$by.rain)
gather.bcs$site <- as.factor(gather.bcs$site)

#find first row for 2nd rain value
first_idx <- which(gather.bcs$by.rain == 2)[1] # 4500 values of rain 1

# unscale and uncenter rain.sim
annual.by.sim.1 <- (annual.by.sim * sd(data$annual.by)) + mean(data$annual.by)
# rain.by.sim1 <- (rain.by.sim * sd(data$annual.by)) + mean(data$annual.by)

#create vector containing simulated rainfall data but in the format to sync up with gather
vector1 <- numeric(0)
rain.sim3 <- for (i in annual.by.sim.1) {
  rep_i <- rep(i, times = 4500) #adjust based upon first_idx
  vector1 <- c(vector1,rep_i)

}

gather.bcs$annual.by.sim <- vector1



plot1<- gather.bcs  %>%
  ggplot(aes(x=annual.by.sim, y=.value, color = site, fill = site)) +
  stat_lineribbon(.width = 0.95)+ #statline ribbon takes posterior estimates and calculates CRI
  scale_fill_viridis_d(alpha = .2, labels = c('DMP', 'EY', 'WY')) + #this allowed me to opacify the ribbon but not the line
  scale_color_viridis_d(labels = c('DMP', 'EY', 'WY'))+ #color of line but no opacification
  labs(x = "BIRTHYEAR ANNUAL RAINFALL", y = "ANTLER SCORE (IN)", title = "ANTLERS BY SITE AND BIRTH RAIN")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.15,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        plot.title = element_text(face = 'bold', size = 21, hjust = 0.5),
        axis.title = element_text(face = 'bold',size = 32, hjust = 0.5),
        axis.text = element_text(face='bold',size = 28),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) #transparent plot bg)

ggsave('.//figures/antlers.rain.bs.intrx.jpg', plot1, width = 10, height = 8)

