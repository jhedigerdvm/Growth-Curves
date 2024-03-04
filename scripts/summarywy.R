wy <- data %>% filter(birthsite=="wy")
wy<- drop_na(wy)
wy %>% group_by(age) %>%  summarise(mean = mean(bcsin), sd = sd(bcsin), count = n())
wy %>% group_by(age) %>%  summarise(mean = mean(weight), sd = sd(weight), count = n())
