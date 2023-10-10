library(ggplot2)
data<-read.csv('./raw/testosteron.csv', header = T)

plot_base_lme <-
  ggplot(data = data, aes(x=month, y=testosterone))+
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


plot<- plot_base_lme +
  # stat_halfeye()+
  geom_line(linewidth = 3, color = 'red')+
  scale_x_discrete(limits=c('1', '2', '3', '4' ,'5' ,'6' ,'7','8','9','10','11','12'))+
                   # labels = c('1.5', '2.5', '3.5', '4.5' ,'5.5' ,'6.5' ,'7.5','8.5',
                   #            '9.5','10.5','11.5','12.5'))+
  labs(x = "MONTH", y = "TESTOSTERONE (ng/dL)",
       title = "TESTOSTERONE BY MONTH")

ggsave('./figures/testo.png', plot, bg='transparent', width = 16, height = 12)
