setwd("~/SUSPOLL/Code")
#install.packages("ggplot2")
library(ggplot2)

masses = seq(from=0.01, by=0.0005, to=0.16)   #full range

MaxAirTemp_HB = read.csv('MaxAirTemp_HB.csv',header=F)
MaxAirTemp_BB = read.csv('MaxAirTemp_BB.csv',header=F)

Bees = data.frame(mass = masses, slowBB = MaxAirTemp_BB[,1], 
                 defaultBB = MaxAirTemp_BB[,2], 
                 fastBB = MaxAirTemp_BB[,3],
                 slowHB = MaxAirTemp_HB[,1], 
                 defaultHB = MaxAirTemp_HB[,2], 
                 fastHB = MaxAirTemp_HB[,3])



p = ggplot(Bees) + 
#  geom_point(aes(x=mass,y=slow)) +
  geom_smooth(aes(x=mass,y=slowBB),colour="blue") +
  geom_smooth(aes(x=mass,y=defaultBB),colour="blue") +
  geom_smooth(aes(x=mass,y=fastBB),colour="blue") +
  geom_smooth(aes(x=mass,y=slowHB),colour="red") +
  geom_smooth(aes(x=mass,y=defaultHB),colour="red") +
  geom_smooth(aes(x=mass,y=fastHB),colour="red") +
  geom_ribbon(aes(x = mass, ymin=slowHB,ymax=fastHB), fill="red", alpha=0.3) +
  geom_ribbon(aes(x = mass, ymin=slowBB,ymax=fastBB), fill="blue", alpha=0.3) +
  labs(x='mass of thorax (g)', 
       y='Maximum air temperature for sustained flying') + 
  theme_bw()
# p + annotate("text", x = 0.14, y = 30, label = "v = 4.5 m/s")
# p  + annotate("text", x = 0.14, y = 19, label = "v = 1 m/s")
# p  + annotate("text", x = 0.14, y = 34, label = "v = 5.5 m/s")
# p + annotate("text", x = 0.14, y = 30, label = "v = 4.5 m/s")
# p  + annotate("text", x = 0.14, y = 19, label = "v = 1 m/s")
# p  + annotate("text", x = 0.14, y = 34, label = "v = 5.5 m/s")


p + annotate("text", x = c(0.14,0.14,0.14,0.145,0.14,0.145), y = c(19,30,34,23,30,39), label = c("v=1 m/s","v = 4.5 m/s","v=5.5 m/s","v=3.6 m/s","v = 4.5 m/s","v=9 m/s"))

             