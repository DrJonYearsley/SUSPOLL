setwd("~/SUSPOLL/Code")


LookupTable_Bumblebee_min = read.csv('LookupTable_Bumblebee_min.csv',header=FALSE)
LookupTable_Bumblebee_default = read.csv('LookupTable_Bumblebee_default.csv',header=FALSE)
LookupTable_Bumblebee_max = read.csv('LookupTable_Bumblebee_max.csv',header=FALSE)

plot(LookupTable_Bumblebee_default[ ,6],type='l')
points(0:50,LookupTable_Bumblebee_min[ ,6],type='l')
points(0:50,LookupTable_Bumblebee_max[ ,6],type='l')


x_temp = c(0:50,50:0)
y_temp = c(LookupTable_Bumblebee_min[ ,6],rev(LookupTable_Bumblebee_max[ ,6]))


plot(x_temp,y_temp,type='l',xlab = 'Air Temperature (C)',ylab='Equilibrium Thorax Temperature (C)',ylim=c(0,50))
abline(h=45,col='red')
abline(h=42,col='red')
abline(h=30,col='blue')
polygon(x_temp,y_temp, col = rgb(0.53, 0.79, 0.88, alpha = 0.5))
points(0:50,LookupTable_Bumblebee_default[ ,6],type='l',lty=4)

