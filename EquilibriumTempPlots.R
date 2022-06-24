CooperAirTemps = c(13, 18, 22, 28, 32, 36, 40)
CooperTempsIn = c(35.5, 34, 37, 40, 42, 43, 44)
EquilibriaResting = c(298.7166,  304.3341,  308.8684,  315.7475,  320.3936,  325.0957,  329.8629)

EquilibriaFlying = c(305.0420,  314.2029,  322.8167,  338.8385,  352.4270,  369.2817,  390.5608)


plot(CooperAirTemps,CooperTempsIn,pch=20,col=4,type='b',ylim=c(24,118),
     xlab='Air temperature', ylab='Returning (equilibrium) thorax temperature')
#abline(h=33,lty='dotted') #minimum thoracic temp for flight
points(CooperAirTemps,EquilibriaResting-273.15,pch=20,col=6,type='b')
points(CooperAirTemps,EquilibriaFlying-273.15,pch=20,col=2,type='b')
legend(13,120,legend=c('Cooper 1985 Data','Heat Exchange ODE (resting)','Heat Exchange ODE (flying)'),col=c(4,6,2),pch=20)
abline(h=39,col='gray')  #starting thorax temp
#abline(a=0,b=1) T_a = T_th line
