T_thC = 45
T_thK = T_thC+273.15
C_l = 2.429809*10^(-7)   #fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485       ##fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
l_th = 0.005467   #characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
v = 4.1
A_th = 9.3896*10^(-5)   #thorax surface area in m^2, from Church1960
#s = 0.9141875
s = 1
M_th = 0.057  #mass of thorax in g, Joos1991
c = 3.349  #specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
Pr = 1.013*10^5 #atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = 287.058  #J/kg/K for dry air

plot(0,0,ylim=c(-0.1,0.1),xlim=c(0,50),pch=".",
     main=('Convection \n as a function of air temperature'),
     ylab = 'Q_C', xlab = 'Air temperature (C)')
abline(h=0)
abline(v = s*(T_thK)-273.15,lty="dotted")
abline(v = 0)

for(i in 1:51){
  T_aC = i-1    #varying air temp
  T_aK = T_aC+273.15     #air temp in K
  mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4) #dynamic viscosity of air 
#  rho = 1.2256  #in dry air air at sea level at 15C
  rho = Pr/(R_specific*T_aK)  #in humid air air at sea level 
  #nu = 2.791*10^(-7)*T_aK^(0.7355)/1.2256   #kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
  nu = mu/rho #kinematic viscosity #The Shock Absorber Handbook
  kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK))

  h  = (C_l*kappa/l_th)*(v*l_th/nu)^n
  C1 = (h*A_th*s*T_thK)           #thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
  C2 = -(h*A_th*T_aK)           #thorax
  C = h*A_th*(s*T_thK-T_aK)
  points(T_aK-273.15,C1+C2,pch=20,col=2)

}
