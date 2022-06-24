T_aK = 15+273.15
y = 39+273.15   #thorax temp
r = 0.0367/9

h_fg = 2.33*10^6  #latent heat of vaporization of water, J/kg 
A = 9.1496*10^10  #clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3  #clausius-clapyron constant B for water in K
rh =  0.7   #ambient relative humidity - actually calculate this from Arrian's data
Pr = 1.013*10^5 #atmospheric pressure in N/m^2
R_0 = 0.001  #radius of nectar droplet, 1mm in m
v=3.1
MM_air = 0.0289652   #molar mass of dry air in kg/mol
MM_vapor = 0.018016  #molar mass of water vapor in kg/mol
delta_T_h = 2.9
R_specific = 287.058  #J/kg/K for dry air

air_h = rh
sfc_h = 1  #humidity at surface
#sfc_h = rh  #humidity at surface

Ev = numeric(51)
Ab = numeric(51)
Temps = 0:50

###################################
##### Estimating D_A

#using the Peclet number, UL/D = 1
U = v   #scale for speed, average flight speed
L = 2*R_0   #approximate domain length or size of release location 
Pe = (U*L)/(2.06*10^-5)


D = U*L
D_A = D  #diffusion coefficient of air into itself in m^2/s
#D_A = v/6  #diffusion coefficient of air into itself in m^2/s
#D_A = 2.06*10^-5  #diffusion coefficient of air into humid air in m^2/s



for(i in 0:50){
  
  y = i+273.15
  
#### abdomen cooling ####
Ab[i+1] = (y-T_aK)*r

#### evaporative cooling #####


p_sat_air = A*exp(B/T_aK)  #partial pressure of saturated vapor (humid air at saturation)
p_air = air_h*p_sat_air   #partial pressure of humid air at temperature T_aK
X_air = p_air/Pr  #mole fraction of ambient air
Y_air = 1/( 1 + ((1-X_air)/X_air)*(MM_air/MM_vapor) )

p_sat_sfc = A*exp(B/(y-delta_T_h))
p_sfc = sfc_h*p_sat_sfc
X_sfc = p_sfc/Pr  #mole fraction of surface air, at 100% humidity
Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) )

rho_humid = ((Pr-p_air)*MM_air + p_air*MM_vapor)/(R_specific*T_aK)   #density of humid air

m_evap_coef = 2*pi*R_0*rho_humid*D_A
Ev1 = h_fg*(m_evap_coef*log(1-Y_air))  #constant for a fixed T_air
Ev2 = -h_fg*m_evap_coef*log(1-Y_sfc)   #varies with T_th (or T_nectar)

Ev[i+1] = Ev1+Ev2
}

plot(Temps,Ab,pch=20,type="b",col=2,
     main="Comparison of Cooling Methods \n as a funtion of thorax temperature",
     xlab = "Thorax Temp (C)",ylab="Heat Flux (W)",
     ylim=c(-0.5,0.5))
points(Temps,Ev,pch=20,type="b",col=4)
legend(0,0.4,c("Abdomen cooling","Evaporative cooling"),pch=20,col=c(2,4))
abline(h=0,lty='dotted')



