### From textbook example fig 13.7 ###

# Input parameters 
D = 0.0893 #mug inner diameter in m
Pinf = 1.013*10^5 #atmospheric pressure in N/m^2
Tinf = 25+273 #ambient temperature in K
rh = 0.7  #ambient relative humidity 

# Physical constants
Dv = 2.56*10^-5  #diffusion coefficient of water vapor into air in m^2/s
DA = 2.06*10^-5  #diffusion coefficient of air into itself in m^2/s
hfg = 2.33*10^6  #heat of vaporization of water in J/kg
A = 9.1496*10^10  #clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3  #clausius-clapyron constant B for water in K

#physical constants added by me
MW_air = 28.97  #molar mass of dry air in g/mol
MW_vapor = 18.01528  #molar mass of water vapor in g/mol

#my intermediate calculations
p_sat = A*exp(B/Tinf)
p_v = rh*p_sat

p_v_sfc = p_sat   #at the surface, relative humidity is 100#

# Derived parameters
S = pi*(D/2)^2  #free surface area 0.00626315m^2
Xinf = p_v/Pinf #mol fraction of ambient 0.0222
Yinf = 1/( 1 + ((1-Xinf)/Xinf)*(MW_air/MW_vapor) )  #mass fraction of ambient 0.0138
rho
X0_sat = p_v_sfc/Pinf  #mol fraction at surface when coffee is ambient temp
Y0_sat = 1/( 1 + ((1-X0_sat)/X0_sat)*(MW_air/MW_vapor) )  #mass fraction at surface when coffee is at air temp
DewPoint




#### Exploring evaporative cooling #####

h_fg = 2.33*10^6  #latent heat of vaporization of water, J/kg 
A = 9.1496*10^10  #clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3  #clausius-clapyron constant B for water in K
rh =  0.7   #relative humidity - actually calculate this from Arrian's data
Pr = 1.013*10^5 #atmospheric pressure in N/m^2
R_0 = 0.001  #radius of nectar droplet, 1mm in m
D_A = 2.06*10^-5  #diffusion coefficient of air into itself in m^2/s
MM_air = 0.0289652   #molar mass of dry air in g/mol
MM_vapor = 0.018016  #molar mass of water vapor in kg/mol
delta_T_h = 2.9
R_specific = 287.058  #J/kg/K for dry air

T_aK = 15+273.15
y = 39+273.15   #thorax temp
sfc_h = 1  #humidity at surface

p_sat_air = A*exp(B/T_aK)  #partial pressure of saturated vapor (humid air at saturation)
p_air = rh*p_sat_air   #partial pressure of humid air at temperature T_aK
rho_humid = ((Pr-p_air)*MM_air + p_air*MM_vapor)/(R_specific*T_aK)   #density of humid air

X_air = p_air/Pr  #mole fraction of ambient air
Y_air = 1/( 1 + ((1-X_air)/X_air)*(MM_air/MM_vapor) )

X_sfc = (sfc_h*A*exp(B/(y-delta_T_h)))/Pr  #mole fraction of surface air, at 100% humidity
Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) )

m_evap_coef = 2*pi*R_0*rho_humid*D_A
Ev1 = h_fg*m_evap_coef*log(1-Y_air)   #constant for a fixed T_air
Ev2 = -h_fg*m_evap_coef*log(1-Y_sfc)

Ev = Ev1+Ev2


#calculate Ev as a function of sfc_h
Ev = numeric(10)
for(i in 1:10){
  sfc_h = i/10
  X_sfc = (sfc_h*A*exp(B/(y-delta_T_h)))/Pr  #mole fraction of surface air, at 100% humidity
  Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) )
  Ev2 = -h_fg*m_evap_coef*log(1-Y_sfc)
  Ev[i] = Ev1+Ev2
}
plot((1:10/10),Ev)

# calculate Ev as a function of D_A
Ev = numeric(10)
for(i in 1:10){
  D_A = i/10  #but make this something reasonable
  m_evap_coef = 2*pi*R_0*rho_humid*D_A
  Ev1 = h_fg*m_evap_coef*log(1-Y_air)   #constant for a fixed T_air
  Ev2 = -h_fg*m_evap_coef*log(1-Y_sfc)
  Ev[i] = Ev1+Ev2
}
plot((1:10/10),Ev)
