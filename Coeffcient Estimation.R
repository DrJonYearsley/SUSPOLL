Power_Input = c(1.84,1.64,1.54,1.54,  1.54,1.54,1.54,1.54,1.37,  0.53,0.45,0.40,0.40,  0.40,0.40,  
                0.53,0.40,0.40,  0.40,0.40,0.40,0.24,0.24)*0.0697333333   #convert from cal/min to Watts; power applied to bees (to produce heat)A = 4*pi*(c()/2)^2  #m^2; surface area of thorax (assuming spherical, based on average thorax diam)
# T_s = c() +273.15   #bee surface temp in K
# T_e = c() +273.15   #air temp in K
deltaT = c(11.3,10.2,9.1,11.4,  10.8,12.2,11.1,10.6,9.2,  7.2,6.5,7.2,6.9,  6.3,5.6,  5.5,5.8,7.3,  3.9,5.6,5.1,4.9,5.3) +273.15
D = c(8.8,8.5,8.3,8.3,  8.3,8.3,8.3,8.3,8.0,  5.8,5.5,5.3,5.3,5.3,5.3,  5.8,5.3,5.3,  5.3,5.3,5.3,4.5,4.5)/1000 #avg thorax diam, converted to m
#The preceding all from Church1960

r = D/2
A_cyl = 2*pi*r*D+2*pi*r^2
A = 4*pi*r^2  #m^2; bee surface area assuming a sphere, based on average thorax diam measurements
h_c = (Power_Input/A)/(deltaT) 
mean(h_c)

# Estimate n and C_l
V = 3 #m/s; velocity
k = 0.025  #W/mK; assuming air temp 20C, thermal conductivity of air from https://www.me.psu.edu/cimbala/me433/Links/Table_A_9_CC_Properties_of_Air.pdf
nu =  1.470*10^(-5)  #m^2/s; kinematic viscosity of air at 20C, same source

# T_e_mean = mean(T_e)    #use these average temperatures to calculate the fluid properties k and mu
# T_s_mean = mean(T_s)

Nus = h_c*D/k    #Nusselt number
Rey = V*D/nu     #Reynolds number


plot(log(Rey),log(Nus))
mod1=lm(log(Nus)~log(Rey))
abline(mod1)
summary(mod1)
C_l = exp(mod1$coefficients[1])
se_logC_l = 0.054062
C_l_lower = exp(mod1$coefficients[1]-se_logC_l)
C_l_upper = exp(mod1$coefficients[1]+se_logC_l)

n = mod1$coefficients[2]
se_n = 0.007548  #from summary
h  = (C_l*k/D)*(V*D/nu)^n

################################################
##### Metabolic rate #####

#Kammer1974: 
  # 55-66 ml O2/g_body/hr tethered flight
  # estimated 166-188 ml O2/g_th/hr free flight (based on mlO2/spike and spike rate)
  #resting 1.3 ml O2/g_body/hr

#Heinrich1975 Thermoregulation in BB II
  #150-350 ml O2/g_th/hr for B.vosnesenskii (increasing abdomen weight)

#Dzialowski
  #mean 11.25ish, range 9.5-13 mlO2/hr
  #not flying, but no recording of whether shivering


#I_flying = (21.117J/mlO2)*(ml_O2/g?/hr)*(1hr/60min)*(1min/60s)*g?  where? is thorax or body
I_Kammer_tetheredflight = 21.117*((55+66)/2)*(1/60)*(1/60)*0.149 
I_Kammer_freeflight = 21.117*((166+188)/2)*(1/60)*(1/60)*0.057
I_Heinrich = 21.117*((150+350)/2)*(1/60)*(1/60)*0.057
I_Kammer_resting = 21.117*1.3*(1/60)*(1/60)*0.149 
I_Dzialowski = 21.117*11.25*(1/60)*(1/60)
I_flying_min = 21.117*150*(1/60)*(1/60)*0.057
I_flying_max = 21.117*350*(1/60)*(1/60)*0.057
I_flying_min_Kammer = 21.117*166*(1/60)*(1/60)*0.057
I_flying_max_Kammer = 21.117*188*(1/60)*(1/60)*0.057

E = 0.63
k = 8.617*10^(-5)
i0_flying = exp(log(I_Kammer_freeflight) - (3/4)*log(0.149) + E/(k*(25+273.15)))
