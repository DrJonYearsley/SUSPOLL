#T_th = ( (epsilon_a/epsilon_e)*(0.25*csc(beta)*P + alpha*a*P + alpha*delta*T_aK^6 + alpha*sigma*T_gK^4) )^(1/4) + 273.15 #in celsius

#default values
epsilon_a = 0.935
epsilon_e = 0.97
delta = 5.31*10^(-13)
sigma = 5.67*10^(-8)
alpha = 0.5

#Varying environmental values
#1
P = mean(ArrianDataSet$meansolarstation)      #solar radiations observed by Arrian = a good range of typical Irish weather
P_min = min(ArrianDataSet$meansolarstation)
P_max = max(ArrianDataSet$meansolarstation)

#2
T_aK = mean(ArrianDataSet$thermtemp)+273.15   #0-30C should cover it
T_aK_min = 0+273.15
T_aK_max = 30+273.15

#3
T_gK = 11+273.15        #from https://www.met.ie/forecasts/farming/agricultural-data-report May 25 2021
T_gK_min = 6+273.15     #https://www.farmersjournal.ie/soil-temperature-rise-to-help-growth-151968, https://www.farmersjournal.ie/dairy-management-return-to-winter-weather-612799
T_gK_max = 15+273.15       #http://edepositireland.ie/bitstream/handle/2262/71180/Agromet%20Memo%20No.%203.pdf?sequence=1&isAllowed=y

#4
beta_min = 0.1*(pi/180)        #sunrise... 
beta_max = 60.14*(pi/180)    #max reached on solstice (in Dublin) https://www.suncalc.org/#/53.3066,-6.2223,13/2021.06.22/13:35/1/3
beta = (48)*(pi/180)

#5
a = 0.25
a_min = 0.9*a   #unknown, so do +- 10%
a_max = 1.1*a
 
T_th_default = ( (epsilon_a/epsilon_e)*(0.25*(1/sin(beta))*P + alpha*a*P + alpha*delta*T_aK^6 + alpha*sigma*T_gK^4 ) )^(1/4) - 273.15
T_th_min = ( (epsilon_a/epsilon_e)*(0.25*(1/sin(beta_min))*P_min + alpha*a_min*P_min + alpha*delta*T_aK_min^6 + alpha*sigma*T_gK_min^4 ) )^(1/4) - 273.15
T_th_max = ( (epsilon_a/epsilon_e)*(0.25*(1/sin(beta_max))*P_max + alpha*a_max*P_max + alpha*delta*T_aK_max^6 + alpha*sigma*T_gK_max^4 ) )^(1/4) - 273.15