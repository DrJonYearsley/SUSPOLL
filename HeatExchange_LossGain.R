
delta = 5.31*10^(-13)
sigma = 5.67*10^(-8)
k = 8.617333262145*10^-5    #Boltzmann's constant
E = 0.63     #activation energy from Brown2004

T_aC = seq(13,40,1)    #air temp in C
T_sC = 9                  #ground surface temp in C, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
T_aK = T_aC+273.15    #air temp in K
T_sK = T_sC+273.15                  #ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report

T_thC = 39      ##lowest value for departing bees in Cooper1985
T_thK = T_thC+273.15


A_th = 4.5*10^(-5)     #in m^2, from Cooper1985
M_b = 0.100   #mass of the bee in g (the default used by Cooper1985 is 100mg)
M_th = 0.0407  #mass of thorax in g (mean reported in Cooper1985)

v = 3.1    #using M data for 100mg bee flying at 3m/s

c = 4.184*0.8   #specific heat (in cal/g*degC converted to J/g*degC), cited in May1976

### Calculate Some things ###

#solar radiation
P=850   #from Cooper1985 
S = (0.25*0.91*A_th*P) + (0.5*0.91*A_th*0.25*P)   #%solar radiation

#thermal radiation
R1 = 0.5*0.91*A_th*(delta*T_aK^6 + sigma*T_sK^4)
R2 = 0.97*A_th*sigma*T_thK^4  #temps need to be in K

#convective heat loss - theoretical 
kin.viscosity =  function(T){       #function of temp in K (nu)
  #read a table to extract the right entry?
  return(1.45*10^(-5))   #ish for 10-15C
}
therm.conduct =  function(T){       #function of temp in K (k)
  #read a table to extract the right entry?
  return(0.024)   #ish for 10-15C
}
C_l = 0.0749
n = 0.78	
l_th = 0.0053
h_flying = ((therm.conduct(T_aK)*C_l)/l_th)*((v*l_th)/(kin.viscosity(T_aK)))^n
h_resting = ((therm.conduct(T_aK)*C_l)/l_th)*((v*l_th)/(kin.viscosity(T_aK)))^n #assumes forced airflow, so probably shouldn't actually be used for "resting" bee
C_flying = h_flying*A_th*(T_thK-T_aK)
C_resting = h_resting*A_th*(T_thK-T_aK) 


#metabolic heat production 
#M = 0.037 from Cooper1985
I_resting = 1.31*10^(-3)     #From Rothe1989, in mW/g
i0_resting = exp(log(I_resting) - (3/4)* log(M_b) + E/(k*(10+273.15)))
I_flying = 0.4     #From Nachtigall1989, in W/g
i0_flying = exp(log(I_flying) - (3/4)* log(M_b) + E/(k*(20+273.15)))

#M = i_0*M_b^(3/4)*exp(-E/(k*T_aK))    #mass in g and temp in K

#M_resting = I_resting*M_b
#M_flying = I_flying*M_b
M_resting = i0_resting*M_b^(3/4)*exp(-E/(k*T_aK))*M_b    #mass in g and temp in K
M_flying = i0_flying*M_b^(3/4)*exp(-E/(k*T_aK))*M_b    #mass in g and temp in K


plot(T_aC,M_flying+rep(S,length(T_aC))+R1,ylim=c(-0.05,0.25),col='green',pch=20)  #gains
points(T_aC,C_flying+R2,col='red',pch=20)  #losses
abline(h=0,col='grey')

plot(T_aC,M_flying,col='green',pch=20)  #gains


# plot(T_aC,(M_flying+rep(S,length(T_aC))+R1)/(M_th*c),col='green',pch=20)  #gains
# points(T_aC,(C_flying+R2)/(M_th*c),col='red',pch=20)  #losses
# abline(h=0,col='grey')
