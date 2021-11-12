library(lhs)
setwd("~/SUSPOLL/Code")


# date = "SensitivityAnalysis"           #fill in with today's date (year_month_day) and create this folder in Data Files
n_samples = 1000


##set the parameter ranges/distributions
# Matches to Table 2: Variability in all input/parameter values based on field measurements or experiments.

#1
# beta = 50.98      #vary this in the environmental variables section, because it's a time of day measure
# beta_min = 49.04    #use range corresponding to times of observations in data set
# beta_max = 52.92
alpha_si = 0.25
alpha_si_min = 0.9*alpha_si    #unknown, so do +- 10%
alpha_si_max = 1.1*alpha_si

#2
epsilon_a = 0.945
epsilon_a_min = 0.94       #given range
epsilon_a_max = 0.95

#3
A_th = 9.218*10^(-5)
A_th_min = 0.9*A_th
A_th_max = 1.1*A_th

#4
A_h = 2.46*10^(-5)          #this is the HB value
A_h_mean = A_h
A_h_sd = 0.43*10^(-5)

#P = 2   #use a fixed value for now

#5  
alpha_so = 0.5
alpha_so_min = 0.9*alpha_so    #unknown, so do +- 10%
alpha_so_max = 1.1*alpha_so

# #5
# a = 0.25            #include in environmental variables 
# a_min = 0.9*a   #unknown, so do +- 10%
# a_max = 1.1*a

#6
alpha_th = 0.5
alpha_th_min = 0.9*alpha_th    #unknown, so do +- 10%
alpha_th_max = 1.1*alpha_th

#T_aK = 2+273.15   #use a fixed value for now

#T_gK = 2+273.15   #use a fixed value for now

#7
epsilon_e = 0.97
epsilon_e_min = 0.9*epsilon_e    #unknown, so do +- 10%
epsilon_e_max = 1.1*epsilon_e

#8
C_l = 2.429809*10^(-7)
C_l_min = 2.301936*10^(-7)    #+-1sd from Coefficient Estimation.R
C_l_max = 2.564785*10^(-7)

#9
n = 1.975485
n_mean = 1.975485    #from model fit
n_sd = 0.007548

#10
l_th = 0.005467 
l_th_min = 0.0053    #range of workers in Church1960
l_th_max = 0.0058

#11
v = 3.1
v_min = 1    #check these later - flight speeds
v_max = 7
#note: v ranges should be different for resting vs. flying - wind speed vs. flight speed 



#12

i0_resting = 21.117*1.3*(1/60)*(1/60)*0.177 #Kammer1974
i0_resting_min = 0.9*i0_resting    #unknown, so do +- 10% (fill in later if found)
i0_resting_max = 1.1*i0_resting

i0_flying  = 21.117*((150+350)/2)*(1/60)*(1/60)*0.143 #Heinrich1975
i0_flying_Kammer  = 21.117*((166+188)/2)*(1/60)*(1/60)*0.06 #Kammer1974
i0_flying_min  = 21.117*(150)*(1/60)*(1/60)*0.143
i0_flying_max  = 21.117*(350)*(1/60)*(1/60)*0.143
i0_flying_min_Kammer = 21.117*166*(1/60)*(1/60)e*0.06
i0_flying_max_Kammer = 21.117*188*(1/60)*(1/60)*0.06


#13
M_b = 0.149
M_b_min = 0.9*M_b
M_b_max = 1.1*M_b

#14
M_th = 0.057
M_th_min = 0.9*M_th
M_th_max = 1.1*M_th

# #20
# M_h = 0.00997
# M_h_mean = M_h
# M_h_sd = 0.00136

#15
E = 0.63
E_min = 0.6
E_max = 0.7

#delta = 5.31*10^(-13)   #treat as a fixed constant for now

#16
c = 3.3472
c_min = 0.9*c
c_max = 1.1*c

#17
delta_T_h = 3
delta_T_h_min = 2
delta_T_h_max = 4


#T_0 = 39+273.15  #treat as fixed for now

#T_thK = 39+273.15 #treat as fixed for now




## Take the latin hypercube sample - n samples from k parameters
k_length = 17
LHS_sample = randomLHS(n_samples,k_length)

## Use the latin hypercube sample to get the parameter values from the ranges
#Param_values_resting_BB = matrix(NA, nrow=n_samples,ncol=k_length)
Param_values_flying_BB = matrix(NA, nrow=n_samples,ncol=k_length)

for(i in 1:n_samples){ #for each sample

  Param_values_flying_BB[i,1] = qunif(LHS_sample[i,1],min=alpha_si_min,max=alpha_si_max)
  Param_values_flying_BB[i,2] = qunif(LHS_sample[i,2],min=epsilon_a_min,max=epsilon_a_max)
  Param_values_flying_BB[i,3] = qunif(LHS_sample[i,3],min=A_th_min,max=A_th_max)
  Param_values_flying_BB[i,4] = qnorm(LHS_sample[i,4],mean=A_h_mean,sd=A_h_sd)
  Param_values_flying_BB[i,5] = qunif(LHS_sample[i,5],min=alpha_so_min,max=alpha_so_max)
  Param_values_flying_BB[i,6] = qunif(LHS_sample[i,6],min=alpha_th_min,max=alpha_th_max)
  Param_values_flying_BB[i,7] = qunif(LHS_sample[i,7],min=epsilon_e_min,max=epsilon_e_max)
  Param_values_flying_BB[i,8] = qunif(LHS_sample[i,8],min=C_l_min,max=C_l_max)
  Param_values_flying_BB[i,9] = qnorm(LHS_sample[i,9],mean=n_mean,sd=n_sd)
  Param_values_flying_BB[i,10] = qunif(LHS_sample[i,10],min=l_th_min,max=l_th_max)
  Param_values_flying_BB[i,11] = qunif(LHS_sample[i,11],min=v_min,max=v_max)
  Param_values_flying_BB[i,12] = qunif(LHS_sample[i,12],min=i0_flying_min_Kammer,max=i0_flying_max_Kammer)
  Param_values_flying_BB[i,13] = qunif(LHS_sample[i,13],min=M_b_min,max=M_b_max)
  Param_values_flying_BB[i,14] = qunif(LHS_sample[i,14],min=M_th_min,max=M_th_max)
  Param_values_flying_BB[i,15] = qunif(LHS_sample[i,15],min=E_min,max=E_max)
  Param_values_flying_BB[i,16] = qunif(LHS_sample[i,16],min=c_min,max=c_max)
  Param_values_flying_BB[i,17] = qunif(LHS_sample[i,17],min=delta_T_h_min,max=delta_T_h_max)
  
  #note, possibly make so resting and flying samples are the same except for i_0 and v
  
}

#write.csv(Param_values_resting_BB,file="Parameter_Sample_resting_BB_1000.csv")
write.csv(Param_values_flying_BB,file="ParameterSample_FlyingBee_BB_1000.csv")











####################################################################################
################ Now do the environmental variables ################################
ArrianDataSet = read.csv('HiveActivityWeatherdataset.csv',header=TRUE)
n_samples = 1000

#1
P = mean(ArrianDataSet$meansolarstation)      #solar radiations observed by Arrian = a good range of typical Irish weather
P_min = min(ArrianDataSet$meansolarstation)
P_max = max(ArrianDataSet$meansolarstation)

#2
T_aK = mean(ArrianDataSet$thermtemp)+273.15   #0-30C should cover it
T_aK_min = 0+273.15
T_aK_max = 50+273.15

#3
T_gK = 11+273.15        #from https://www.met.ie/forecasts/farming/agricultural-data-report May 25 2021
T_gK_min = 6+273.15     #https://www.farmersjournal.ie/soil-temperature-rise-to-help-growth-151968, https://www.farmersjournal.ie/dairy-management-return-to-winter-weather-612799
T_gK_max = 15+273.15       #http://edepositireland.ie/bitstream/handle/2262/71180/Agromet%20Memo%20No.%203.pdf?sequence=1&isAllowed=y

#4
beta_min = 0        #sunrise... 
beta_max = 60.14    #max reached on solstice (in Dublin) https://www.suncalc.org/#/53.3066,-6.2223,13/2021.06.22/13:35/1/3

#5
a = 0.25
a_min = 0.9*a   #unknown, so do +- 10%
a_max = 1.1*a

#6
kappa = 0.024   #approximation, update later
kappa_min = 0.9*kappa    #unknown, so do +- 10%
kappa_max = 1.1*kappa

#7
nu = 2.791*10^(-7)*285.7708^(0.7355)/1.225
nu_min = 0.99*nu    #can be off by up to 1%
nu_max = 1.01*nu

k_length = 7
LHS_sample = randomLHS(n_samples,k_length)

## Use the latin hypercube sample to get the parameter values from the ranges
Envir_values = matrix(NA, nrow=n_samples,ncol=k_length)

for(i in 1:n_samples){ #for each sample
  Envir_values[i,1] = qunif(LHS_sample[i,1],min=P_min,max=P_max)
  Envir_values[i,2] = qunif(LHS_sample[i,2],min=T_aK_min,max=T_aK_max)
  Envir_values[i,3] = qunif(LHS_sample[i,3],min=T_gK_min,max=T_gK_max)
  Envir_values[i,4] = qunif(LHS_sample[i,4],min=beta_min,max=beta_max)
  Envir_values[i,5] = qunif(LHS_sample[i,5],min=a_min,max=a_max)
  Envir_values[i,6] = qunif(LHS_sample[i,6],min=kappa_min,max=kappa_max)
  Envir_values[i,7] = qunif(LHS_sample[i,7],min=nu_min,max=nu_max)
}

write.csv(Envir_values,file="ParameterSample_Enviro_1000.csv")

