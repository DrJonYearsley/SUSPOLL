library(lhs)
setwd("~/SUSPOLL/Code")


# date = "SensitivityAnalysis"           #fill in with today's date (year_month_day) and create this folder in Data Files
n_samples = 10000


##set the parameter ranges/distributions
# Matches to Table 2: Variability in all input/parameter values based on field measurements or experiments.

#1
delta_T_h = 3
delta_T_h_min = 1.6
delta_T_h_max = 4.2

#metabolic rates

#Heinrich - overall body weight for queens 0.25-0.6g midpoint = 0.4125

# i0_resting = 21.117*1.3*(1/60)*(1/60)*0.177 #Kammer1974 - bdy wt
# i0_resting_min = 21.117*0.3*(1/60)*(1/60)*0.06 #Kammer1974 - thx wt, extracted from fig 1 with imageJ
# i0_resting_max = 21.117*23.4*(1/60)*(1/60)*0.06 #Kammer1974 - thx wt, extracted from fig 1 with imageJ
# 
# i0_flying_Heinrich_1  = 21.117*(76.6)*(1/60)*(1/60)*0.4125 #Heinrich1975 - bdy wt
# i0_flying_Heinrich_2  = 21.117*((150+350)/2)*(1/60)*(1/60)*0.143 #Heinrich1975 - thx wt
# i0_flying_Kammer  = 21.117*((166+188)/2)*(1/60)*(1/60)*0.06 #Kammer1974 - thx wt
# i0_flying_min_Heinrich_1  = 21.117*(47.2)*(1/60)*(1/60)*0.4125 #- bdy wt
# i0_flying_max_Heinrich_1  = 21.117*(106.7)*(1/60)*(1/60)*0.4125 #- bdy wt
# i0_flying_min_Heinrich_2  = 21.117*(150)*(1/60)*(1/60)*0.143 #- thx wt
# i0_flying_max_Heinrich_2  = 21.117*(350)*(1/60)*(1/60)*0.143 #- thx wt
# i0_flying_min_Kammer = 21.117*166*(1/60)*(1/60)*0.06 #- thx wt
# i0_flying_max_Kammer = 21.117*188*(1/60)*(1/60)*0.06 #- thx wt

#2
i0_resting = 0.001349728
i0_resting_min = 0.00011
i0_resting_max = 0.0082356

#3 

#default i0 range
# i0_flying = 0.06229515
# i0_flying_min = 0.0584237
# i0_flying_max = 0.0661666


#low i0 range
# i0_flying = 0.025
# i0_flying_min = 0.001349728
# i0_flying_max = 0.05

#combined i0 range
i0_flying_min = 0.001349728
i0_flying_max = 0.0661666


#4
M_b = 0.149
M_b_min = 0.035
M_b_max = 0.351

#5
E = 0.63
# E_min = 0.6  #default
E_min = 0      #combined
E_max = 0.7

#6
M_th = 0.057
M_th_min = 0.014
M_th_max = 0.132

#7
c = 3.349
c_min = 0.9*c
c_max = 1.1*c

#8
r = 0.0367/9
r_min = 0.002
r_max = 0.008

#9
T_mK = 42+273.15
T_mK_min = 40+273.15
T_mK_max = 44+273.15

#10
alpha_si = 0.25
alpha_si_min = 0.9*alpha_si    #unknown, so do +- 10%
alpha_si_max = 1.1*alpha_si

#11
epsilon_a = 0.935
epsilon_a_min = 0.92       #given range
epsilon_a_max = 0.95

#12
A_th = 9.218*10^(-5)
A_th_min = 8.8247*10^(-5)
A_th_max = 10.5683*10^(-5)

#13
A_h = 2.46*10^(-5)          #this is the HB value
A_h_mean = A_h
A_h_sd = 0.43*10^(-5)

#14
alpha_so = 0.5
alpha_so_min = 0.9*alpha_so    #unknown, so do +- 10%
alpha_so_max = 1.1*alpha_so

#15
alpha_th = 0.5
alpha_th_min = 0.9*alpha_th    #unknown, so do +- 10%
alpha_th_max = 1.1*alpha_th

#16
f = 0.25  #ground reflectance/albedo
f_min = 0.17
f_max = 0.32

#17
P=332.3878
P_min = 29.87
P_max = 1041

#18
T_g = 17.1
T_g_min = 3.1
T_g_max = 19.6

#19
T_air_min = 0
T_air_max = 50

#20
epsilon_e = 0.97
epsilon_e_min = 0.955
epsilon_e_max = 0.99

#21
s = 0.9
s_min = s*0.9
s_max = s*1.1

#22
C_l = 2.429809*10^(-7)
C_l_min = 2.301936*10^(-7)    #+-1sd from Coefficient Estimation.R
C_l_max = 2.564785*10^(-7)

#23
n = 1.975485
n_mean = 1.975485    #from model fit
n_sd = 0.007548

#24
l_th = 0.005467 
l_th_min = 0.0053    #range of workers in Church1960
l_th_max = 0.0058

#25
v = 4.1
v_min = 1    # flight speeds
v_max = 5.5
#note: v ranges should be different for resting vs. flying - wind speed vs. flight speed 

#26
T_0 = 30
T_0_min = 20
T_0_max = 39








## Take the latin hypercube sample - n samples from k parameters
k_length = 26
LHS_sample = randomLHS(n_samples,k_length)

## Use the latin hypercube sample to get the parameter values from the ranges
#Param_values_resting_BB = matrix(NA, nrow=n_samples,ncol=k_length)
Param_values_BB = matrix(NA, nrow=n_samples,ncol=k_length)

for(i in 1:n_samples){ #for each sample

  Param_values_BB[i,1] = qunif(LHS_sample[i,1],min=delta_T_h_min,max=delta_T_h_max)
  Param_values_BB[i,2] = qunif(LHS_sample[i,2],min=i0_resting_min,max=i0_resting_max)
  Param_values_BB[i,3] = qunif(LHS_sample[i,3],min=i0_flying_min,max=i0_flying_max)
  Param_values_BB[i,4] = qunif(LHS_sample[i,4],min=M_b_min,max=M_b_max)
  Param_values_BB[i,6] = qunif(LHS_sample[i,6],min=M_th_min,max=M_th_max)
  Param_values_BB[i,5] = qunif(LHS_sample[i,5],min=E_min,max=E_max)
  Param_values_BB[i,7] = qunif(LHS_sample[i,7],min=c_min,max=c_max)
  Param_values_BB[i,8] = qunif(LHS_sample[i,8],min=r_min,max=r_max)
  Param_values_BB[i,9] = qunif(LHS_sample[i,9],min=T_mK_min,max=T_mK_max)
  Param_values_BB[i,10] = qunif(LHS_sample[i,10],min=alpha_si_min,max=alpha_si_max)
  Param_values_BB[i,11] = qunif(LHS_sample[i,11],min=epsilon_a_min,max=epsilon_a_max)
  Param_values_BB[i,12] = qunif(LHS_sample[i,12],min=A_th_min,max=A_th_max)
  Param_values_BB[i,13] = qnorm(LHS_sample[i,13],mean=A_h_mean,sd=A_h_sd)
  Param_values_BB[i,14] = qunif(LHS_sample[i,14],min=alpha_so_min,max=alpha_so_max)
  Param_values_BB[i,15] = qunif(LHS_sample[i,15],min=alpha_th_min,max=alpha_th_max)
  Param_values_BB[i,16] = qunif(LHS_sample[i,16],min=f_min,max=f_max)
  Param_values_BB[i,17] = qunif(LHS_sample[i,17],min=P_min,max=P_max)
  Param_values_BB[i,18] = qunif(LHS_sample[i,18],min=T_g_min,max=T_g_max)
  Param_values_BB[i,19] = qunif(LHS_sample[i,19],min=T_air_min,max=T_air_max)
  Param_values_BB[i,20] = qunif(LHS_sample[i,20],min=epsilon_e_min,max=epsilon_e_max)
  Param_values_BB[i,21] = qunif(LHS_sample[i,21],min=s_min,max=s_max)
  Param_values_BB[i,22] = qunif(LHS_sample[i,22],min=C_l_min,max=C_l_max)
  Param_values_BB[i,23] = qnorm(LHS_sample[i,23],mean=n_mean,sd=n_sd)
  Param_values_BB[i,24] = qunif(LHS_sample[i,24],min=l_th_min,max=l_th_max)
  Param_values_BB[i,25] = qunif(LHS_sample[i,25],min=v_min,max=v_max)
  Param_values_BB[i,26] = qunif(LHS_sample[i,26],min=T_0_min,max=T_0_max)

  #note, possibly make so resting and flying samples are the same except for i_0 and v
  
}

#write.csv(Param_values_BB,file="Parameter_Sample_BB_10000.csv")
#write.csv(Param_values_BB,file="ParameterSample_BB_10000_lowi0.csv")
write.csv(Param_values_BB,file="ParameterSample_BB_10000_combined.csv")











# ####################################################################################
# ################ Now do the environmental variables ################################
# ArrianDataSet = read.csv('HiveActivityWeatherdataset.csv',header=TRUE)
# n_samples = 1000
# 
# #1
# P = mean(ArrianDataSet$meansolarstation)      #solar radiations observed by Arrian = a good range of typical Irish weather
# P_min = min(ArrianDataSet$meansolarstation)
# P_max = max(ArrianDataSet$meansolarstation)
# 
# #2
# T_aK = mean(ArrianDataSet$thermtemp)+273.15   #0-30C should cover it
# T_aK_min = 0+273.15
# T_aK_max = 50+273.15
# 
# #3
# T_gK = 11+273.15        #from https://www.met.ie/forecasts/farming/agricultural-data-report May 25 2021
# T_gK_min = 6+273.15     #https://www.farmersjournal.ie/soil-temperature-rise-to-help-growth-151968, https://www.farmersjournal.ie/dairy-management-return-to-winter-weather-612799
# T_gK_max = 15+273.15       #http://edepositireland.ie/bitstream/handle/2262/71180/Agromet%20Memo%20No.%203.pdf?sequence=1&isAllowed=y
# 
# #4
# beta_min = 0        #sunrise... 
# beta_max = 60.14    #max reached on solstice (in Dublin) https://www.suncalc.org/#/53.3066,-6.2223,13/2021.06.22/13:35/1/3
# 
# #5
# a = 0.25
# a_min = 0.9*a   #unknown, so do +- 10%
# a_max = 1.1*a
# 
# #6
# kappa = 0.024   #approximation, update later
# kappa_min = 0.9*kappa    #unknown, so do +- 10%
# kappa_max = 1.1*kappa
# 
# #7
# nu = 2.791*10^(-7)*285.7708^(0.7355)/1.225
# nu_min = 0.99*nu    #can be off by up to 1%
# nu_max = 1.01*nu
# 
# k_length = 7
# LHS_sample = randomLHS(n_samples,k_length)
# 
# ## Use the latin hypercube sample to get the parameter values from the ranges
# Envir_values = matrix(NA, nrow=n_samples,ncol=k_length)
# 
# for(i in 1:n_samples){ #for each sample
#   Envir_values[i,1] = qunif(LHS_sample[i,1],min=P_min,max=P_max)
#   Envir_values[i,2] = qunif(LHS_sample[i,2],min=T_aK_min,max=T_aK_max)
#   Envir_values[i,3] = qunif(LHS_sample[i,3],min=T_gK_min,max=T_gK_max)
#   Envir_values[i,4] = qunif(LHS_sample[i,4],min=beta_min,max=beta_max)
#   Envir_values[i,5] = qunif(LHS_sample[i,5],min=a_min,max=a_max)
#   Envir_values[i,6] = qunif(LHS_sample[i,6],min=kappa_min,max=kappa_max)
#   Envir_values[i,7] = qunif(LHS_sample[i,7],min=nu_min,max=nu_max)
# }
# 
# write.csv(Envir_values,file="ParameterSample_Enviro_1000.csv")
# 
