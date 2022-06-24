##Example of code written to efficiently use repeated chunks 
##Previously I had one version of this code for honeybees and one for bumblebees
##Then I re-wrote it to have a single version for either/or
##Which has a few parts that are specific to the bee type 
##And a few parts that are the same no matter what bee type



library(lhs)   #package used to make the latin hypercube sample
setwd("~/SUSPOLL/Code")  #set the directory where the file will be saved


n_samples = 10000
k_length = 30
Bumblebee = TRUE    #only one of Bumblebee and Honeybee can be true
Honeybee = FALSE

### Storage vector to be filled in 
Param_values = matrix(NA, nrow=n_samples,ncol=k_length)

## Take the latin hypercube sample - n samples from k parameters
LHS_sample = randomLHS(n_samples,k_length)



##set the parameter ranges/distributions

####### these ones are for both HB and BB
#5
E_min = 0      
E_max = 0.7

#7
c_min = 0.9*c
c_max = 1.1*c

#10
alpha_si_min = 0.9*alpha_si   
alpha_si_max = 1.1*alpha_si

#14
alpha_so_min = 0.9*alpha_so    
alpha_so_max = 1.1*alpha_so

#15
alpha_th_min = 0.9*alpha_th    #
alpha_th_max = 1.1*alpha_th

#16
f_min = 0.17
f_max = 0.32

#17
P_min = 29.87
P_max = 1041

#18
T_g_min = 3.1
T_g_max = 19.6

#19
T_air_min = 0
T_air_max = 50



####### these ones are BB only, but included for both as required in matlab code
#8
r_min = 0.002
r_max = 0.008

####### these ones are HB only, but included for both as required in matlab code
#27 
R_0_min = 0.00025  
R_0_max = 0.002  

#28
D_A_min = 0.5*D_A
D_A_max = 3*D_A

#29
h_fg_min = 0.9*h_fg
h_fg_max = 1.1*h_fg

#30
rh_min = 0.3920
rh_max = 0.9349


if(Honeybe==TRUE){  ###these are the HB specific values
  #1
  delta_T_h_min = 1.6
  delta_T_h_max = 4.2
  
  #2
  i0_resting_min = 0.9*i0_resting
  i0_resting_max = 1.1*i0_resting
  
  #3 
  i0_flying_min = 4.52*10^(-4)
  i0_flying_max = 1.5*i0_flying
  
  #4
  M_b_min = 0.100-0.0202
  M_b_max = 0.1202
  
  #6
  M_th_min = M_th-0.0029
  M_th_max = M_th+0.0029
  
  #9
  T_mK_min = 42+273.15
  T_mK_max = 46+273.15
  
  #11
  epsilon_a_min = 0.90       
  epsilon_a_max = 0.92
  
  #12
  A_th_min = A_th-0.29*10^(-5)
  A_th_max = A_th+0.29*10^(-5)
  
  #13
  A_h_mean = A_h
  A_h_sd = 0.43*10^(-5)
  
  #20
  epsilon_e_min = 0.955
  epsilon_e_max = 0.99
  
  #21
  s_min = 0.99   
  s_max = 1  
  
  #22
  C_l_min = 2.301936*10^(-7)    
  C_l_max = 2.564785*10^(-7)
  
  #23
  n_mean = 1.975485    
  n_sd = 0.007548
  
  #24
  l_th_min = 0.0053    
  l_th_max = 0.0058
  
  #25
  v_min = 1    
  v_max = 5.5

  #26
  T_0 = 39
  T_0_min = 30
  T_0_max = 45
}

if(Bumblebee==TRUE){  ###these are the BB-specific values
  #1
  delta_T_h_min = 1.6
  delta_T_h_max = 4.2
  
  #2
  i0_resting_min = 0.00011
  i0_resting_max = 0.0082356
  
  #3 
  i0_flying_min = 0.001349728
  i0_flying_max = 0.0661666
  
  #4
  M_b_min = 0.035
  M_b_max = 0.351
  
  #6
  M_th_min = 0.014
  M_th_max = 0.132
  
  #8
  r_min = 0.002
  r_max = 0.008
  
  #9
  T_mK_min = 40+273.15
  T_mK_max = 44+273.15
  
  #11
  epsilon_a_min = 0.92       
  epsilon_a_max = 0.95
  
  #12
  A_th_min = 8.8247*10^(-5)
  A_th_max = 10.5683*10^(-5)
  
  #13
  A_h_mean = A_h
  A_h_sd = 0.43*10^(-5)
  
  #20
  epsilon_e_min = 0.955
  epsilon_e_max = 0.99
  
  #21
  s_min = 0.99   
  s_max = 1    
  
  #22
  C_l_min = 2.301936*10^(-7)    
  C_l_max = 2.564785*10^(-7)
  
  #23
  n_mean = 1.975485    
  n_sd = 0.007548
  
  #24
  l_th_min = 0.0053    
  l_th_max = 0.0058
  
  #25
  v_min = 1    
  v_max = 5.5

  #26
  T_0_min = 20
  T_0_max = 39
}



## Use the latin hypercube sample to get the parameter values from the ranges
for(i in 1:n_samples){ #for each sample
  
  Param_values[i,1] = qunif(LHS_sample[i,1],min=delta_T_h_min,max=delta_T_h_max)
  Param_values[i,2] = qunif(LHS_sample[i,2],min=i0_resting_min,max=i0_resting_max)
  Param_values[i,3] = qunif(LHS_sample[i,3],min=i0_flying_min,max=i0_flying_max)
  Param_values[i,4] = qunif(LHS_sample[i,4],min=M_b_min,max=M_b_max)
  Param_values[i,6] = qunif(LHS_sample[i,6],min=M_th_min,max=M_th_max)
  Param_values[i,5] = qunif(LHS_sample[i,5],min=E_min,max=E_max)
  Param_values[i,7] = qunif(LHS_sample[i,7],min=c_min,max=c_max)
  Param_values[i,8] = qunif(LHS_sample[i,8],min=r_min,max=r_max)
  Param_values[i,9] = qunif(LHS_sample[i,9],min=T_mK_min,max=T_mK_max)
  Param_values[i,10] = qunif(LHS_sample[i,10],min=alpha_si_min,max=alpha_si_max)
  Param_values[i,11] = qunif(LHS_sample[i,11],min=epsilon_a_min,max=epsilon_a_max)
  Param_values[i,12] = qunif(LHS_sample[i,12],min=A_th_min,max=A_th_max)
  Param_values[i,13] = qnorm(LHS_sample[i,13],mean=A_h_mean,sd=A_h_sd)
  Param_values[i,14] = qunif(LHS_sample[i,14],min=alpha_so_min,max=alpha_so_max)
  Param_values[i,15] = qunif(LHS_sample[i,15],min=alpha_th_min,max=alpha_th_max)
  Param_values[i,16] = qunif(LHS_sample[i,16],min=f_min,max=f_max)
  Param_values[i,17] = qunif(LHS_sample[i,17],min=P_min,max=P_max)
  Param_values[i,18] = qunif(LHS_sample[i,18],min=T_g_min,max=T_g_max)
  Param_values[i,19] = qunif(LHS_sample[i,19],min=T_air_min,max=T_air_max)
  Param_values[i,20] = qunif(LHS_sample[i,20],min=epsilon_e_min,max=epsilon_e_max)
  Param_values[i,21] = qunif(LHS_sample[i,21],min=s_min,max=s_max)
  Param_values[i,22] = qunif(LHS_sample[i,22],min=C_l_min,max=C_l_max)
  Param_values[i,23] = qnorm(LHS_sample[i,23],mean=n_mean,sd=n_sd)
  Param_values[i,24] = qunif(LHS_sample[i,24],min=l_th_min,max=l_th_max)
  Param_values[i,25] = qunif(LHS_sample[i,25],min=v_min,max=v_max)
  Param_values[i,26] = qunif(LHS_sample[i,26],min=T_0_min,max=T_0_max)
  Param_values[i,27] = qunif(LHS_sample[i,27],min=R_0_min,max=R_0_max)
  Param_values[i,28] = qunif(LHS_sample[i,28],min=D_A_min,max=D_A_max)
  Param_values[i,29] = qunif(LHS_sample[i,29],min=h_fg_min,max=h_fg_max)
  Param_values[i,30] = qunif(LHS_sample[i,30],min=rh_min,max=rh_max)
}

if(Honeybee==TRUE){
  write.csv(Param_values,file="ParameterSample_10000_combined_HB.csv")
}

if(Bumblebee==TRUE){
  write.csv(Param_values,file="ParameterSample_10000_combined_BB.csv")
}









