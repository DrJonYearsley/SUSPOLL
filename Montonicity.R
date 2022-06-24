setwd("~/SUSPOLL/Code")



##set the parameter ranges/distributions
# Matches to Table 2: Variability in all input/parameter values based on field measurements or experiments.

#1
alpha_si = 0.25
alpha_si_min = 0.9*alpha_si    #unknown, so do +- 10%
alpha_si_max = 1.1*alpha_si

#2
epsilon_a = 0.91
epsilon_a_min = 0.90       #given range
epsilon_a_max = 0.92

#3
A_th = 4.5*10^(-5)
A_th_mean = A_th
A_th_sd = 0.43*10^(-5)

#P = 2   #use a fixed value for now

#4  
alpha_so = 0.5
alpha_so_min = 0.9*alpha_so    #unknown, so do +- 10%
alpha_so_max = 1.1*alpha_so

#5
a = 0.25
a_min = 0.9*a   #unknown, so do +- 10%
a_max = 1.1*a

#6
alpha_th = 0.5
alpha_th_min = 0.9*alpha_th    #unknown, so do +- 10%
alpha_th_max = 1.1*alpha_th

T_aK = 288.5838   #use a fixed value for now - mean of one site, one day

T_gK = 282.1500   #use a fixed value for now

#7
epsilon_e = 0.97
epsilon_e_min = 0.9*epsilon_e    #unknown, so do +- 10%
epsilon_e_max = 1.1*epsilon_e

#8
C_l = 0.0749
C_l_min = 0.9*C_l    #unknown, so do +- 10%
C_l_max = 1.1*C_l

#9
kappa = 0.024   #approximation, update later
kappa_min = 0.9*kappa    #unknown, so do +- 10%
kappa_max = 1.1*kappa

#10
l_th = 0.004  
l_th_min = 0.9*l_th    #unknown, so do +- 10%
l_th_max = 1.1*l_th

#11
v = 3.1
v_min = 1    #check these later - flight speeds
v_max = 7

#12
nu = 2.791*10^(-7)*285.7708^(0.7355)/1.225
nu_min = 0.99*nu    #can be off by up to 1%
nu_max = 1.01*nu

#13
n = 0.78
n_min = 0.9*n    #unknown, so do +- 10%
n_max = 1.1*n

#14
k = 8.617333262145*10^(-5)
E = 0.63   #maybe should use the random value of E?
M_b = 0.1008   #maybe should use the random value of M_b?
I_resting = 1.31*10^(-3)     #%From Rothe1989, in mW/g, at T_a = 10C
#i0_resting = exp(log(I_resting) - (3/4)*log(M_b) + E/(k*(10+273.15)))  #fit i_0 to the data
i0_resting = 1.1968*10^9
i0_resting_min = 0.9*i0_resting    #unknown, so do +- 10% (fill in later if found)
i0_resting_max = 1.1*i0_resting

I_flying = 0.4     #%From Nachtigall1989, in W/g
#i0_flying = exp(log(I_flying) - (3/4)*log(M_b) + E/(k*(20+273.15))) #fit i_0 to the data
i0_flying = 1.5146*10^11
i0_flying_min = 0.9*i0_flying    #unknown, so do +- 10% (fill in later if found)
i0_flying_max = 1.1*i0_flying

#15
M_b = 0.1008
M_b_mean = M_b
M_b_sd = 0.0202

#16
E = 0.63
E_min = 0.6
E_max = 0.7

delta = 5.31*10^(-13)   #treat as a fixed constant for now

#17
c = 3.3472
c_min = 0.9*c
c_max = 1.1*c

T_0K = 39+273.15  #treat as fixed for now



### Heat flux as fn of each non-linear parameter ###
heat_flux_l_th = function(l_th){
 heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
  (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
  (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
  (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
  (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
  (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
 return(heatflux)
}
curve(heat_flux_l_th,from=l_th_min,to=l_th_max)


heat_flux_v = function(v){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_v,from=v_min,to=v_max)


heat_flux_nu = function(nu){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_nu,from=nu_min,to=nu_max)


heat_flux_n = function(n){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_n,from=n_min,to=n_max)


heat_flux_M_b = function(M_b){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_M_b,from=M_b-2*M_b_sd,to=M_b+2*M_b_sd)

heat_flux_E = function(E){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_E,from=E_min,to=E_max)


heat_flux_T_aK = function(T_aK){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_T_aK,from=0+273.15,to=30+273.15)

heat_flux_T_gK = function(T_gK){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_T_gK,from=0+273.15,to=30+273.15)

heat_flux_T_thK = function(T_thK){
  heatflux =  (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
    (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) -      #R1
    (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
    (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
    (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
    (i0_flying*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I
  return(heatflux)
}
curve(heat_flux_T_thK,from=0+273.15,to=40+273.15)

