# Take the data generated in HeatExchangeODE_Variability_resting.m or 
# HeatExchangeODE_Variability_resting.m and do the regression as described in 
# Saltelli2006

setwd("~/SUSPOLL/Code")
# install.packages("dismo")
# install.packages("gbm")
# library(dismo)
# library(gbm)

### Read in the data for parameter samples and hear flux ###

#Parameter_Sample_resting = read.csv(file="Parameter_Sample_resting_1000.csv")
Parameter_Values_resting = read.csv(file="Parameter_Sample_resting_1000.csv")
Parameter_Values_resting = Parameter_Values_resting[,-1]
colnames(Parameter_Values_resting) = c('alpha_si','epsilon_a','A_th','alpha_so','a',
                                      'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                                      'i0_resting','M_b','E','c')

#Heatflux_resting = read.csv(file="heatflux_resting.csv")
Heatflux_resting = read.csv(file="heatflux_resting.csv",header=FALSE)
Equilibria_resting = read.csv(file="equilibria_resting.csv",header=FALSE)
T_eq_default = 293.7593-273.15
Q_default = -0.1039

hist(Heatflux_resting[,1],main='Heat flux of a resting bee',xlab='Q')
points(Q_default,0,pch=20,col='red',cex=2)
hist(Equilibria_resting[,1]-273.15,main='Equilibrium thorax temp of a resting bee',xlab='T_th')
points(T_eq_default,0,pch=20,col='red',cex=2)

### assign each parameter to its name ###
#each column is a particular parameter; fixed values are fixed
alpha_si = Parameter_Values_resting[ ,1] 
epsilon_a = Parameter_Values_resting[ ,2] 
A_th = Parameter_Values_resting[ ,3] 
P = 662.4983
alpha_so = Parameter_Values_resting[ ,4] 
a = Parameter_Values_resting[ ,5]
alpha_th = Parameter_Values_resting[ ,6]
T_aK = 288.5838
T_gK = 282.1500
epsilon_e = Parameter_Values_resting[ ,7] 
C_l = Parameter_Values_resting[ ,8] 
kappa = Parameter_Values_resting[ ,9]
l_th = Parameter_Values_resting[ ,10] 
v = Parameter_Values_resting[ ,11] 
nu = Parameter_Values_resting[ ,12] 
n = Parameter_Values_resting[ ,13] 
i0_resting = Parameter_Values_resting[ ,14] 
M_b = Parameter_Values_resting[ ,15]
E = Parameter_Values_resting[ ,16] 
delta = 5.31*10^(-13)
k = 8.6173*10^(-5)
sigma = 5.67*10^(-8)
c = Parameter_Values_resting[ ,17] 
T_thK = 39+273.15


### Calculte the heat flux??? ###
# heat_flux = (alpha_si*epsilon_a*A_th*P + alpha_so*epsilon_a*A_th*a*P) +  #S
# (alpha_th*epsilon_a*A_th*(delta*T_aK^6+sigma*T_gK^4)) +      #R1
# (epsilon_e*A_th*sigma)*T_thK^4 -                             #R2
# (((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th)*T_thK -              #C1
# (-((C_l*kappa/l_th)*(v*l_th/nu)^n)*A_th*T_aK) +              #C2
# (i0_resting*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)               #I


### combine heat flux/equilibrium temp and parameter values into a data frame ###
reg_data_heatflux = cbind(Heatflux_resting,Parameter_Values_resting[1:16])
colnames(reg_data_heatflux) = c('heatflux','alpha_si','epsilon_a','A_th','alpha_so','a',
                                'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                                'i0_resting','M_b','E')

reg_data_equilibria = cbind(Equilibria_resting,Parameter_Values_resting[1:17])
colnames(reg_data_equilibria) = c('heatflux','alpha_si','epsilon_a','A_th','alpha_so','a',
                                  'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                                  'i0_resting','M_b','E','c')

### fit the least squares regression model 
ls_mod_heatflux = lm(heatflux~.,data=reg_data_heatflux)
summary(ls_mod_heatflux)

ls_mod_equilibria = lm(heatflux~.,data=reg_data_equilibria)
summary(ls_mod_equilibria)


### calculate the standardized regression coefficients
b_heatflux <- summary(ls_mod_heatflux)$coef[-1, 1]
sx_heatflux <- apply(Parameter_Values_resting[1:16],2,sd)
sy_heatflux <- apply(Heatflux_resting,2,sd)
beta_heatflux <- b_heatflux * (sx_heatflux/sy_heatflux)

b_equilibria <- summary(ls_mod_equilibria)$coef[-1, 1]
sx_equilibria <- apply(Parameter_Values_resting[1:17],2,sd)
sy_equilibria <- apply(Equilibria_resting,2,sd)
beta_equilibria <- b_equilibria * (sx_equilibria/sy_equilibria)


### Fit the boosted regression trees 
interaction_level = 1   #1 = no interactions

BRT_heatflux = gbm.step(data=reg_data_heatflux,gbm.x=2:17,gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_heatflux,n.plots=15,write.title=FALSE)    #plots
BRT_heatflux$contributions      #the influence of each parameter
BRT_heatflux$cv.statistics$correlation.mean^2 #might be the R^2 value... 

# Interactions_heatflux = gbm.interactions(BRT_heatflux)
# Interactions_heatflux$interactions
# Interactions_heatflux$rank.list

BRT_equilibria = gbm.step(data=reg_data_equilibria,gbm.x=2:18,gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria,n.plots=15,write.title=FALSE)    #plots
BRT_equilibria$contributions      #the influence of each parameter
BRT_equilibria$cv.statistics$correlation.mean^2 #might be the R^2 value... 

# Interactions_equilibria = gbm.interactions(BRT_equilibria)
# Interactions_equilibria$interactions
# Interactions_equilibria$rank.list
