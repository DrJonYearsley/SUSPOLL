# Take the data generated in HeatExchangeODE_Variability_flying.m or 
# HeatExchangeODE_Variability_resting.m and do the regression as described in 
# Saltelli2006

setwd("~/SUSPOLL/Code")
# install.packages("dismo")
# install.packages("gbm")
library(dismo)
library(gbm)

interaction_level = 2   #level of interactions to consider in the BRT fitting (tree complexity)

### Read in the data for parameter samples and hear flux ###

#Parameter_Sample_resting = read.csv(file="Parameter_Sample_resting_1000.csv")
Parameter_Values_flying = read.csv(file="Parameter_Sample_flying_1000.csv")
Parameter_Values_flying = Parameter_Values_flying[,-1]
colnames(Parameter_Values_flying) = c('alpha_si','epsilon_a','A_th','alpha_so','a',
                       'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                       'i0_flying','M_b','E','c','A_h','M_th','M_h','delta_T_h')

#Heatflux_resting = read.csv(file="heatflux_resting.csv")
Heatflux_flying = read.csv(file="heatflux_flying.csv",header=FALSE)
Equilibria_flying = read.csv(file="equilibria_flying.csv",header=FALSE)
T_eq_default = 297.18-273.15
Q_default = -0.0687

hist(Heatflux_flying[,1],main='Heat flux of a flying bee',xlab='Q')
points(Q_default,0,pch=20,col='red',cex=2)
hist(Equilibria_flying[,1]-273.15,main='Equilibrium thorax temp of a flying bee',xlab='T_th')
points(T_eq_default,0,pch=20,col='red',cex=2)



### combine heat flux/equilibrium temp and parameter values into a data frame ###
reg_data_heatflux = cbind(Heatflux_flying,Parameter_Values_flying[,-17])
colnames(reg_data_heatflux) = c('heatflux','alpha_si','epsilon_a','A_th','alpha_so','a',
                       'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                       'i0_flying','M_b','E','A_h','M_th','M_h','delta_T_h')

reg_data_equilibria = cbind(Equilibria_flying,Parameter_Values_flying)
colnames(reg_data_equilibria) = c('heatflux','alpha_si','epsilon_a','A_th','alpha_so','a',
                                'alpha_th','epsilon_e','C_l','kappa','l_th','v' ,'nu','n',
                                'i0_flying','M_b','E','c','A_h','M_th','M_h','delta_T_h')

### fit the least squares regression model 
ls_mod_heatflux = lm(heatflux~.,data=reg_data_heatflux)
summary(ls_mod_heatflux)

ls_mod_equilibria = lm(heatflux~.,data=reg_data_equilibria)
summary(ls_mod_equilibria)


### calculate the standardized regression coefficients
b_heatflux <- summary(ls_mod_heatflux)$coef[-1, 1]
sx_heatflux <- apply(Parameter_Values_flying[1:16],2,sd)
sy_heatflux <- apply(Heatflux_flying,2,sd)
beta_heatflux <- b_heatflux * (sx_heatflux/sy_heatflux)

b_equilibria <- summary(ls_mod_equilibria)$coef[-1, 1]
sx_equilibria <- apply(Parameter_Values_flying[1:17],2,sd)
sy_equilibria <- apply(Equilibria_flying,2,sd)
beta_equilibria <- b_equilibria * (sx_equilibria/sy_equilibria)


### Fit the boosted regression trees
BRT_heatflux = gbm.step(data=reg_data_heatflux,gbm.x=2:21,gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_heatflux,n.plots=15,write.title=FALSE)    #plots
BRT_heatflux$contributions      #the influence of each parameter
BRT_heatflux$cv.statistics$correlation.mean^2 #might be the R^2 value... 

Interactions_heatflux = gbm.interactions(BRT_heatflux)
Interactions_heatflux$interactions
Interactions_heatflux$rank.list

BRT_equilibria = gbm.step(data=reg_data_equilibria,gbm.x=2:22,gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria,n.plots=15,write.title=FALSE)    #plots
BRT_equilibria$contributions      #the influence of each parameter
BRT_equilibria$cv.statistics$correlation.mean^2 #might be the R^2 value... 

Interactions_equilibria = gbm.interactions(BRT_equilibria)
Interactions_equilibria$interactions
Interactions_equilibria$rank.list





######################################################################################################
########### Check for sufficient sample size ###################################
