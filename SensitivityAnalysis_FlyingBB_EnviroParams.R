# Take the data generated in HeatExchangeODE_Variability_flying.m or 
# HeatExchangeODE_Variability_resting.m and do the regression as described in 
# Saltelli2006

setwd("~/SUSPOLL/Code")
# install.packages("dismo")
# install.packages("gbm")
library(dismo)
library(gbm)

interaction_level = 2   #level of interactions to consider in the BRT fitting (tree complexity)
n_params = 5

### Read in the data for parameter samples and hear flux ###

#Parameter_Sample_resting = read.csv(file="Parameter_Sample_resting_1000.csv")
Parameter_Values_flying = read.csv(file="ParameterSample_Enviro_1000.csv")
Parameter_Values_flying = Parameter_Values_flying[,-1]
colnames(Parameter_Values_flying) = c('P','T_a','T_g','rh','a')

Equilibria_flying = read.csv(file="Thorax_Equilibria_EnviroVariability_resting_1000_HB.csv",header=FALSE)

hist(Equilibria_flying[,1],main='Equilibrium thorax temp of flying bumblebee',xlab='T_th')



### combine heat flux/equilibrium temp and parameter values into a data frame ###
# reg_data_heatflux = cbind(Heatflux_flying,Parameter_Values_flying[,-17])
# colnames(reg_data_heatflux) = c('heatflux','alpha_si','epsilon_a','A_th','A_h','alpha_so',
#                                 'alpha_th','epsilon_e','C_l','n','l_th','v',
#                                 'i0','M_b','M_th','E','c','delta_T_h')

reg_data_equilibria = cbind(Equilibria_flying,Parameter_Values_flying)
colnames(reg_data_equilibria) = c('equilibria','P','T_a','T_g','rh','a')

# ### fit the least squares regression model 
# ls_mod_heatflux = lm(heatflux~.,data=reg_data_heatflux)
# summary(ls_mod_heatflux)
# 
# ls_mod_equilibria = lm(heatflux~.,data=reg_data_equilibria)
# summary(ls_mod_equilibria)
# 
# 
# ### calculate the standardized regression coefficients
# b_heatflux <- summary(ls_mod_heatflux)$coef[-1, 1]
# sx_heatflux <- apply(Parameter_Values_flying[1:16],2,sd)
# sy_heatflux <- apply(Heatflux_flying,2,sd)
# beta_heatflux <- b_heatflux * (sx_heatflux/sy_heatflux)
# 
# b_equilibria <- summary(ls_mod_equilibria)$coef[-1, 1]
# sx_equilibria <- apply(Parameter_Values_flying[1:17],2,sd)
# sy_equilibria <- apply(Equilibria_flying,2,sd)
# beta_equilibria <- b_equilibria * (sx_equilibria/sy_equilibria)


# ### Fit the boosted regression trees
# BRT_heatflux = gbm.step(data=reg_data_heatflux,gbm.x=2:18,gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
# #par(mfrow=c(4,4))
# #gbm.plot(BRT_heatflux,n.plots=15,write.title=FALSE)    #plots
# BRT_heatflux$contributions      #the influence of each parameter
# BRT_heatflux$cv.statistics$correlation.mean^2 #might be the R^2 value... 
# 
# Interactions_heatflux = gbm.interactions(BRT_heatflux)
# Interactions_heatflux$interactions
# Interactions_heatflux$rank.list

BRT_equilibria = gbm.step(data=reg_data_equilibria,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
gbm.plot(BRT_equilibria,n.plots=16,write.title=FALSE)    #plots
BRT_equilibria$contributions      #the influence of each parameter
BRT_equilibria$cv.statistics$correlation.mean^2 #might be the R^2 value... 

# Interactions_equilibria = gbm.interactions(BRT_equilibria)
# Interactions_equilibria$interactions
# Interactions_equilibria$rank.list





######################################################################################################
########### Check for sufficient sample size ###################################
