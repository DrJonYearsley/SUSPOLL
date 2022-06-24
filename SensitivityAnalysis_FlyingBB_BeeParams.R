# Take the data generated in HeatExchangeODE_Variability_flying.m or 
# HeatExchangeODE_Variability_resting.m and do the regression as described in 
# Saltelli2006

setwd("~/SUSPOLL/Code")
# install.packages("dismo")
# install.packages("gbm")
library(dismo)
library(gbm)

interaction_level = 2   #level of interactions to consider in the BRT fitting (tree complexity)
n_params = 25 #25 because only using one i0 param

### Read in the data for parameter samples and hear flux ###

#Parameter_Sample_resting = read.csv(file="Parameter_Sample_resting_1000.csv")
Parameter_Values = read.csv(file="ParameterSample_BB_1000.csv")
Parameter_Values = read.csv(file="ParameterSample_BB_2000_lowi0.csv")
Parameter_Values_flying = Parameter_Values[,c(-1,-3)]  #remove the resting i0
Parameter_Values_shivering = Parameter_Values[,c(-1,-3)]  #remove the resting i0
Parameter_Values_resting = Parameter_Values[,c(-1,-4)]  #remove the flying i0
Parameter_Values_resting[ ,24]  = 0  #no flight
Parameter_Values_shivering[ ,24] = 0    #no flight
colnames(Parameter_Values_resting) = c('delta_T_h','i0','M_b','E','M_th',
                                      'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                      'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                      'f_s','c_l','n','l_th','v','T_0')
colnames(Parameter_Values_flying) = c('delta_T_h','i0','M_b','E','M_th',
                                      'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                      'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                      'f_s','c_l','n','l_th','v','T_0')
colnames(Parameter_Values_shivering) = c('delta_T_h','i0','M_b','E','M_th',
                                      'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                      'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                      'f_s','c_l','n','l_th','v','T_0')



Equilibria_flying = read.csv(file="BB_Thorax_Equilibria_Variability_flying_case2.csv",header=FALSE)
Equilibria_resting = read.csv(file="BB_Thorax_Equilibria_Variability_resting_case2.csv",header=FALSE)
Equilibria_shivering = read.csv(file="BB_Thorax_Equilibria_Variability_shivering_case2.csv",header=FALSE)
#T_eq_default = 297.18-273.15   #equlibrium T_th from default params (not updated)
#Q_default = -0.0687


#hist(Equilibria_flying[,1],main='Equilibrium thorax temp of flying bumblebee',xlab='T_th')
#points(T_eq_default,0,pch=20,col='red',cex=2)
#hist(Equilibria_resting[,1],main='Equilibrium thorax temp of flying bumblebee',xlab='T_th')




reg_data_equilibria_flying = cbind(Equilibria_flying,Parameter_Values_flying)
reg_data_equilibria_resting = cbind(Equilibria_resting,Parameter_Values_resting)
reg_data_equilibria_shivering = cbind(Equilibria_resting,Parameter_Values_shivering)
colnames(reg_data_equilibria_flying) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                         'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                         'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                         'f_s','c_l','n','l_th','v','T_0')
colnames(reg_data_equilibria_resting) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                         'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                         'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                         'f_s','c_l','n','l_th','v','T_0')
colnames(reg_data_equilibria_shivering) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                         'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                         'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                         'f_s','c_l','n','l_th','v','T_0')
reg_data_equilibria_flying_edited = na.omit(reg_data_equilibria_flying)
reg_data_equilibria_resting_edited = na.omit(reg_data_equilibria_resting)
reg_data_equilibria_shivering_edited = na.omit(reg_data_equilibria_shivering)

par(mfrow=c(3,1))
hist(reg_data_equilibria_flying_edited[ ,1], main='equilibrium thorax temps while flying',xlab='temperature')
hist(reg_data_equilibria_shivering_edited[ ,1], main='equilibrium thorax temps while shivering',xlab='temperature')
hist(reg_data_equilibria_resting_edited[ ,1], main='equilibrium thorax temps while resting',xlab='temperature')

maxouts_flying=which(reg_data_equilibria_flying_edited[ ,1]==100)
maxouts_resting=which(reg_data_equilibria_resting_edited[ ,1]==100)
maxouts_shivering=which(reg_data_equilibria_shivering_edited[ ,1]==100)

if(length(maxouts_flying)>0){
  reg_data_equilibria_flying_edited = reg_data_equilibria_flying_edited[-maxouts_flying, ]  #remove the rows with temp 100
}
if(length(maxouts_resting)>0){
  reg_data_equilibria_resting_edited = reg_data_equilibria_resting_edited[-maxouts_resting, ]  #remove the rows with temp 100
}
if(length(maxouts_shivering)>0){
  reg_data_equilibria_shivering_edited = reg_data_equilibria_shivering_edited[-maxouts_shivering, ]  #remove the rows with temp 100
}
par(mfrow=c(3,1))
hist(reg_data_equilibria_resting_edited[ ,1], main='equilibrium thorax temps while resting',xlab='temperature')
hist(reg_data_equilibria_shivering_edited[ ,1], main='equilibrium thorax temps while shivering',xlab='temperature')
hist(reg_data_equilibria_flying_edited[ ,1], main='equilibrium thorax temps while flying',xlab='temperature')

# ### fit the least squares regression model 
# ls_mod_heatflux = lm(heatflux~.,data=reg_data_heatflux)
# summary(ls_mod_heatflux)
# 
# ls_mod_equilibria = lm(equilibria~.,data=reg_data_equilibria)
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

BRT_equilibria_flying = gbm.step(data=reg_data_equilibria_flying_edited,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
par(mfrow=c(4,4))
gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
flying_contributions = BRT_equilibria_flying$contributions      #the influence of each parameter
flying_R2 = BRT_equilibria_flying$cv.statistics$correlation.mean^2 #might be the R^2 value... 

BRT_equilibria_resting = gbm.step(data=reg_data_equilibria_resting_edited,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
resting_contributions = BRT_equilibria_resting$contributions      #the influence of each parameter
resting_R2 = BRT_equilibria_resting$cv.statistics$correlation.mean^2 #might be the R^2 value... 

BRT_equilibria_shivering = gbm.step(data=reg_data_equilibria_shivering_edited,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
shivering_contributions = BRT_equilibria_shivering$contributions      #the influence of each parameter
shivering_R2 = BRT_equilibria_shivering$cv.statistics$correlation.mean^2 #might be the R^2 value... 

Interactions_equilibria_flying = gbm.interactions(BRT_equilibria_flying)
Interactions_equilibria_flying$interactions
Interactions_equilibria_flying$rank.list





######################################################################################################
########### Check for sufficient sample size ###################################
