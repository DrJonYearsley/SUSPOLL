# Take the data generated in HeatExchangeODE_Variability_flying.m or 
# HeatExchangeODE_Variability_resting.m and do the regression as described in 
# Saltelli2006

setwd("~/SUSPOLL/Code")
# install.packages("dismo")
# install.packages("gbm")
library(dismo)
library(gbm)

interaction_level = 2   #level of interactions to consider in the BRT fitting (tree complexity)
n_params_resting = 24 #24 because only using one i0 param and not v
n_params_shivering = 24 #24 because only using one i0 param and not v
n_params_flying = 25 #25 because only using one i0 param

### Read in the data for parameter samples and hear flux ###

Parameter_Values = read.csv(file="Parameter_Sample_BB_10000.csv")
Parameter_Values_resting = Parameter_Values[,c(-1,-4,-25)]  #remove the flying i0 and v
Parameter_Values_shivering = Parameter_Values[,c(-1,-3,-25)]  #remove the resting i0 and v
Parameter_Values_flying = Parameter_Values[,c(-1,-3)]  #remove the resting i0
colnames(Parameter_Values_resting) = c('delta_T_h','i0','M_b','E','M_th',
                                       'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                       'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                       'f_s','c_l','n','l_th','T_0')
colnames(Parameter_Values_shivering) = c('delta_T_h','i0','M_b','E','M_th',
                                         'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                         'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                         'f_s','c_l','n','l_th','T_0')
colnames(Parameter_Values_flying) = c('delta_T_h','i0','M_b','E','M_th',
                                      'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                      'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                      'f_s','c_l','n','l_th','v','T_0')



Equilibria_resting = read.csv(file="BB_Thorax_Equilibria_Variability_resting_case1_10000.csv",header=FALSE)
Equilibria_shivering = read.csv(file="BB_Thorax_Equilibria_Variability_shivering_case1_10000.csv",header=FALSE)
Equilibria_flying = read.csv(file="BB_Thorax_Equilibria_Variability_flying_case1_10000.csv",header=FALSE)


reg_data_equilibria_resting = cbind(Equilibria_resting,Parameter_Values_resting)
reg_data_equilibria_shivering = cbind(Equilibria_resting,Parameter_Values_shivering)
reg_data_equilibria_flying = cbind(Equilibria_flying,Parameter_Values_flying)
colnames(reg_data_equilibria_resting) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                          'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                          'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                          'f_s','c_l','n','l_th','T_0')
colnames(reg_data_equilibria_shivering) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                            'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                            'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                            'f_s','c_l','n','l_th','T_0')
colnames(reg_data_equilibria_flying) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
                                         'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
                                         'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
                                         'f_s','c_l','n','l_th','v','T_0')
reg_data_equilibria_resting_edited = na.omit(reg_data_equilibria_resting)
reg_data_equilibria_shivering_edited = na.omit(reg_data_equilibria_shivering)
reg_data_equilibria_flying_edited = na.omit(reg_data_equilibria_flying)

par(mfrow=c(3,1))
# hist(reg_data_equilibria_resting_edited[ ,1], main='equilibrium thorax temps while resting',xlab='temperature')
# hist(reg_data_equilibria_shivering_edited[ ,1], main='equilibrium thorax temps while shivering',xlab='temperature')
# hist(reg_data_equilibria_flying_edited[ ,1], main='equilibrium thorax temps while flying',xlab='temperature')

maxouts_resting=which(reg_data_equilibria_resting_edited[ ,1]==100)
maxouts_shivering=which(reg_data_equilibria_shivering_edited[ ,1]==100)
maxouts_flying=which(reg_data_equilibria_flying_edited[ ,1]==100)

if(length(maxouts_resting)>0){
  reg_data_equilibria_resting_edited = reg_data_equilibria_resting_edited[-maxouts_resting, ]  #remove the rows with temp 100
}
if(length(maxouts_shivering)>0){
  reg_data_equilibria_shivering_edited = reg_data_equilibria_shivering_edited[-maxouts_shivering, ]  #remove the rows with temp 100
}
if(length(maxouts_flying)>0){
  reg_data_equilibria_flying_edited = reg_data_equilibria_flying_edited[-maxouts_flying, ]  #remove the rows with temp 100
}
par(mfrow=c(3,1))
hist(reg_data_equilibria_resting_edited[ ,1], main='equilibrium thorax temps while resting',xlab='temperature')
hist(reg_data_equilibria_shivering_edited[ ,1], main='equilibrium thorax temps while shivering',xlab='temperature')
hist(reg_data_equilibria_flying_edited[ ,1], main='equilibrium thorax temps while flying',xlab='temperature')

     
     
     
######################## Fitting #################################

BRT_equilibria_resting = gbm.step(data=reg_data_equilibria_resting_edited,gbm.x=2:(n_params_resting+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
resting_contributions = BRT_equilibria_resting$contributions      #the influence of each parameter
resting_R2 = BRT_equilibria_resting$cv.statistics$correlation.mean^2 #might be the R^2 value... 

BRT_equilibria_shivering = gbm.step(data=reg_data_equilibria_shivering_edited,gbm.x=2:(n_params_shivering+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
shivering_contributions = BRT_equilibria_shivering$contributions      #the influence of each parameter
shivering_R2 = BRT_equilibria_shivering$cv.statistics$correlation.mean^2 #might be the R^2 value... 

BRT_equilibria_flying = gbm.step(data=reg_data_equilibria_flying_edited,gbm.x=2:(n_params_flying+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
# par(mfrow=c(4,4))
# gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
flying_contributions = BRT_equilibria_flying$contributions      #the influence of each parameter
flying_R2 = BRT_equilibria_flying$cv.statistics$correlation.mean^2 #might be the R^2 value... 

resting_contributions
shivering_contributions
flying_contributions

write.csv(resting_contributions,file="case1_contributions_resting.csv")
write.csv(shivering_contributions,file="case1_contributions_shivering.csv")
write.csv(flying_contributions,file="case1_contributions_flying.csv")

resting_R2
shivering_R2
flying_R2

write.csv(c(resting_R2,shivering_R2,flying_R2),file="case1_R2_values.csv")



# Interactions_equilibria_resting = gbm.interactions(BRT_equilibria_resting)
# interactions_resting = Interactions_equilibria_resting$interactions
# InteractionsRank_resting = Interactions_equilibria_resting$rank.list
# 
# Interactions_equilibria_shivering = gbm.interactions(BRT_equilibria_shivering)
# Interactions_shivering = Interactions_equilibria_shivering$interactions
# InteractionsRank_shivering = Interactions_equilibria_shivering$rank.list
# 
# Interactions_equilibria_flying = gbm.interactions(BRT_equilibria_flying)
# Interactions_flying = Interactions_equilibria_flying$interactions
# InteractionsRank_flying = Interactions_equilibria_flying$rank.list





######################################################################################################
########### Check for sufficient sample size ###################################

### Check the sampling efficiency
resting_length = dim(reg_data_equilibria_resting_edited)[1]
shivering_length = dim(reg_data_equilibria_shivering_edited)[1]
flying_length = dim(reg_data_equilibria_flying_edited)[1]

sample_1000_resting = sample(1:resting_length,size=1000,replace=FALSE)
sample_1000_shivering = sample(1:shivering_length,size=1000,replace=FALSE)
sample_100_flying = sample(1:flying_length,size=100,replace=FALSE)
sample_2500_resting = sample(1:resting_length,size=2500,replace=FALSE)
sample_2500_shivering = sample(1:shivering_length,size=2500,replace=FALSE)
sample_250_flying = sample(1:flying_length,size=250,replace=FALSE)
sample_5000_resting = sample(1:resting_length,size=5000,replace=FALSE)
sample_5000_shivering = sample(1:shivering_length,size=5000,replace=FALSE)
sample_500_flying = sample(1:flying_length,size=500,replace=FALSE)
sample_7500_resting = sample(1:resting_length,size=7500,replace=FALSE)
sample_7500_shivering = sample(1:shivering_length,size=7500,replace=FALSE)
sample_750_flying = sample(1:flying_length,size=750,replace=FALSE)

reg_data_equilibria_resting_edited_1000 = reg_data_equilibria_resting_edited[sample_1000_resting,]
BRT_resting_1000 = gbm.step(data=reg_data_equilibria_resting_edited_1000,gbm.x=2:(n_params_resting+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_resting_1000 = BRT_resting_1000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_shivering_edited_1000 = reg_data_equilibria_shivering_edited[sample_1000_shivering,]
BRT_shivering_1000 = gbm.step(data=reg_data_equilibria_shivering_edited_1000,gbm.x=2:(n_params_shivering+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_shivering_1000 = BRT_shivering_1000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_flying_edited_100 = reg_data_equilibria_flying_edited[sample_100_flying,]
BRT_flying_100 = gbm.step(data=reg_data_equilibria_flying_edited_100,gbm.x=2:(n_params_flying+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_flying_100 = BRT_flying_100$contributions[,2]      #the influence of each parameter

write.csv(Influences_resting_1000,file="case1_Influences_resting_1000.csv")
write.csv(Influences_shivering_1000,file="case1_Influences_shivering_1000.csv")
write.csv(Influences_flying_100,file="case1_Influences_flying_100.csv")




reg_data_equilibria_resting_edited_2500 = reg_data_equilibria_resting_edited[sample_2500_resting,]
BRT_resting_2500 = gbm.step(data=reg_data_equilibria_resting_edited_2500,gbm.x=2:(n_params_resting+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_resting_2500 = BRT_resting_2500$contributions[,2]      #the influence of each parameter

reg_data_equilibria_shivering_edited_2500 = reg_data_equilibria_shivering_edited[sample_2500_shivering,]
BRT_shivering_2500 = gbm.step(data=reg_data_equilibria_shivering_edited_2500,gbm.x=2:(n_params_shivering+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_shivering_2500 = BRT_shivering_2500$contributions[,2]      #the influence of each parameter

reg_data_equilibria_flying_edited_250 = reg_data_equilibria_flying_edited[sample_250_flying,]
BRT_flying_250 = gbm.step(data=reg_data_equilibria_flying_edited_250,gbm.x=2:(n_params_flying+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_flying_250 = BRT_flying_250$contributions[,2]      #the influence of each parameter

write.csv(Influences_resting_2500,file="case1_Influences_resting_2500.csv")
write.csv(Influences_shivering_2500,file="case1_Influences_shivering_2500.csv")
write.csv(Influences_flying_250,file="case1_Influences_flying_250.csv")



reg_data_equilibria_resting_edited_5000 = reg_data_equilibria_resting_edited[sample_5000_resting,]
BRT_resting_5000 = gbm.step(data=reg_data_equilibria_resting_edited_5000,gbm.x=2:(n_params_resting+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_resting_5000 = BRT_resting_5000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_shivering_edited_5000 = reg_data_equilibria_shivering_edited[sample_5000_shivering,]
BRT_shivering_5000 = gbm.step(data=reg_data_equilibria_shivering_edited_5000,gbm.x=2:(n_params_shivering+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_shivering_5000 = BRT_shivering_5000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_flying_edited_500 = reg_data_equilibria_flying_edited[sample_500_flying,]
BRT_flying_500 = gbm.step(data=reg_data_equilibria_flying_edited_500,gbm.x=2:(n_params_flying+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_flying_500 = BRT_flying_500$contributions[,2]      #the influence of each parameter

write.csv(Influences_resting_5000,file="case1_Influences_resting_5000.csv")
write.csv(Influences_shivering_5000,file="case1_Influences_shivering_5000.csv")
write.csv(Influences_flying_500,file="case1_Influences_flying_500.csv")



reg_data_equilibria_resting_edited_7500 = reg_data_equilibria_resting_edited[sample_7500_resting,]
BRT_resting_7500 = gbm.step(data=reg_data_equilibria_resting_edited_7500,gbm.x=2:(n_params_resting+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_resting_7500 = BRT_resting_7500$contributions[,2]      #the influence of each parameter

reg_data_equilibria_shivering_edited_7500 = reg_data_equilibria_shivering_edited[sample_7500_shivering,]
BRT_shivering_7500 = gbm.step(data=reg_data_equilibria_shivering_edited_7500,gbm.x=2:(n_params_shivering+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_shivering_7500 = BRT_shivering_7500$contributions[,2]      #the influence of each parameter

reg_data_equilibria_flying_edited_750 = reg_data_equilibria_flying_edited[sample_750_flying,]
BRT_flying_750 = gbm.step(data=reg_data_equilibria_flying_edited_750,gbm.x=2:(n_params_flying+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_flying_750 = BRT_flying_750$contributions[,2]      #the influence of each parameter

write.csv(Influences_resting_7500,file="case1_Influences_resting_7500.csv")
write.csv(Influences_shivering_7500,file="case1_Influences_shivering_7500.csv")
write.csv(Influences_flying_750,file="case1_Influences_flying_750.csv")



Influences_resting_10000 = BRT_equilibria_resting$contributions[,2]
Influences_shivering_10000 = BRT_equilibria_shivering$contributions[,2]
Influences_flying_2000 = BRT_equilibria_flying$contributions[,2]

write.csv(Influences_resting_10000,file="case1_Influences_resting_10000.csv")
write.csv(Influences_shivering_10000,file="case1_Influences_shivering_10000.csv")
write.csv(Influences_flying_2000,file="case1_Influences_flying_2000.csv")





# For 100-250
Average_Influences_resting_1000_2500 = apply(cbind(Influences_resting_1000,Influences_resting_2500),1,mean)
Inner_sum_1000 = sum(Influences_resting_1000*log(Influences_resting_1000)/2,na.rm=TRUE)
Inner_sum_2500 = sum(Influences_resting_2500*log(Influences_resting_2500)/2,na.rm=TRUE)
D_resting_1000_2500 = Inner_sum_1000+Inner_sum_2500 - sum(Average_Influences_resting_1000_2500*log(Average_Influences_resting_1000_2500),na.rm=TRUE)

Average_Influences_shivering_1000_2500 = apply(cbind(Influences_shivering_1000,Influences_shivering_2500),1,mean)
Inner_sum_1000 = sum(Influences_shivering_1000*log(Influences_shivering_1000)/2,na.rm=TRUE)
Inner_sum_2500 = sum(Influences_shivering_2500*log(Influences_shivering_2500)/2,na.rm=TRUE)
D_shivering_1000_2500 = Inner_sum_1000+Inner_sum_2500 - sum(Average_Influences_shivering_1000_2500*log(Average_Influences_shivering_1000_2500),na.rm=TRUE)

Average_Influences_flying_100_250 = apply(cbind(Influences_flying_100,Influences_flying_250),1,mean)
Inner_sum_100 = sum(Influences_flying_100*log(Influences_flying_100)/2,na.rm=TRUE)
Inner_sum_250 = sum(Influences_flying_250*log(Influences_flying_250)/2,na.rm=TRUE)
D_flying_100_250 = Inner_sum_100+Inner_sum_250 - sum(Average_Influences_flying_100_250*log(Average_Influences_flying_100_250),na.rm=TRUE)

# For 250-500
Average_Influences_resting_2500_5000 = apply(cbind(Influences_resting_2500,Influences_resting_5000),1,mean)
Inner_sum_2500 = sum(Influences_resting_2500*log(Influences_resting_2500)/2,na.rm=TRUE)
Inner_sum_5000 = sum(Influences_resting_5000*log(Influences_resting_5000)/2,na.rm=TRUE)
D_resting_2500_5000 = Inner_sum_2500+Inner_sum_5000 - sum(Average_Influences_resting_2500_5000*log(Average_Influences_resting_2500_5000),na.rm=TRUE)

Average_Influences_shivering_2500_5000 = apply(cbind(Influences_shivering_2500,Influences_shivering_5000),1,mean)
Inner_sum_2500 = sum(Influences_shivering_2500*log(Influences_shivering_2500)/2,na.rm=TRUE)
Inner_sum_5000 = sum(Influences_shivering_5000*log(Influences_shivering_5000)/2,na.rm=TRUE)
D_shivering_2500_5000 = Inner_sum_2500+Inner_sum_5000 - sum(Average_Influences_shivering_2500_5000*log(Average_Influences_shivering_2500_5000),na.rm=TRUE)

Average_Influences_flying_250_500 = apply(cbind(Influences_flying_250,Influences_flying_500),1,mean)
Inner_sum_250 = sum(Influences_flying_250*log(Influences_flying_250)/2,na.rm=TRUE)
Inner_sum_500 = sum(Influences_flying_500*log(Influences_flying_500)/2,na.rm=TRUE)
D_flying_250_500 = Inner_sum_250+Inner_sum_500 - sum(Average_Influences_flying_250_500*log(Average_Influences_flying_250_500),na.rm=TRUE)

# For 500-750 
Average_Influences_resting_5000_7500 = apply(cbind(Influences_resting_5000,Influences_resting_7500),1,mean)
Inner_sum_5000 = sum(Influences_resting_5000*log(Influences_resting_5000)/2,na.rm=TRUE)
Inner_sum_7500 = sum(Influences_resting_7500*log(Influences_resting_7500)/2,na.rm=TRUE)
D_resting_5000_7500 = Inner_sum_5000+Inner_sum_7500 - sum(Average_Influences_resting_5000_7500*log(Average_Influences_resting_5000_7500),na.rm=TRUE)

Average_Influences_shivering_5000_7500 = apply(cbind(Influences_shivering_5000,Influences_shivering_7500),1,mean)
Inner_sum_5000 = sum(Influences_shivering_5000*log(Influences_shivering_5000)/2,na.rm=TRUE)
Inner_sum_7500 = sum(Influences_shivering_7500*log(Influences_shivering_7500)/2,na.rm=TRUE)
D_shivering_5000_7500 = Inner_sum_5000+Inner_sum_7500 - sum(Average_Influences_shivering_5000_7500*log(Average_Influences_shivering_5000_7500),na.rm=TRUE)

Average_Influences_flying_500_750 = apply(cbind(Influences_flying_500,Influences_flying_750),1,mean)
Inner_sum_500 = sum(Influences_flying_500*log(Influences_flying_500)/2,na.rm=TRUE)
Inner_sum_750 = sum(Influences_flying_750*log(Influences_flying_750)/2,na.rm=TRUE)
D_flying_500_750 = Inner_sum_500+Inner_sum_750 - sum(Average_Influences_flying_500_750*log(Average_Influences_flying_500_750),na.rm=TRUE)

# For 750-2000
Average_Influences_resting_7500_10000 = apply(cbind(Influences_resting_7500,Influences_resting_10000),1,mean)
Inner_sum_7500 = sum(Influences_resting_7500*log(Influences_resting_7500)/2,na.rm=TRUE)
Inner_sum_10000 = sum(Influences_resting_10000*log(Influences_resting_10000)/2,na.rm=TRUE)
D_resting_7500_10000 = Inner_sum_7500+Inner_sum_10000 - sum(Average_Influences_resting_7500_10000*log(Average_Influences_resting_7500_10000),na.rm=TRUE)

Average_Influences_shivering_7500_10000 = apply(cbind(Influences_shivering_7500,Influences_shivering_10000),1,mean)
Inner_sum_7500 = sum(Influences_shivering_7500*log(Influences_shivering_7500)/2,na.rm=TRUE)
Inner_sum_10000 = sum(Influences_shivering_10000*log(Influences_shivering_10000)/2,na.rm=TRUE)
D_shivering_7500_10000 = Inner_sum_7500+Inner_sum_10000 - sum(Average_Influences_shivering_7500_10000*log(Average_Influences_shivering_7500_10000),na.rm=TRUE)

Average_Influences_flying_750_2000 = apply(cbind(Influences_flying_750,Influences_flying_2000),1,mean)
Inner_sum_750 = sum(Influences_flying_750*log(Influences_flying_750)/2,na.rm=TRUE)
Inner_sum_2000 = sum(Influences_flying_2000*log(Influences_flying_2000)/2,na.rm=TRUE)
D_flying_750_2000 = Inner_sum_750+Inner_sum_2000 - sum(Average_Influences_flying_750_2000*log(Average_Influences_flying_750_2000),na.rm=TRUE)


## All together
D_resting = exp(c(D_resting_1000_2500,D_resting_2500_5000,D_resting_5000_7500,D_resting_7500_10000))
D_shivering = exp(c(D_shivering_1000_2500,D_shivering_2500_5000,D_shivering_5000_7500,D_shivering_7500_10000))
D_flying = exp(c(D_flying_100_250,D_flying_250_500,D_flying_500_750,D_flying_750_2000))

write.csv(cbind(D_resting,D_shivering,D_flying),file="case1_D_values.csv")

par(mfrow=c(3,1))
plot(D_resting,main="Resting Bee Case 1",type='b')
plot(D_shivering, main="Shivering Bee Case 1",type='b')
plot(D_flying, main="Flying Bee Case 1",type='b')

