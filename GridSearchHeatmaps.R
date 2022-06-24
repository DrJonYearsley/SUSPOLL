setwd("~/SUSPOLL/Code")
library(fields) #used of image.plot to get legends

Bumblebee=F   #only one can be true
Honeybee=T 


# read in the model results

if(Bumblebee==TRUE){
  # read in the model results
  Fly_P1_P2 = read.csv('Fly_P1_P2_Bumblebee.csv',header=FALSE)
  Cool_P1_P2 = read.csv('Cool_P1_P2_Bumblebee.csv',header=FALSE)
  # target for ability to fly is 10\dC
  # target for cooling is 20-25\dC,
  fly_min = 9
  fly_max = 11
  cool_min = 20
  cool_max = 25
  
  ## #these are the parameter values used in the simulations
  i0_min = 0.001349728  #shouldn't go lower than resting i0...
  i0_max = 0.15  #Heinrich max number (rounded up?)
  i0_step = 0.0005
  i0_axis = seq(from=i0_min,by=i0_step,to=i0_max)
  E_min = 0
  E_max = 0.7
  E_step = 0.002
  E_axis = seq(from=E_min,by=E_step,to=E_max)
  
  shiver_target = 10
  diverge_target = 22.5
  
}

if(Honeybee==TRUE){
  # read in the model results
  Fly_P1_P2 = read.csv('Fly_P1_P2_Honeybee.csv',header=FALSE)
  Cool_P1_P2 = read.csv('Cool_P1_P2_Honeybee.csv',header=FALSE)
  # target for ability to fly is 10\dC
  # target for cooling is 35\dC,
  fly_min = 9
  fly_max = 11
  cool_min = 34
  cool_max = 36
  
  i0_default = 3.20*10^(-2)     #I_flying
  i0_min = 4.52*10^(-4)    #shouldn't go lower than resting P1... 
  # i0_max = 1.5*i0_default  #some maximum
  i0_max = 0.02  #some maximum
  i0_step = 0.00005
  i0_axis = seq(from=i0_min,by=i0_step,to=i0_max)
  
  E_min = 0
  E_max = 0.7
  E_step = 0.002
  E_axis = seq(from=E_min,by=E_step,to=E_max)
  
  shiver_target = 10
  diverge_target = 35
  
}

# 
par(mfrow=c(2,1))
image.plot(i0_axis,E_axis,as.matrix(Fly_P1_P2),main='air temp able to warm to 30C')
image.plot(i0_axis,E_axis,as.matrix(Cool_P1_P2),main='air temp abdomen cooling starts')
points(i0_default,E_default,0.95,col='red',pch=20)



##########################################################################################
############ Calculate the euclidean distance from the target values

shiver_difference = Fly_P1_P2-shiver_target
diverge_difference = Cool_P1_P2-diverge_target

par(mfrow=c(2,1))
image.plot(i0_axis,E_axis,abs(as.matrix(Fly_P1_P2)-shiver_target), 
           main='able to fly (difference from target)')
image.plot(i0_axis,E_axis,abs(as.matrix(Cool_P1_P2)-diverge_target),
           main='cooling starts (difference from target)')

Euclid_dist = sqrt(shiver_difference^2+diverge_difference^2)
par(mfrow=c(1,1))
image.plot(i0_axis,E_axis,as.matrix(Euclid_dist),
           main='Euclidean distance from target')


# par(mfrow=c(3,1))
# candidates0 = which(Euclid_dist==0,arr.ind=TRUE)
# #rows for i0 and colums for E
# i0_values0 = i0_axis[candidates0[,1]]
# f_values0 = E_axis[candidates0[,2]]
# plot(i0_values0,f_values0, main = 'candidates that meet the criteria exactly')
# 
# 
# candidates1 = which(Euclid_dist<=1,arr.ind=TRUE)
# #rows for i0 and colums for E
# i0_values1 = i0_axis[candidates1[,1]]
# f_values1 = E_axis[candidates1[,2]]
# plot(i0_values1,f_values1, main = 'candidates within a distance of 1 from the criteria')
# 
# candidates5 = which(Euclid_dist<=5,arr.ind=TRUE)
# #rows for i0 and colums for E
# i0_values5 = i0_axis[candidates5[,1]]
# f_values5 = E_axis[candidates5[,2]]
# plot(i0_values5,f_values5, main = 'candidates within a distance of 5 from the criteria')
# 
