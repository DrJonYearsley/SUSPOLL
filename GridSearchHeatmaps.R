setwd("~/SUSPOLL/Code")
library(fields) #used of image.plot to get legends

# Diverge = read.csv('Diverge.csv',header=FALSE)
# Diverge_fine = read.csv('Diverge_fine.csv',header=FALSE)
# Diverge_fine2 = read.csv('Diverge_fine2.csv',header=FALSE)
Diverge_i0_f_1 = read.csv('Diverge_i0_f_1.csv',header=FALSE)

# Diverge_dist = read.csv('Diverge_dist.csv',header=FALSE)
# Diverge_dist_fine = read.csv('Diverge_dist_fine.csv',header=FALSE)
# Diverge_dist_fine2 = read.csv('Diverge_dist_fine.csv',header=FALSE)
Diverge_dist_i0_f_1 = read.csv('Diverge_dist_i0_f_1.csv',header=FALSE)

# ShiverTo30 = read.csv('ShiverTo30.csv',header=FALSE)
# ShiverTo30_fine = read.csv('ShiverTo30_fine.csv',header=FALSE)
# ShiverTo30_fine2 = read.csv('ShiverTo30_fine.csv',header=FALSE)
ShiverTo30_i0_f_1 = read.csv('ShiverTo30_i0_f_1.csv',header=FALSE)

# ShiverTo30_dist = read.csv('ShiverTo30_dist.csv',header=FALSE)
# ShiverTo30_dist_fine = read.csv('ShiverTo30_dist_fine.csv',header=FALSE)
# ShiverTo30_dist_fine2 = read.csv('ShiverTo30_dist_fine2.csv',header=FALSE)
ShiverTo30_dist_i0_f_1 = read.csv('ShiverTo30_dist_i0_f_1.csv',header=FALSE)


# i0_min = 0
# i0_max = 0.06229515
# r_min = 0
# r_max = 0.005
i0_min = 0.001349728  #shouldn't go lower than resting i0... 
i0_max = 0.15  #Heinrich max number (rounded up?)
f_min = 0.7
f_max = 1


# i0_axis = seq(from=i0_min,by=0.0002,to=i0_max)
# r_axis = seq(from=r_min,by=0.00002,to=r_max)
# 
# i0_axis_fine = seq(from=0,by=0.00001,to=0.0045);
# r_axis_fine = seq(from=r_min,by=0.00002,to=r_max);
# 
# i0_axis_fine2 = seq(from=0,by=0.00001,to=0.0045);
# r_axis_fine2 = seq(from=0,by=0.000005,to=0.001);

i0_axis = seq(from=i0_min,by=0.0005,to=i0_max)
f_axis = seq(from=f_min,by=0.001,to=f_max)


#row for i0 value
#column for r value

##need to figure out the transposition....??
# par(mfrow=c(2,1))
# image(i0_axis_fine,r_axis,as.matrix(ShiverTo30_dist))
# image(i0_axis_fine,r_axis,as.matrix(Diverge_dist))

# 
par(mfrow=c(2,1))
#If the z component is a matrix then the user should be aware that this function locates the matrix element z[i,j] 
#at the grid locations (x[i], y[j]) this is very different than simply listing out the matrix in the usual 
#row column tabular form. 

# image.plot(i0_axis,r_axis,as.matrix(ShiverTo30_dist),main='difference from temp able to warm to 30C')
# image.plot(i0_axis,r_axis,as.matrix(Diverge_dist),main='difference from temp abdomen cooling starts')
# image.plot(i0_axis,r_axis,as.matrix(ShiverTo30),main='temp able to warm to 30C')
# image.plot(i0_axis,r_axis,as.matrix(Diverge),main='temp abdomen cooling starts')
# 
# 
# image.plot(i0_axis_fine,r_axis_fine,as.matrix(ShiverTo30_dist_fine),main='difference from temp able to warm to 30C')
# image.plot(i0_axis_fine,r_axis_fine,as.matrix(Diverge_dist_fine),main='difference from temp abdomen cooling starts')
# image.plot(i0_axis_fine,r_axis_fine,as.matrix(ShiverTo30_fine),main='temp able to warm to 30C')
# image.plot(i0_axis_fine,r_axis_fine,as.matrix(Diverge_fine),main='temp abdomen cooling starts')
# 
# 
# image.plot(i0_axis_fine2,r_axis_fine2,as.matrix(ShiverTo30_dist_fine2),main='difference from temp able to warm to 30C')
# image.plot(i0_axis_fine2,r_axis_fine2,as.matrix(Diverge_dist_fine2),main='difference from temp abdomen cooling starts')
# image.plot(i0_axis_fine2,r_axis_fine2,as.matrix(ShiverTo30_fine2),main='temp able to warm to 30C')
# image.plot(i0_axis_fine2,r_axis_fine2,as.matrix(Diverge_fine2),main='temp abdomen cooling starts')


# image.plot(i0_axis,f_axis,as.matrix(ShiverTo30_dist_i0_f_1),main='difference from temp able to warm to 30C')
# image.plot(i0_axis,f_axis,as.matrix(Diverge_dist_i0_f_1),main='difference from temp abdomen cooling starts')
image.plot(i0_axis,f_axis,as.matrix(ShiverTo30_i0_f_1),main='air temp able to warm to 30C')
image.plot(i0_axis,f_axis,as.matrix(Diverge_i0_f_1),main='air temp abdomen cooling starts')


image.plot(i0_axis,f_axis,as.matrix(ShiverTo30_i0_f_1),xlim=c(i0_min,0.04), ylim=c(0.8,1),
           main='air temp able to warm to 30C')
points(0.02,0.95,col='red',pch=20)
image.plot(i0_axis,f_axis,as.matrix(Diverge_i0_f_1),xlim=c(i0_min,0.04), ylim=c(0.8,1),
           main='air temp abdomen cooling starts')
points(0.02,0.95,col='red',pch=20)

ShiverQ_min = which(ShiverTo30_i0_f_1>6,arr.ind=TRUE)
ShiverW_min = which(ShiverTo30_i0_f_1>10,arr.ind=TRUE)
ShiverQ_max = which(ShiverTo30_i0_f_1<10,arr.ind=TRUE)
ShiverW_max = which(ShiverTo30_i0_f_1<15,arr.ind=TRUE)

DivergeQ_min = which(Diverge_i0_f_1>20,arr.ind=TRUE)
DivergeQ_max = which(Diverge_i0_f_1<25,arr.ind=TRUE)
DivergeW_min = which(Diverge_i0_f_1>20,arr.ind=TRUE)
DivergeW_max = which(Diverge_i0_f_1<25,arr.ind=TRUE)

DivergeW = which(Diverge_i0_f_1>20 & Diverge_i0_f_1<25,arr.ind=TRUE)
ShiverW = which(ShiverTo30_i0_f_1<10 & ShiverTo30_i0_f_1<15,arr.ind=TRUE)


##########################################################################################
############ Calculate the euclidean distance from the target values

# target for warming is 10\dC
shiver_target = 10
shiver_difference = ShiverTo30_i0_f_1-shiver_target

# target for divergence is 20-25\dC, so use midpoint
diverge_target = 22.5
diverge_difference = Diverge_i0_f_1-diverge_target

Euclid_dist = sqrt(shiver_difference^2+diverge_difference^2)
par(mfrow=c(1,1))
image.plot(i0_axis,f_axis,as.matrix(Euclid_dist),
           main='Euclidean distance from target  (10C, 22.5C)')
image.plot(i0_axis,f_axis,as.matrix(Euclid_dist),xlim=c(i0_min,0.04), ylim=c(0.8,1),
           main='Euclidean distance from target  (10C, 22.5C)')

candidates0 = which(Euclid_dist==0,arr.ind=TRUE)
#rows for i0 and colums for f
i0_values0 = i0_axis[candidates0[,1]]
f_values0 = f_axis[candidates0[,2]]
plot(i0_values0,f_values0, main = 'candidates that meet the criteria exactly')

candidates1 = which(Euclid_dist<=1,arr.ind=TRUE)
#rows for i0 and colums for f
i0_values1 = i0_axis[candidates1[,1]]
f_values1 = f_axis[candidates1[,2]]
plot(i0_values1,f_values1, main = 'candidates within a distance of 1 from the criteria')

candidates5 = which(Euclid_dist<=5,arr.ind=TRUE)
#rows for i0 and colums for f
i0_values5 = i0_axis[candidates5[,1]]
f_values5 = f_axis[candidates5[,2]]
plot(i0_values5,f_values5, main = 'candidates within a distance of 5 from the criteria')
