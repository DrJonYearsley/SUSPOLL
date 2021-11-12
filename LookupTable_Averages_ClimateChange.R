#library(ggplot2)

LookupTable=read.csv("C:/Users/UCD/Documents/SUSPOLL/Code/LookupTable_Bumblebee_Active.csv")
colnames(LookupTable)=c("resting","shivering","flying")


# MeraData_1 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_01.RData")
# MeraData_2 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_02.RData")
# MeraData_3 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_03.RData")
# MeraData_4 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_04.RData")
# MeraData_5 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_05.RData")
# MeraData_6 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_06.RData")
# MeraData_7 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_07.RData")
# MeraData_8 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_08.RData")
# MeraData_9 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_09.RData")
# MeraData_10 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_10.RData")
# MeraData_11 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_11.RData")
# MeraData_12 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_12.RData")

#MeraData = mera_hourly
#Day = 20160115
#Location = c(53.457,-6.140)
#mode = 3

CalculateAverages_CC = function(MeraData,LookupTable,Location,mode,increase){
  #takes in Mera climate data, lookup table of temps, a day, and a location 
  #plots the equilibrium thorax temp as a function of time over the day
  #Day should be in format yyyymmdd with no spaces
  #Location should be in format c(latitude,longitude)
  #mode is 1,2,3 for "resting","shivering","flying"
  #increase is the projected mean annual temp increase
  
  
  
  
  SubsetLat = MeraData[which(round(MeraData$Latitude,3)==Location[1]),]
  SubsetLong = SubsetLat[which(round(SubsetLat$Longitude,3)==Location[2]),]
  
  
  Index0 = which(SubsetLong$dataTime.==0)
  Index300 = which(SubsetLong$dataTime.==300)
  Index600 = which(SubsetLong$dataTime.==600)
  Index900 = which(SubsetLong$dataTime.==900)
  Index1200 = which(SubsetLong$dataTime.==1200)
  Index1500 = which(SubsetLong$dataTime.==1500)
  Index1800 = which(SubsetLong$dataTime.==1800)
  Index2100 = which(SubsetLong$dataTime.==2100)
  
  
  Subset0 = SubsetLong[Index0,] 
  Subset300 = SubsetLong[Index300,] 
  Subset600 = SubsetLong[Index600,] 
  Subset900 = SubsetLong[Index900,] 
  Subset1200 = SubsetLong[Index1200,] 
  Subset1500 = SubsetLong[Index1500,] 
  Subset1800 = SubsetLong[Index1800,] 
  Subset2100 = SubsetLong[Index2100,] 
  
  Avg_SD = matrix(NA, nrow=8,ncol=2)
  Avg_SD[1,1] = mean(Subset0$Value.)+increase      #column for mean/sd
  Avg_SD[2,1] = mean(Subset300$Value.)+increase    #column for 
  Avg_SD[3,1] = mean(Subset600$Value.)+increase    #column for 
  Avg_SD[4,1] = mean(Subset900$Value.)+increase    #column for 
  Avg_SD[5,1] = mean(Subset1200$Value.)+increase    #column for 
  Avg_SD[6,1] = mean(Subset1500$Value.)+increase    #column for 
  Avg_SD[7,1] = mean(Subset1800$Value.)+increase    #column for 
  Avg_SD[8,1] = mean(Subset2100$Value.)+increase    #column for 
  Avg_SD[1,2] = sd(Subset0$Value.)      #column for 
  Avg_SD[2,2] = sd(Subset300$Value.)   #column for 
  Avg_SD[3,2] = sd(Subset600$Value.)   #column for 
  Avg_SD[4,2] = sd(Subset900$Value.)   #column for 
  Avg_SD[5,2] = sd(Subset1200$Value.)   #column for 
  Avg_SD[6,2] = sd(Subset1500$Value.)   #column for 
  Avg_SD[7,2] = sd(Subset1800$Value.)   #column for 
  Avg_SD[8,2] = sd(Subset2100$Value.)   #column for 
  
  AirTemps = round(cbind(Avg_SD[,1]-Avg_SD[,2],Avg_SD[,1],Avg_SD[,1]+Avg_SD[,2])-273.15+increase)
  
  ThoraxTemps_MinAvgMax = matrix(NA, nrow=8,ncol=3)
  for(i in 1:8){
    ThoraxTemps_MinAvgMax[i,1] = LookupTable[AirTemps[i,1]+6,mode]   #column for mode
    ThoraxTemps_MinAvgMax[i,2] = LookupTable[AirTemps[i,2]+6,mode]   #column for mode
    ThoraxTemps_MinAvgMax[i,3] = LookupTable[AirTemps[i,3]+6,mode]   #column for mode
    #In lookup table, row n corresponds to air temp n-6
    #so if AirTemps[i]=1, need to look in row 7 of lookup table, i.e. AirTemps[i]+6
  }
  return(cbind(ThoraxTemps_MinAvgMax,AirTemps))
}

Find_Heatwaves_CC = function(MeraData,LookupTable,mode,increase){
  
  ULocations = unique(MeraData[,c('Latitude','Longitude')])
  ULocations$Location = paste0('L',seq(nrow(ULocations)))
  UDays = data.frame(dataDate.=unique(MeraData$dataDate.))
  UDays$Day = paste0('D',seq(nrow(UDays)))  
  MeraData=merge(MeraData,ULocations)
  MeraData=merge(MeraData,UDays)
  
  Heatwave_Time_Indices = which((MeraData$Value.-273.15)>25)
  temp = MeraData[Heatwave_Time_Indices,]
  
  HeatwaveDays = subset(MeraData, 
                        Location %in% temp$Location & 
                          Day %in% temp$Day)
  
  # Heatwave_LocationDays = MeraData$Day[Heatwave_Time_Indices,]
  
  
  
  
  Index0 = which(HeatwaveDays$dataTime.==0)
  Index300 = which(HeatwaveDays$dataTime.==300)
  Index600 = which(HeatwaveDays$dataTime.==600)
  Index900 = which(HeatwaveDays$dataTime.==900)
  Index1200 = which(HeatwaveDays$dataTime.==1200)
  Index1500 = which(HeatwaveDays$dataTime.==1500)
  Index1800 = which(HeatwaveDays$dataTime.==1800)
  Index2100 = which(HeatwaveDays$dataTime.==2100)
  
  
  Subset0 = HeatwaveDays[Index0,] 
  Subset300 = HeatwaveDays[Index300,] 
  Subset600 = HeatwaveDays[Index600,] 
  Subset900 = HeatwaveDays[Index900,] 
  Subset1200 = HeatwaveDays[Index1200,] 
  Subset1500 = HeatwaveDays[Index1500,] 
  Subset1800 = HeatwaveDays[Index1800,] 
  Subset2100 = HeatwaveDays[Index2100,] 
  
  Avg_SD = matrix(NA, nrow=8,ncol=2)
  Avg_SD[1,1] = mean(Subset0$Value.)+increase      #column for mean/sd
  Avg_SD[2,1] = mean(Subset300$Value.)+increase    #column for 
  Avg_SD[3,1] = mean(Subset600$Value.)+increase    #column for 
  Avg_SD[4,1] = mean(Subset900$Value.)+increase    #column for 
  Avg_SD[5,1] = mean(Subset1200$Value.)+increase   #column for 
  Avg_SD[6,1] = mean(Subset1500$Value.)+increase    #column for 
  Avg_SD[7,1] = mean(Subset1800$Value.)+increase    #column for 
  Avg_SD[8,1] = mean(Subset2100$Value.)+increase    #column for 
  Avg_SD[1,2] = sd(Subset0$Value.)      #column for 
  Avg_SD[2,2] = sd(Subset300$Value.)   #column for 
  Avg_SD[3,2] = sd(Subset600$Value.)   #column for 
  Avg_SD[4,2] = sd(Subset900$Value.)   #column for 
  Avg_SD[5,2] = sd(Subset1200$Value.)   #column for 
  Avg_SD[6,2] = sd(Subset1500$Value.)   #column for 
  Avg_SD[7,2] = sd(Subset1800$Value.)   #column for 
  Avg_SD[8,2] = sd(Subset2100$Value.)   #column for 
  
  AirTemps = round(cbind(Avg_SD[,1]-Avg_SD[,2],Avg_SD[,1],Avg_SD[,1]+Avg_SD[,2])-273.15)
  
  ThoraxTemps_MinAvgMax = matrix(NA, nrow=8,ncol=3)
  for(i in 1:8){
    ThoraxTemps_MinAvgMax[i,1] = LookupTable[AirTemps[i,1]+6,mode]   #column for mode
    ThoraxTemps_MinAvgMax[i,2] = LookupTable[AirTemps[i,2]+6,mode]   #column for mode
    ThoraxTemps_MinAvgMax[i,3] = LookupTable[AirTemps[i,3]+6,mode]   #column for mode
    #In lookup table, row n corresponds to air temp n-6
    #so if AirTemps[i]=1, need to look in row 7 of lookup table, i.e. AirTemps[i]+6
  }
  return(cbind(ThoraxTemps_MinAvgMax,AirTemps))
}

RCPno = 1.6 #degrees warming predicted
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_01.RData")
MeraData1 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_02.RData")
MeraData2 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_03.RData")
MeraData3 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_04.RData")
MeraData4 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_05.RData")
MeraData5 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_06.RData")
MeraData6 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_07.RData")
MeraData7 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_08.RData")
MeraData8 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_09.RData")
MeraData9 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_10.RData")
MeraData10 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_11.RData")
MeraData11 = mera_hourly
load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_12.RData")
MeraData12 = mera_hourly

# ThoraxTemps_MinAvgMax_1 = CalculateAverages_CC(MeraData1,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_2 = CalculateAverages_CC(MeraData2,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_3 = CalculateAverages_CC(MeraData3,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_4 = CalculateAverages_CC(MeraData4,LookupTable,c(53.457,-6.140),3,RCPno)
ThoraxTemps_MinAvgMax_5_CC = CalculateAverages_CC(MeraData5,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_6 = CalculateAverages_CC(MeraData6,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_7 = CalculateAverages_CC(MeraData7,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_8 = CalculateAverages_CC(MeraData8,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_9 = CalculateAverages_CC(MeraData9,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_10 = CalculateAverages_CC(MeraData10,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_11 = CalculateAverages_CC(MeraData11,LookupTable,c(53.457,-6.140),3,RCPno)
# ThoraxTemps_MinAvgMax_12 = CalculateAverages_CC(MeraData12,LookupTable,c(53.457,-6.140),3,RCPno)

MeraDataGrand = rbind(MeraData1,MeraData2,MeraData3,MeraData4,MeraData5,MeraData6,
                      MeraData7,MeraData8,MeraData9,MeraData10,MeraData11,MeraData12)
Heatwave_Averages_CC = Find_Heatwaves_CC(MeraDataGrand,LookupTable,3,RCPno)

# par(mfrow=c(3,4))
x=c(0,300,600,900,1200,1500,1800,2100)/100
epsilon = .5
y_range = c(0,50)
hours=3:8

# y=ThoraxTemps_MinAvgMax_1[,2]
# y_plus = ThoraxTemps_MinAvgMax_1[,3]
# y_minus = ThoraxTemps_MinAvgMax_1[,1]
# #points(x, y, ylim=c(38,43),col=rainbow(12)[1])
# plot(x, y, ylim=y_range,col=rainbow(12)[1],xlab='time',ylab='thorax temp', main='Jan')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[1])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[1])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[1])
# points(x,ThoraxTemps_MinAvgMax_1[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_1[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_1[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_2[,2]
# y_plus = ThoraxTemps_MinAvgMax_2[,3]
# y_minus = ThoraxTemps_MinAvgMax_2[,1]
# #points(x, y,col=rainbow(12)[2])
# plot(x, y, ylim=y_range,col=rainbow(12)[2],xlab='time',ylab='thorax temp', main='Feb')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[2])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[2])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[2])
# points(x,ThoraxTemps_MinAvgMax_2[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_2[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_2[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_3[,2]
# y_plus = ThoraxTemps_MinAvgMax_3[,3]
# y_minus = ThoraxTemps_MinAvgMax_3[,1]
# #points(x, y,col=rainbow(12)[3])
# plot(x, y, ylim=y_range,col=rainbow(12)[3],xlab='time',ylab='thorax temp', main='Mar')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[3])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[3])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[3])
# points(x,ThoraxTemps_MinAvgMax_3[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_3[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_3[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_4[,2]
# y_plus = ThoraxTemps_MinAvgMax_4[,3]
# y_minus = ThoraxTemps_MinAvgMax_4[,1]
# #points(x, y,col=rainbow(12)[4])
# plot(x, y, ylim=y_range,col=rainbow(12)[4],xlab='time',ylab='thorax temp', main='Apr')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[4])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[4])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[4])
# points(x,ThoraxTemps_MinAvgMax_4[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_4[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_4[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)

y=ThoraxTemps_MinAvgMax_5_CC[,2]
y_plus = ThoraxTemps_MinAvgMax_5_CC[,3]
y_minus = ThoraxTemps_MinAvgMax_5_CC[,1]
#points(x, y,col=rainbow(12)[5])
plot(x[hours], y[hours], ylim=y_range,col=rainbow(12)[5],xlab='time',ylab='temperature (C)', main='2041-2060 May Average')
segments(x[hours], y_minus[hours],x[hours], y_plus[hours],col=rainbow(12)[5])
segments(x[hours]-epsilon,y_plus[hours],x[hours]+epsilon,y_plus[hours],col=rainbow(12)[5])
segments(x[hours]-epsilon,y_minus[hours],x[hours]+epsilon,y_minus[hours],col=rainbow(12)[5])
points(x[hours],ThoraxTemps_MinAvgMax_5_CC[hours,5])      #add the air temperature
y_plus_air = ThoraxTemps_MinAvgMax_5_CC[,6]
y_minus_air = ThoraxTemps_MinAvgMax_5_CC[,4]
segments(x[hours], y_minus_air[hours],x[hours], y_plus_air[hours])
segments(x[hours]-epsilon,y_plus_air[hours],x[hours]+epsilon,y_plus_air[hours])
segments(x[hours]-epsilon,y_minus_air[hours],x[hours]+epsilon,y_minus_air[hours])

# y=ThoraxTemps_MinAvgMax_6[,2]
# y_plus = ThoraxTemps_MinAvgMax_6[,3]
# y_minus = ThoraxTemps_MinAvgMax_6[,1]
# #points(x, y,col=rainbow(12)[6])
# plot(x, y, ylim=y_range,col=rainbow(12)[6],xlab='time',ylab='thorax temp', main='Jun')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[6])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[6])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[6])
# points(x,ThoraxTemps_MinAvgMax_6[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_6[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_6[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_7[,2]
# y_plus = ThoraxTemps_MinAvgMax_7[,3]
# y_minus = ThoraxTemps_MinAvgMax_7[,1]
# #points(x, y,col=rainbow(12)[7])
# plot(x, y, ylim=y_range,col=rainbow(12)[7],xlab='time',ylab='thorax temp', main='Jul')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[7])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[7])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[7])
# points(x,ThoraxTemps_MinAvgMax_7[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_7[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_7[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_8[,2]
# y_plus = ThoraxTemps_MinAvgMax_8[,3]
# y_minus = ThoraxTemps_MinAvgMax_8[,1]
# #points(x, y,col=rainbow(12)[8])
# plot(x, y, ylim=y_range,col=rainbow(12)[8],xlab='time',ylab='thorax temp', main='Aug')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[8])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[8])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[8])
# points(x,ThoraxTemps_MinAvgMax_8[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_8[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_8[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_9[,2]
# y_plus = ThoraxTemps_MinAvgMax_9[,3]
# y_minus = ThoraxTemps_MinAvgMax_9[,1]
# #points(x, y,col=rainbow(12)[9])
# plot(x, y, ylim=y_range,col=rainbow(12)[9],xlab='time',ylab='thorax temp', main='Sep')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[9])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[9])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[9])
# points(x,ThoraxTemps_MinAvgMax_9[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_9[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_9[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_10[,2]
# y_plus = ThoraxTemps_MinAvgMax_10[,3]
# y_minus = ThoraxTemps_MinAvgMax_10[,1]
# #points(x, y,col=rainbow(12)[10])
# plot(x, y, ylim=y_range,col=rainbow(12)[10],xlab='time',ylab='thorax temp', main='Oct')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[10])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[10])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[10])
# points(x,ThoraxTemps_MinAvgMax_10[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_10[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_10[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_11[,2]
# y_plus = ThoraxTemps_MinAvgMax_11[,3]
# y_minus = ThoraxTemps_MinAvgMax_11[,1]
# #points(x, y,col=rainbow(12)[11])
# plot(x, y, ylim=y_range,col=rainbow(12)[11],xlab='time',ylab='thorax temp', main='Nov')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[11])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[11])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[11])
# points(x,ThoraxTemps_MinAvgMax_11[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_11[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_11[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)
# 
# y=ThoraxTemps_MinAvgMax_12[,2]
# y_plus = ThoraxTemps_MinAvgMax_12[,3]
# y_minus = ThoraxTemps_MinAvgMax_12[,1]
# #points(x, y,col=rainbow(12)[12])
# plot(x, y, ylim=y_range,col=rainbow(12)[12],xlab='time',ylab='thorax temp', main='Dec')
# segments(x, y_minus,x, y_plus,col=rainbow(12)[12])
# segments(x-epsilon,y_plus,x+epsilon,y_plus,col=rainbow(12)[12])
# segments(x-epsilon,y_minus,x+epsilon,y_minus,col=rainbow(12)[12])
# points(x,ThoraxTemps_MinAvgMax_12[,5])      #add the air temperature
# y_plus_air = ThoraxTemps_MinAvgMax_12[,6]
# y_minus_air = ThoraxTemps_MinAvgMax_12[,4]
# segments(x, y_minus_air,x, y_plus_air)
# segments(x-epsilon,y_plus_air,x+epsilon,y_plus_air)
# segments(x-epsilon,y_minus_air,x+epsilon,y_minus_air)










x=c(0,300,600,900,1200,1500,1800,2100)/100
epsilon = .5
hours=3:8

y=Heatwave_Averages_CC[,2]
y_plus = Heatwave_Averages_CC[,3]
y_minus = Heatwave_Averages_CC[,1]
#points(x, y, ylim=c(38,43),col=rainbow(12)[1])
plot(x[hours], y[hours], ylim=y_range,col=1,xlab='time (hours)',ylab='temperature (C)', main='2041-2060 Heatwave Day Average')
segments(x[hours], y_minus[hours],x[hours], y_plus[hours],col=1)
segments(x[hours]-epsilon,y_plus[hours],x[hours]+epsilon,y_plus[hours],col=1)
segments(x[hours]-epsilon,y_minus[hours],x[hours]+epsilon,y_minus[hours],1)
points(x[hours],Heatwave_Averages_CC[hours,5],col='blue')      #add the air temperature
y_plus_air = Heatwave_Averages_CC[,6]
y_minus_air = Heatwave_Averages_CC[,4]
segments(x[hours], y_minus_air[hours],x[hours], y_plus_air[hours],col='blue')
segments(x[hours]-epsilon,y_plus_air[hours],x[hours]+epsilon,y_plus_air[hours],col='blue')
segments(x[hours]-epsilon,y_minus_air[hours],x[hours]+epsilon,y_minus_air[hours],col='blue')







