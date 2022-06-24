LookupTable=read.csv("C:/Users/UCD/Documents/SUSPOLL/Code/LookupTable_Bumblebee_Active.csv")
colnames(LookupTable)=c("resting","shivering","flying")


# MeraData_3 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_03.RData")
# MeraData_4 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_04.RData")
# MeraData_5 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_05.RData")
# MeraData_6 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_06.RData")
# MeraData_7 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_07.RData")
# MeraData_8 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_08.RData")
# MeraData_9 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_09.RData")
# MeraData_10 = load("C:/Users/UCD/Downloads/Temp_subset_hourly_2016_10.RData")

#MeraData = mera_hourly
#Day = 20160115
#Location = c(53.457,-6.140)
#mode = 3

HighDay = function(MeraData,LookupTable,Location,mode){
  #takes in Mera climate data, lookup table of temps, a day, and a location 
  #plots the equilibrium thorax temp as a function of time over the day
  #Day should be in format yyyymmdd with no spaces
  #Location should be in format c(latitude,longitude)
  #mode is 1,2,3 for "resting","shivering","flying"
  
  
  SubsetLat = MeraData[which(round(MeraData$Latitude,3)==Location[1]),]
  SubsetLong = SubsetLat[which(round(SubsetLat$Longitude,3)==Location[2]),]
  HighIndex = which(SubsetLong$Value.==max(SubsetLong$Value.))
  SubsetDaytime = SubsetLong[which(SubsetLong$dataTime.>300),]
  LowIndex = which(SubsetDaytime$Value.==min(SubsetDaytime$Value.))
  HighDay = SubsetLong[HighIndex,4]
  LowDay = SubsetLong[LowIndex,4]
  SubsetDate_High = SubsetLong[which(SubsetLong$dataDate.==HighDay),]  
  SubsetDate_Low = SubsetLong[which(SubsetLong$dataDate.==LowDay),]  
  AirTemps_High = round(SubsetDate_High$Value.-273.15)
  AirTemps_Low = round(SubsetDate_Low$Value.-273.15)
  ThoraxTemps = matrix(NA, nrow=length(AirTemps_High),ncol=2)
  for(i in 1:length(AirTemps_High)){
    ThoraxTemps[i,1] = LookupTable[AirTemps_Low[i]+6,mode]   #column for mode
    ThoraxTemps[i,2] = LookupTable[AirTemps_High[i]+6,mode]   #column for mode
    #In lookup table, row n corresponds to air temp n-6
    #so if AirTemps[i]=1, need to look in row 7 of lookup table, i.e. AirTemps[i]+6
  }
  return(cbind(SubsetDate_Low$dataTime.,ThoraxTemps))
}

YearData_ThoraxTemps_High = matrix(NA, nrow=8,ncol=12)  #matrix to store the bee temps in 
YearData_ThoraxTemps_Low = matrix(NA, nrow=8,ncol=12)  #matrix to store the bee temps in 
# mera data seems to have recorded 8 times in the day

for(j in 1:9){
  MeraName = paste('C:/Users/UCD/Downloads/Temp_subset_hourly_2016_0',j,'.RData',sep='')
  load(MeraName)
  BeeTemps=HighDay(mera_hourly,LookupTable,c(53.457,-6.140),3)
  YearData_ThoraxTemps_Low[,j] = BeeTemps[,2]
  YearData_ThoraxTemps_High[,j] = BeeTemps[,3]
}
for(j in 10:12){
  MeraName = paste('C:/Users/UCD/Downloads/Temp_subset_hourly_2016_',j,'.RData',sep='')
  load(MeraName)
  BeeTemps=HighDay(mera_hourly,LookupTable,c(53.457,-6.140),3)
  YearData_ThoraxTemps_Low[,j] = BeeTemps[,2]
  YearData_ThoraxTemps_High[,j] = BeeTemps[,3]
}  

par(mfrow=c(1,2))
matplot(BeeTemps[,1]/100,as.data.frame(YearData_ThoraxTemps_Low),
        type='b',lty=1,pch=20,col=rainbow(12), ylim = c(38,43),
        xlab = "Time (hrs)", ylab = "Equlibrium thorax temp (C)",main='Monthly Low')
matplot(BeeTemps[,1]/100,as.data.frame(YearData_ThoraxTemps_High),
        type='b',lty=1,pch=20,col=rainbow(12), ylim = c(38,43),
        xlab = "Time (hrs)", ylab = "Equlibrium thorax temp (C)",main='Monthly High')
#plot(0,0,pch='.',axes='n')
#legend(-1,1,legend=format(ISOdate(2004,1:12,1),"%B"),col=rainbow(12),lty=1,pch=20)
 