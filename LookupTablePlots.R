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

SampleDay = function(MeraData,LookupTable,Day,Location,mode){
  #takes in Mera climate data, lookup table of temps, a day, and a location 
  #plots the equilibrium thorax temp as a function of time over the day
  #Day should be in format yyyymmdd with no spaces
  #Location should be in format c(latitude,longitude)
  #mode is 1,2,3 for "resting","shivering","flying"
  
  SubsetDate = MeraData[which(MeraData$dataDate.==Day),]
  SubsetLat = SubsetDate[which(SubsetDate$Latitude==Location[1]),]
  SubsetLong = SubsetLat[which(round(SubsetLat$Longitude,3)==Location[2]),]
  AirTemps = round(SubsetLong$Value.-273.15)
  ThoraxTemps = numeric(length(AirTemps))
  for(i in 1:length(AirTemps)){
    ThoraxTemps[i] = LookupTable[AirTemps[i]+6,mode]   #column for mode
    #In lookup table, row n corresponds to air temp n-6
    #so if AirTemps[i]=1, need to look in row 7 of lookup table, i.e. AirTemps[i]+6
  }
  return(cbind(SubsetLong$dataTime.,ThoraxTemps))
}

YearData_ThoraxTemps = matrix(NA, nrow=8,ncol=12)  #matrix to store the temps in 
# mera data seems to have recorded 8 times in the day

for(j in 1:9){
  MeraName = paste('C:/Users/UCD/Downloads/Temp_subset_hourly_2016_0',j,'.RData',sep='')
  load(MeraName)
  Day = paste('20160',j,'15',sep='')
  BeeTemps=SampleDay(mera_hourly,LookupTable,Day,c(53.457,-6.140),3)
  YearData_ThoraxTemps[,j] = BeeTemps[,2]
}
for(j in 10:12){
  MeraName = paste('C:/Users/UCD/Downloads/Temp_subset_hourly_2016_',j,'.RData',sep='')
  load(MeraName)
  Day = paste('2016',j,'15',sep='')
  BeeTemps=SampleDay(mera_hourly,LookupTable,Day,c(53.457,-6.140),3)
  YearData_ThoraxTemps[,j] = BeeTemps[,2]
}  

par(mfrow=c(1,2))
matplot(BeeTemps[,1]/100,as.data.frame(YearData_ThoraxTemps),
        type='b',lty=1,pch=20,col=rainbow(12),
        xlab = "Time (hrs)", ylab = "Equlibrium thorax temp (C)")
plot(0,0,pch='.',axes='n')
legend(-1,1,legend=format(ISOdate(2004,1:12,1),"%B"),col=rainbow(12),lty=1,pch=20)
