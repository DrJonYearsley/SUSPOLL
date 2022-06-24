Data = read.csv("~/SUSPOLL/Code/HiveActivityWeatherdataset.csv")
attach(Data)

StartH = as.numeric(substr(timestart,start=0,stop=2))    #extract the starting hour
StartM = as.numeric(substr(timestart,start=5,stop=6))    #extract the starting minutes

Start = 60*StartH+StartM         #convert hours and minutes into minutes

windspeed = wind
windspeed[which(windspeed==1)]= mean(c(0.5,1.5))
windspeed[which(windspeed==2)]= mean(c(1.6,3.3))
windspeed[which(windspeed==3)]= mean(c(3.4,5.5))
windspeed[which(windspeed==4)]= mean(c(5.5,7.9))
windspeed[which(windspeed==5)]= mean(c(8,10.7))

datefactor = as.factor(date)
nlevels = length(levels(datefactor))

### Establish plotting colors ###
rainbow_ids = rainbow(nlevels)
colors = numeric(length(date))
for(i in 1:nlevels){
  colors[which(date==levels(datefactor)[i])]=rainbow_ids[i]
}

par(mfrow=c(5,5))
par(mar=c(2,3,3,1))
plot(Start,meantempstation,col=site, xlab="Time (minutes from midnight)", 
     ylab = "Air Temperature",pch=20, main="All",ylim=c(7,21))
for(i in 1:nlevels){
#  plot(Start[which(date==levels(datefactor)[i])],meantempstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=site+10, col=colors[which(date==levels(datefactor)[i])], main=levels(datefactor)[i])
  plot(Start[which(date==levels(datefactor)[i])],meantempstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=20, col=site[which(date==levels(datefactor)[i])],
       main=date[which(date==levels(datefactor)[i])[1]],cex=2,
       ylim=c(7,21))
}

par(mfrow=c(5,5))
par(mar=c(2,3,3,1))
plot(Start,meansolarstation,col=colors, xlab="Time (minutes from midnight)", ylab = "Air Temperature",pch=site+10, main="All")
for(i in 1:nlevels){
  plot(Start[which(date==levels(datefactor)[i])],meansolarstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=site+10, col=colors[which(date==levels(datefactor)[i])], main=levels(datefactor)[i])
}


plot(Start[which(date==levels(datefactor)[1])],meantempstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=20, col=site)
