###Process the Data Set###
#Data = read.csv("~/SUSPOLL/Code/SUSPOLLdata.csv")
Data = read.csv("~/SUSPOLL/Code/HiveActivityWeatherdataset.csv")
attach(Data)

# StartH = as.numeric(substr(obsstart,start=0,stop=2))    #extract the starting hour
# StartM = as.numeric(substr(obsstart,start=5,stop=6))    #extract the starting minutes

StartH = as.numeric(substr(timestart,start=0,stop=2))    #extract the starting hour
StartM = as.numeric(substr(timestart,start=5,stop=6))    #extract the starting minutes

Start = 60*StartH+StartM         #convert hours and minutes into minutes


# ### Establish plotting colors ###
# rainbow_ids = rainbow(10)
# colors = numeric(length(date))
# for(i in 1:10){
# colors[which(date==levels(as.factor(date))[i])]=rainbow_ids[i]
# }


### Set the values and functions ###
C_th = function(v){
	(1.93*sqrt(v)+1.4)*10^(-3)
}
C_h = function(v){
	(0.81*sqrt(v)+0.34)*10^(-3)
}
C_ab = function(v){
	C_h(0)*C_th(v)/(6.3*10^(-4))
}

delta = 5.31*10^(-13)
sigma = 5.67*10^(-8)
k = 8.617333262145*10^-5    #Boltzmann's constant
i_0 = exp(19.75)      #normalization constant from Brown2004
E = 0.63     #activation energy from Brown2004

T_aC = meantempstation    #air temp in C
T_sC = 9                  #ground surface temp in C, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
T_aK = T_aC+273.15    #air temp in K
T_sK = T_sC+273.15                  #ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report

T_th = 39      ##lowest value for departing bees in Cooper1985
T_ab = 29.5      #lowest value for departing bees in Cooper1985
T_h = 34      #lowest value for departing bees in Cooper1985
T_bC = mean(c(T_th,T_ab,T_h))      
T_bK = T_bC+273.15    #bee temp in K
T_thK = 39+273.15

l_b = 0.0010 # body length of bee in m

A_b = 2.46*10^(-5)+4.15*10^(-5)+9*10^(-5)     #sum of h, th, ab reported in Cooper1985
A_th = 4.5*10^(-5)     #in m^2, from Cooper1985
M_b = 0.100   #mass of the bee in g (the default used by Cooper1985 is 100mg)
M_th = 0.0407  #mass of thorax in g (mean reported in Cooper1985)

v = 3.1    #using M data for 100mg bee flying at 3m/s
windspeed = wind
windspeed[which(windspeed==1)]= mean(c(0.5,1.5))
windspeed[which(windspeed==2)]= mean(c(1.6,3.3))
windspeed[which(windspeed==3)]= mean(c(3.4,5.5))
windspeed[which(windspeed==4)]= mean(c(5.5,7.9))
windspeed[which(windspeed==5)]= mean(c(8,10.7))

### Calculate Some things ###

#solar radiation
#S = 0.25*0.91*A_th*meansolarstation
S = (0.25*0.91*A_th*P) + (0.5*0.91*A_th*0.25*P)   #%solar radiation

#thermal radiation
R = 0.5*A_th*(delta*T_aK^6 + sigma*T_sK^4) - 0.97*A_th*sigma*T_bK^4  #temps need to be in K

#convective heat loss - theoretical 
kin.viscosity =  function(T){       #function of temp in K (nu)
	#read a table to extract the right entry?
	return(1.45*10^(-5))   #ish for 10-15C
	}
therm.conduct =  function(T){       #function of temp in K (k)
	#read a table to extract the right entry?
	return(0.024)   #ish for 10-15C
	}
C_l = 0.0749
n = 0.78	
l = 5.5*10^(-3)   #characteristic dimension (diameter) in m
l_th = 0.0053
h_flying = ((therm.conduct(T_aK)*C_l)/l_th)*((v*l_th)/(kin.viscosity(T_aK)))^n
h_resting = ((therm.conduct(T_aK)*C_l)/l_th)*((v*l_th)/(kin.viscosity(T_aK)))^n #assumes forced airflow, so probably shouldn't actually be used for "resting" bee
C_flying = h_flying*A_th*(T_thK-T_aK)
C_resting = h_resting*A_th*(T_thK-T_aK) 

# #convective heat loss - empirical, using conductance
# C_flying = C_th(v)*(T_th-T_aC) + C_ab(v)*(T_ab-T_aC)+C_h(v)*(T_h-T_aC)  #temps need to be in C
# C_resting = C_th(windspeed)*(T_th-T_aC) + C_ab(windspeed)*(T_ab-T_aC)+C_h(windspeed)*(T_h-T_aC)  #temps need to be in C

#metabolic heat production 
#M = 0.037 from Cooper1985
I_resting = 1.31*10^(-3)     #From Rothe1989, in mW/g
i0_resting = exp(log(I_resting) - (3/4)* log(M_b) + E/(k*(10+273.15)))
I_flying = 0.4     #From Nachtigall1989, in W/g
i0_flying = exp(log(I_flying) - (3/4)* log(M_b) + E/(k*(20+273.15)))

#M = i_0*M_b^(3/4)*exp(-E/(k*T_aK))    #mass in g and temp in K

#M_resting = I_resting*M_b
#M_flying = I_flying*M_b
M_resting = i0_resting*M_b^(3/4)*exp(-E/(k*T_aK))    #mass in g and temp in K
M_flying = i0_flying*M_b^(3/4)*exp(-E/(k*T_aK))    #mass in g and temp in K

Q_flying = S+R-C_flying+M_flying*M_b
Q_resting = S+R-C_resting+M_resting*M_b

min_f = min(Q_flying)-.05
max_f = max(Q_flying)+.05
min_r = min(Q_resting)-.05
max_r = max(Q_resting)+.05
min_rf = min(min_f,min_r)
max_rf = max(max_r,max_f)

par(mfrow=c(1,2))
plot(Start,Q_flying,ylim = c(min_rf,max_rf),abline(h=0,col="gray"), pch=site-3,col=colors,xlab="Time (minutes)", ylab = "Heat Flux (W)", main = "100mg bee \n flying at 3.1m/s")
plot(Start,Q_resting,ylim = c(min_rf,max_rf),abline(h=0,col="gray"), pch=site-3,col=colors,xlab="Time (minutes)", ylab = "Heat Flux (W)", main = "100mg bee at rest")

#dev.new()
par(mfrow=c(1,2))
plot(meantempstation,Q_flying,ylim = c(min_rf,max_rf),abline(h=0,col="gray"), pch=site-3,col=colors,xlab="Temperature (C)", ylab = "Heat Flux (W)", main = "100mg bee \n flying at 3.1m/s")
plot(meantempstation,Q_resting,ylim = c(min_rf,max_rf),abline(h=0,col="gray"), pch=site-3,col=colors,xlab="Temperature (C)", ylab = "Heat Flux (W)", main = "100mg bee at rest")

#dev.new()
par(mfrow=c(1,1))
plot(jitter(wind, factor = 1, amount = NULL)
,Q_resting,ylim = c(min_rf,max_rf),abline(h=0,col="gray"), pch=site-3,col=colors,xlab="Wind speed (Beaufort Scale)", ylab = "Heat Flux (W)", main = "100mg bee at rest")

#dev.new()
par(mfrow=c(1,1))
plot(0,0)
legend(0.1,1,legend=c("site 4","site 5","site 7"), pch=c(1,2,4))
legend(0.1,0.5,legend=levels(as.factor(date)), col=rainbow(10),pch=20)

##################################################################################
datefactor = as.factor(date)
nlevels = length(levels(datefactor))

### Establish plotting colors ###
rainbow_ids = rainbow(nlevels)
colors = numeric(length(date))
for(i in 1:nlevels){
  colors[which(date==levels(datefactor)[i])]=rainbow_ids[i]
}
y_lim = c(-0.13,-.02)

par(mfrow=c(5,5))
par(mar=c(2,3,3,1))
plot(Start,Q_flying,col=site, xlab="Time (minutes from midnight)", 
     ylab = "Air Temperature",pch=20, main="All",
     ylim=y_lim)
for(i in 1:nlevels){
  #  plot(Start[which(date==levels(datefactor)[i])],meantempstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=site+10, col=colors[which(date==levels(datefactor)[i])], main=levels(datefactor)[i])
  plot(Start[which(date==levels(datefactor)[i])],Q_flying[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=20, col=site[which(date==levels(datefactor)[i])],
       main=date[which(date==levels(datefactor)[i])[1]],cex=2,ylim=y_lim)
}

par(mfrow=c(5,5))
par(mar=c(2,3,3,1))
plot(Start,Q_resting,col=site, xlab="Time (minutes from midnight)", 
     ylab = "Air Temperature",pch=20, main="All",ylim=y_lim)
for(i in 1:nlevels){
  #  plot(Start[which(date==levels(datefactor)[i])],meantempstation[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=site+10, col=colors[which(date==levels(datefactor)[i])], main=levels(datefactor)[i])
  plot(Start[which(date==levels(datefactor)[i])],Q_resting[which(date==levels(datefactor)[i])],xlab="", ylab="",pch=20, col=site[which(date==levels(datefactor)[i])],
       main=date[which(date==levels(datefactor)[i])[1]],cex=2,ylim=y_lim)
}

####################################################################################




########## Solving for T_th (why did I want to do this?)
Q = 0.136-0.192
T_a = 20
T_th = (Q-S-R-C_h(3.1)*(T_h-T_a)-C_ab(3.1)*(T_ab-T_a))/C_th(3.1)+T_a