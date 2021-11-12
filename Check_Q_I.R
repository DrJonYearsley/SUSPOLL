#i_0 = 5.1793*10^9
I_flying = 0.05918039;     #Kammer1974, converted to W
m_b = 0.149
m_0 = 0.177  #for Heinrich flying data
#m_b = seq(0.1,0.5,0.01)
E = 0.63
k = 8.617333262145*10^(-5)
T_th = seq(0,40)+273.15
T_0=25+273.15   #Kammer
#T_th = 25+273.15

Q_I = i_0*m_b^(3/4)*exp(-E/(k*T_th))
plot(T_th-273.15,Q_I)
Q_I_refs = I_flying*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/T_th)-(1/T_0)))

Q_I_Kammer = 0.06229515*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/T_th)-(1/(25+273.15))))
Q_I_Heinrich = 0.08358812*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/T_th)-(1/(39.5+273.15))))

curve(x*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/T_th)-(1/T_0))),from=0,to=1,
      xlab='I_flying',ylab='Q_I')

curve(0.08358812*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/T_th)-(1/x))),from=25+273.15,to=40+273.15,
      xlab='T_0',ylab='Q_I')

curve(0.06229515*(m_b/m_0)^(3/4)*exp(-(E/k)*((1/x)-(1/298.15))),from=25+273.15,to=70+273.15,
      xlab='T_th',ylab='Q_I',col=2)
curve((0.02)*(x-(14.5+273.15))*(1/(1+exp(-2*(x - (43+273.15)))) ), from=25+273.15,to=70+273.15,
      add=T,col=4)

#f(T_th) = L/(1+exp(-k*(T_th - T_th_0)))

x=seq(0,60,.1)
curve(r*( (1-(T_air/50))/(1+exp(-k0*(x - T_th_0))) ),from=30,to=60)
curve(r*(-L/(1+exp(-k1*(x - T_th_1)))+L),from=0,to=60)  #goes opposite direction
curve(r*(L/(1+exp(-k0*(x - T_th_0)))-L/(1+exp(-k1*(x - T_th_1)))),from=30,to=60)
curve(r*(L/(1+exp(-k0*(x - T_th_0)))-L/(1+exp(-k1*(x - T_th_1)))),from=35,to=55)
plot(0,0,ylim=c(0,1.1),xlim=c(30,60),ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")
for(T_air in seq(0,50,5)){
  curve(r*( (1-(T_air/50))/(1+exp(-k0*(x - T_th_0))) ),from=30,to=60,add=TRUE)
  }



curve((1-(x/50))/(1+exp(-k0*(x - T_th_0))),from=0,to=60)
curve(-(1-(x/50))/(1+exp(-k1*(x - T_th_1)))+(1-(x/50)),from=0,to=60)  #goes opposite direction
curve((1-(x/50))/(1+exp(-k0*(x - T_th_0)))-(1-(x/50))/(1+exp(-k1*(x - T_th_1))),from=0,to=60)


#Heatsink
curve((r/(1+exp(-3*(x - 43))) ),from=30,to=60,ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")

#Heatsink with T_th-T_air (proxy for ab)
L=1
k0=5
k1=1
T_th_0=42
T_th_1=45
r = 0.0367
r = 0.0367/20
#curve((x-T_air)*(r/(1+exp(-3*(x - 43))) ),from=30,to=60,ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")
plot(0,0,ylim=c(-0.02,0.107),xlim=c(30,60),ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")
abline(h=0)
j=1
for(T_air in seq(0,50,5)){
  curve((x-T_air)*(r/(1+exp(-2*(x - 43))) ),from=30,to=60,add=TRUE,col=j)
  j=j+1
}

#Bothways
curve(r*(1/(1+exp(-3*(x - 43)))-1/(1+exp(-3*(x - 48)))),from=30,to=60,ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")

#Tairs
plot(0,0,ylim=c(0,0.04),xlim=c(30,60),ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")
for(T_air in seq(0,50,5)){
  curve(r*( (1-(T_air/50))/(1+exp(-3*(x - 43))) ),from=30,to=60,add=TRUE)
}

#Tairs with T_th-T_air (proxy for ab)
plot(0,0,ylim=c(0,2.5),xlim=c(35,60),ylab="rate of heat flow (Q_ab)",xlab="thorax temperature")
j=1
for(T_air in seq(0,50,5)){
  curve((x-T_air)*r*( (1-(T_air/50))/(1+exp(-3*(x - 43))) ),from=35,to=60,add=TRUE,col=j)
  j=j+1
}
