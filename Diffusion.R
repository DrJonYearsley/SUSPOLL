D_A = 2.06*10^-5  #diffusion coefficient of air into itself in m^2/s
a = 4
M = 1
x = 1
t = 1/a

curve( ( M/(2*sqrt(pi*D_A*t)) )*exp( -(x^2)/(4*D_A*t) ),from=-1,to=1)


P = ( M/(2*sqrt(pi*D_A*t)) )*exp( -(x^2)/(4*D_A*t) )

#standard normal distribution
#curve( (1/sqrt(2*pi))*exp(-x^2/2), from = -10, to = 10 )
