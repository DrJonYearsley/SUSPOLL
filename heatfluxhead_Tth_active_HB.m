function dydt = heatfluxhead_Tth_active_HB(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ev1,Ev2,rh,A,B,MM_air,MM_vapor,Pr,T_mK)
%This function is for bumblebees
%Q = S + R - E - C + I - T; 
%where y is T_th 
%Has to have all the inputs listed above in the ( )

X_sfc = (rh*A*exp(B/(y-delta_T_h)))/Pr;
Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) );

%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) + (Ev1+Ev2*log( 1 -Y_sfc ) )*(1./(1+exp(-3.*(y - T_mK))))); 
dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))...
    + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h)...
    - Ev1-Ev2*log( 1 - Y_sfc ); %if evaporative cooling is always on

%logistic function for switching: 
%f(T_th) = L/(1+exp(-k*(T_th - T_th_0)))
%f(T_th) = L/(1+exp(-k*(x - T_th_0)))-L/(1+exp(-k*(x - T_th_1)))
%where L is the max value; 
%T_th_0 is the temperature of the increasing curve's midpoint; 
%k is the steepness
%T_th_1 is the midpoint of the decreasing curve's midpoint

