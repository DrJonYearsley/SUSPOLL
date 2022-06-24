function dydt = heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2)
%This function is for bumblebees
%with abdomen cooling always on 

if Honeybee == true 
    X_sfc = (rh*A*exp(B/(y-delta_T_h)))/Pr;
    Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(MM_air/MM_vapor) );
else 
    Y_sfc = 0;
end

if CoolingSwitch_indicator == 1
    f = (1./(1+exp(-3.*(y - T_mK))));
else
    f = 1;
end

dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))...   %thorax baseline
    + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h)... %head baseline
    - Ab*(y-T_aK)*f...   %abdomen 
    - (Ev1+Ev2*log( 1 - Y_sfc ))*f;   %evaporative

