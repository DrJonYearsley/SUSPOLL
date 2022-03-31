function dydt = heatfluxhead_BB_CoolingOn(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK)
%This function is for bumblebees
%with abdomen cooling always on 
dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(y-T_aK); %heatsink version or Tairs version

