function dydt = heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k)
%Q = S + R - E - C + I - T; 
%where y is T_th 
%Has to have all the inputs listed above in the ( )
%for use with constant environmental conditions
dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h);