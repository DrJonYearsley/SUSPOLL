function dydt = heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,Ab,T_aK)
%Q = S + R - E - C + I - T; 
%where y is T_th 
%Has to have all the inputs listed above in the ( )
%for use with constant environmental conditions
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(1./(1+exp(-3.*(y - 316.15)))); %heatsink version or Tairs version
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab./(1+exp(-5.*(y - 315.15)));
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*((1./(1+exp(3.*(y - 316.15))))-(1./(1+exp(3.*(y - 322.15)))) );  %bothways version 
dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(y-T_aK)*(1./(1+exp(-2.*(y - 316.15)))); %heatsink version or Tairs version with T_th-T_air dependency

%logistic function for switching: 
%f(T_th) = L/(1+exp(-k*(T_th - T_th_0)))
%where L is the max value; T_th_0 is the temperature of the curve's
%midpoint; k is the steepness

%to plot the switching function
% x = 0:0.1:50;
% y = 1./(1+exp(-1.*(x - 43)));
% plot(x,y)
% ylabel('fraction of heat transfer to thorax')
% xlabel('T_{th} (celsius)')
% xlim([20,60])

%Using an indicator function for turning on abdomen heat loss
% function dydt = heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h)
% %Q = S + R - E - C + I - T; 
% %where y is T_th 
% %Has to have all the inputs listed above in the ( )
% %for use with constant environmental conditions
% dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - 1*(y>=315.15);

%315.15K = 42C
%321.15K = 48C