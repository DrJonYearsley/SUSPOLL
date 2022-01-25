function dydt = heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK)
%This function is for bumblebees
%Q = S + R - E - C + I - T; 
%where y is T_th 
%Has to have all the inputs listed above in the ( )
%for use with constant environmental conditions
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(1./(1+exp(-3.*(y - 316.15)))); %heatsink version or Tairs version
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - 0.0367.*( 1./(1+exp(-1.*(y - 315.15))) );
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab./(1+exp(-2.*(y - 315.15)) );
dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - (Ab*(y-T_aK)*(1./(1+exp(-3.*(y - T_mK))))); %heatsink version or Tairs version
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(2*(y-T_aK)*(1./(1+exp(-2.*(y - 315.15))))+(y-T_aK)); %heatsink version or Tairs version
%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I)  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(y-T_aK)*(1./(1+exp(-2.*(y - 316.15)))); %just to try a constant I...

%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - Ab*(
%(1./(1+exp(1.*(y - 315.15))))-(1./(1+exp(1.*(y - 322.15)))) );  %bothways
%version 

%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - ((1-y./50)*(1./(1+exp(-3.*(y - 315.15)))-(1-y./50)./(1+exp(-1.*(y - 321.15)))));

%dydt = (S + R1 - R2.*y.^4 - C1.*y - C2 + I*exp(-E/(k.*y)))  + (S_h + R1_h - R2_h.*(y-delta_T_h).^4 - C1_h.*(y-delta_T_h) - C2_h) - (1*(1./(1+exp(-3.*(y - 315.15)))-1./(1+exp(-1.*(y - 321.15)))));

%logistic function for switching: 
%f(T_th) = L/(1+exp(-k*(T_th - T_th_0)))
%f(T_th) = L/(1+exp(-k*(x - T_th_0)))-L/(1+exp(-k*(x - T_th_1)))
%where L is the max value; 
%T_th_0 is the temperature of the increasing curve's midpoint; 
%k is the steepness
%T_th_1 is the midpoint of the decreasing curve's midpoint
%315.15K = 42C
%318.15K = 45C
%321.15K = 48C
%323.1500 = 50C

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

