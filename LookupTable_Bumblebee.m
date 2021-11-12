Temps = zeros(66,3); %rows for temperature, colums for restin/shivering/flying

for i = 1:66  %go through temps -5 to 60 
    for j = 1:3    %go through the three metabolic states
% j = 3;     %temporarily do just flying       
%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%
T_aC = i-6;    %varying air temp
T_aK = T_aC+273.15;     %air temp in K
P = 332.3878; %mean solar irradiance from all of Arrian's data
T_gC = 11;                  %ground surface temp in C
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
kappa = 0.024;   %ish for 10-15C     %thermal conductivity of air (fill in an equation )
nu = 2.791*10^(-7)*T_aK^(0.7355)/1.225;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%
indicator = j;   %is the bee resting/shivering/flying
A_th = 9.218*10^(-5)  ; %thorax surface area in m^2, from Church1960
A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
M_b = 0.149;   %mass of the bee in g, Joos1991
M_th = 0.057;  %mass of thorax in g, Joos1991
M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 4.184*0.8;  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
v_options = [0 0 3.1];   %make the bee be out of wind when resting/shivering
v = v_options(indicator);
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
C_l = 2.43*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.98;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
E = 0.63;    %Brown2004 activation energy
I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
%I_flying = 0.2097035;     %Heinrich1975, converted to W
I_flying = 0.06229515;     %Kammer1974, converted to W
%I_flying = 0.0025;     %temporarily playing around with value of I_flying to see effects
%masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
%masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
%RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
delta_T_h = 3;

r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
%r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour

%%%%%% Solar Radiation %%%%%%%
% Does not depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation thorax
S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

%%%%%% Thermal Radiation %%%%%%%
% Term 2 does depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee thorax (from sky + earth)
R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee thorax
R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head


%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;
C1 = (h*A_th*0.9)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(indicator);
T_ref = RefTemps(indicator);
norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3

% i0_resting = exp(log(I_resting) - (3/4)*log(M_ref) + E/(k*(T_ref)));  %fit i_0 to the data (I_resting has units W)
% i0_flying = exp(log(I_flying) - (3/4)*log(M_ref) + E/(k*(T_ref))); %fit i_0 to the data (I_flying has units W)
% i_0 = i0_flying;
%norm_constants = [i0_resting, i0_flying, i0_flying];   %resting/thermoregulating/flying = 1,2,3

i_0 = norm_constants(indicator);

%indicator = randsample(3,length(StartOrdered),true); %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)

%I = (i_0*(M_b^(3/4))*exp((-E/(k*T_aK))))*(1/(M_th*c)); %uses i0 fit, depends on T_air, ref temp/mass in i0 fit
%I = ( i_0 *  ((M_b/M_ref)^(3/4)) * exp(E/(k*T_ref)) * exp(-E/(k*T_aK)) ) * (1/(M_th*c)); %uses i0=I, depends on T_air
%I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass
%I = (i_0*(M_b)^(3/4))*(1/(M_th*c)); %uses i0 fit, depends on T_th, ref temp/mass in i0 fit
I = i_0;

%%%%%%%%%% Transfer to rest of body %%%%%%%%%%%%%%%
Ab = r*(1/(M_th*c));   %for heatsink version and both ways version
%Ab = r*(1/(M_th*c))*9;      %heatsink and Tairs version %suppose for a minute, T_th-T_ab=9 always
%Ab = r*(1-(T_aK/323.15))*(1/(M_th*c)); %for Tairs version 


%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed
maxy=70+273.15;
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode45 stop solving at T_th=maxy

%the times at which the data for the time dependent funcions was measured

pcolors=['b' 'y' 'r'];
%tiledlayout(2,2);  %to plot multiple panes, run this once before plotting
%any
%nexttile;  %to plot multiple panes, this has to go before each plot
tspan = [0,1000]; 
y0 = 39+273.15;   %initial temperature of the bee's head in K (observed by me! ish)
[t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0,Opt); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_alwaysactive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,Ab,T_aK),tspan,y0,Opt); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, Opt); %for use with constant environmental conditions
plot(t,y-273.15,'color',pcolors(indicator));   %plot in celsius 'k:' creates a black dotted line
hold on;
title('Bumblebee');
subtitle(['P=',num2str(P),'W/m^2'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(50);   %lethal thorax temp Heinrich1976?ish
yline(43,'k:');
legend('resting','shivering','flying','Location','southeast','color',['b' 'y' 'r'])

n_iterations = length(y);
last_its = n_iterations-20;
%Temps(i,j) = mean(y(last_its:n_iterations))-273.15;
Temps(i,j) = y(length(y))-273.15;
[i j]
    end
end
%writematrix(Temps,'LookupTable_Bumblebee_Active.csv');

figure %create a new plot window
plot(-5:60,Temps)
title('Thorax temp')
xlabel('Air Temperature (C)')
ylabel('Equilibrium Thorax Temperature (C)')
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(50);   %lethal thorax temp Heinrich1976?ish
yline(43,'k:');
hold on;
air=-5:60;
thorax=air;
plot(air,thorax,'k:')
hold off;
legend('resting','shivering','flying','T_{th} = T_{air}','Location','southeast')


%% Run this for histogram if overplotting onto variability plot
% ax2 = axes(tplot);  %Create an axes object ax2 by calling the axes function and specifying t as the parent object.
% h = histogram(equilibria_flying_BB-273.15,'orientation','horizontal');  %establish the histogram plot
% set(get(h,'Parent'),'xdir','r')   %put it on the other axis
% ax2.XAxisLocation = 'top'; %Move the x-axis to the top,  
% ax2.YAxisLocation = 'right';   %and move the y-axis to the right.
% set(ax2(1),'XTickLabel','');   %but don't actually show the axis labels
% set(ax2(1),'YTickLabel','');
% xlim([0,800])  %reset the x-limits
% ax2.Color = 'none';  %Set the color of the axes object to 'none' so that the underlying plot is visible.
% ax1.Box = 'off';    %Turn off the plot boxes to prevent the box edges from obscuring the x- and y-axes. 
% ax2.Box = 'off';

% figure
% DayTemps=readmatrix('SampleWeatherDay.csv');
% EquilibTemps = readmatrix('LookupTable_Bumblebee_Active.csv');
% BeeTemps=EquilibTemps(DayTemps(:,2),3);
% plot(BeeTemps)
% title('28 August 2021 (Dublin)')
% subtitle('Flying Bumblebee')
% xlabel('Time')
% ylabel('Equilibrium Thorax Temperature')
