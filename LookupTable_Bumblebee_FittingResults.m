FittedValues = readtable('Fitted_i0_s_5.csv','ReadVariableNames',true);

Temps = zeros(51,6); %rows for temperature, colums for resting/shivering/flying
TimeTo45 = zeros(51,7);
TimeTo30 = zeros(51,7);
solve_TimeTo45 = false;
solve_TimeTo30 = false;

TimeTo45(:,7)=0:50;  %stick air temp in the last column so I don't have to edit

%%% for plotting
air=0:50;
thorax=air;

%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%

P = 332.3878; %mean solar irradiance from all of Arrian's data
T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
%kappa = 0.02534;   %15C dry air at sea level    %thermal conductivity of air (fill in an equation )
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%
%E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
E_opts = [0.63 0 0]; 

A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
M_b = 0.149;   %mass of the bee in g, Joos1991, default
%M_b = 0.035;   %mass of the bee in g, Joos1991, minimum
%M_b = 0.351;   %mass of bee in g, Joos1991, maximum
M_th = 0.057;  %mass of thorax in g, Joos1991
M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
v_options = [0 0 4.1];   %make the bee be out of wind when resting/shivering, default
%v_options = [0 0 5.5];   %make the bee be out of wind when resting/shivering, maximum (Osborne2013)
%v_options = [0 0 1];   %make the bee be out of wind when resting/shivering, minimum (Osborne2013)
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
T_mK = 42+273.15;     %median temp for  abdomen cooling
%s = 0.9;  %fraction of internal temp at surface - default value

I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
%I_flying = 0.06229515;     %Kammer1974, converted to W
%I_flying = 0.2097035;     %Heinrich1975, converted to W
%masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
%masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
%RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)

%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour

%for i = 1:31

for fit_index = 1:height(FittedValues)
s = FittedValues.f_values5(fit_index);  %fraction of internal temp at surface - fitted value
I_flying = FittedValues.i0_values5(fit_index);     %fitted value for active metabolic rate

for i = 1:51  %go through temps 0 to 50 
%    for j = 1:6
    for j = [6]
%    for j = 3:6
%    for j = [1 3 4 6]    %go through the three metabolic states (only passive j=1:3; only active j=4:6)
% j = 3;     %temporarily do just flying       

metabolic_state=[1,2,3,1,2,3]; %go through each metabolic state (1=resting, 2=shivering, 3=flying)
model_state=[1,1,1,2,2,2]; %1 for physiological, 2 for behavioural

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%
T_aC = i-1;    %varying air temp
T_aK = T_aC+273.15;     %air temp in K
%nu = 2.791*10^(-7)*T_aK^(0.7355)/1.2256;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
rho = 1.2256;  %in dry air air at sea level at 15C
nu = mu/rho; %kinematic viscosity %The Shock Absorber Handbook
kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));

%%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%
metabolic_indicator = metabolic_state(j);   %is the bee resting/shivering/flying
model_indicator=model_state(j);       %are we using the active or the passive model

v = v_options(metabolic_indicator);


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
C1 = (h*A_th*s)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(metabolic_indicator);
T_ref = RefTemps(metabolic_indicator);
norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3
i_0 = norm_constants(metabolic_indicator);
E = E_opts(metabolic_indicator);

% i0_resting = exp(log(I_resting) - (3/4)*log(M_ref) + E/(k*(T_ref)));  %fit i_0 to the data (I_resting has units W)
% i0_flying = exp(log(I_flying) - (3/4)*log(M_ref) + E/(k*(T_ref))); %fit i_0 to the data (I_flying has units W)
% i_0 = i0_flying;
%norm_constants = [i0_resting, i0_flying, i0_flying];   %resting/thermoregulating/flying = 1,2,3


%metabolic_indicator = randsample(3,length(StartOrdered),true); %metabolic_indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)

%I = (i_0*(M_b^(3/4))*exp((-E/(k*T_aK))))*(1/(M_th*c)); %uses i0 fit, depends on T_air, ref temp/mass in i0 fit
%I = ( i_0 *  ((M_b/M_ref)^(3/4)) * exp(E/(k*T_ref)) * exp(-E/(k*T_aK)) ) * (1/(M_th*c)); %uses i0=I, depends on T_air
I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass
%I = (i_0*(M_b)^(3/4))*(1/(M_th*c)); %uses i0 fit, depends on T_th, ref temp/mass in i0 fit
%I = i_0;

%%%%%%%%%% Transfer to rest of body %%%%%%%%%%%%%%%
Ab = r*(1/(M_th*c));   %for heatsink version and both ways version
%Ab = r*(1/(M_th*c))*9;      %heatsink and Tairs version %suppose for a minute, T_th-T_ab=9 always
%Ab = r*(1-(T_aK/323.15))*(1/(M_th*c)); %for Tairs version 


%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed
maxy=50+273.15;    %don't solve above 50C because it's not biologically relevant
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode45 stop solving at T_th=maxy


pcolors=['b' 'y' 'r'];
%tiledlayout(2,2);  %to plot multiple panes, run this once before plotting
%any
%nexttile;  %to plot multiple panes, this has to go before each plot
tspan = 0:2000; 
%tspan = 0:1000;  %the time over which I have data
%y0 = 39+273.15;   %initial temperature of the bee's thorax in K (observed by me! ish)
%y0 = T_aK;
y0 = 30+273.15;   %initial temperature of the bee's thorax in K (observed by me! ish)

%solve and plot solution curves separately for passive or active model 

if model_indicator == 1 %passive model
%    figure(1)
    %[t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0,Opt); %for use with constant environmental conditions
    %[t,y] = ode45(@(t,y) heatfluxhead_Tth_alwaysactive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0); %for use with constant environmental conditions
    %[t,y] = ode45(@(t,y) heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,Ab,T_aK),tspan,y0,Opt); %for use with constant environmental conditions
    [t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, Opt); %for use with constant environmental conditions
%    plot(t,y-273.15,'color',pcolors(metabolic_indicator));   %plot in celsius 'k:' creates a black dotted line
    if solve_TimeTo45 == true
        Tmax45=45+273.15;
        OptT45=odeset('Events',@(t,y)myEvent(t,y,Tmax45));
        [t45,y45] = ode45(@(t45,y45) heatfluxhead_Tth_passive(t45,y45,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, OptT45); %for use with constant environmental conditions
        TimeTo45(i,j) = t45(end);
    end
    if solve_TimeTo30 == true
        Tmax30=30+273.15;
        OptT30=odeset('Events',@(t,y)myEvent(t,y,Tmax30));
        [t30,y30] = ode45(@(t30,y30) heatfluxhead_Tth_passive(t30,y30,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, OptT30); %for use with constant environmental conditions
        TimeTo30(i,j) = t30(end);
    end
else  %active model
%    figure(2)
    [t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK), tspan, y0,Opt); %for use with constant environmental conditions
%    plot(t,y-273.15,'color',pcolors(metabolic_indicator));   %plot in celsius 'k:' creates a black dotted line
%    plot(t,y-273.15,'color',0.7*[1,1,1]);   %plot in celsius 'k:' creates a black dotted line
    if solve_TimeTo45 == true
        Tmax45=50+273.15;
        OptT45=odeset('Events',@(t,y)myEvent(t,y,Tmax45));
        [t45,y45] = ode45(@(t45,y45) heatfluxhead_Tth_active(t45,y45,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK), tspan, y0,OptT45); %for use with constant environmental conditions
        TimeTo45(i,j) = t45(end);
    end
    if solve_TimeTo30 == true
        Tmax30=30+273.15;
        OptT30=odeset('Events',@(t,y)myEvent(t,y,Tmax30));
        [t30,y30] = ode45(@(t30,y30) heatfluxhead_Tth_active(t30,y30,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK), tspan, y0,OptT30); %for use with constant environmental conditions
        TimeTo30(i,j) = t30(end);
    end
end

hold on; 

if length(y)<length(tspan)
    y((length(y)+1):length(tspan))=50+273.15;   %if it got cut off with T_th>50, fill in the rest of y
end

% n_iterations = length(y);
% last_its = n_iterations-20;
if max(y) == 50+273.15   %if temp reached 70C, use that
    Temps_median(i,j) = 50;
    Temps_mean(i,j) = 50;
else  %otherwise, take the mean/median of the end
Temps_median(i,j) = median(y(1000:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
Temps_mean(i,j) = mean(y(1000:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
%Temps(i,j) = mean(y(last_its:n_iterations))-273.15;
%Temps(i,j) = y(length(y))-273.15;
end

[i j fit_index]

    end
end
% writematrix(Temps_median,'LookupTable_Bumblebee_default.csv');
% writematrix(TimeTo45,'TimeTo50_Bumblebee_Opt1.csv');

TimeTo45_min = TimeTo45/60;

% %add the labels to the solution curve plots
% figure(1)
% title('Physiological Model');
% subtitle(['P=',num2str(P),'W/m^2'])
% xlabel('Time (s)') ;
% ylabel('Thorax Temperature (C)') ;
% yline(30);  %30 is the min thorax temp for flight, Heinrich1983
% yline(45);   %lethal thorax temp Heinrich1976?ish
% ylim([15 70])
% %yline(43,'k:');
% %legend('resting','','flying','Location','northeast')
% legend('resting','flying','Location','northeast')
% 
% figure(2)
% title('Behavioural Model');
% subtitle(['P=',num2str(P),'W/m^2'])
% xlabel('Time (s)') ;
% ylabel('Thorax Temperature (C)') ;
% yline(30);  %30 is the min thorax temp for flight, Heinrich1983
% yline(45);   %lethal thorax temp Heinrich1976?ish
% ylim([15 70])
% %yline(43,'k:');
% %legend('resting','','flying','Location','northeast')
% legend('resting','flying','Location','northeast')


%plot thorax temp as a function of air temp
figure(1) %create a new plot window for the lookup table plot
hold on
plot(air,Temps_median(:,6),'r-')  %active, flying
alpha(0.1)

end

%plot thorax temp as a function of air temp
figure(1) %create a new plot window for the lookup table plot
ylim([20,45])
%title('Thorax temp - median')
xlabel('Air Temperature (C)')
ylabel('Equilibrium Thorax Temperature (C)')
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(45);   %lethal thorax temp Heinrich1976?ish
II=area([0,50],[45 45],42,'EdgeColor', 'none', 'FaceColor', 'r');
alpha(0.1)
fill(42,45, [0.9 0.9 0.9])
%hold off;
legend('flying','Location','southeast');
