%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ArrianDataSet = readtable('HiveActivityWeatherdataset.csv','ReadVariableNames',true);
WeatherDateIndices = find(ismember(ArrianDataSet.species, 'H') & ArrianDataSet.site==2 & ArrianDataSet.date=='28/05/2019'); 
%use weather data for site 1 on April 26th 2019 (or whenever)
ShortDataSet = ArrianDataSet(WeatherDateIndices,:);  %: indicates all colums

%time of each observation in minutes from midnight
StartH_string = cellfun(@(x) x(1:2),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting hour
StartM_string = cellfun(@(x) x(5:6),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting minutes
StartH = cellfun(@str2double,StartH_string);
StartM = cellfun(@str2double,StartM_string);
Start = 60*StartH+StartM;         %convert hours and minutes into minutes

%Keep just the entries with unique times
[StartOrdered,StartIndex] = unique(Start); 
ShortDataSetOrdered = ShortDataSet(StartIndex,:); 

%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%
%T_aC = ShortDataSetOrdered.meantempstation;    %air temp in C from Arrian's data
%P = ShortDataSetOrdered.meansolarstation; %solar irradiance from Arrian's data
T_aC = mean(ShortDataSetOrdered.meantempstation);    %mean air temp in C from Arrian's data
%T_aC = 5;    %to try different air temps
P = mean(ShortDataSetOrdered.meansolarstation); %mean solar irradiance from Arrian's data
%beta = 48.52; %angle of the sun on average (??? for Dublin, day 5/5/19 between 12 and 1:40, taken from https://www.suncalc.org/#/53.3481,-6.2483,11/2019.05.05/11:19/1/3)
% T_aC = 13;  %mean air temp in C to match Cooper1985 plot
% P = 850;  %mean solar radiation in W/m^2 to match Cooper1985 plot
%v = 3.1; %flight speed or wind speed (3.1m/s from Cooper1985)
T_gC = 11;                  %ground surface temp in C, to match Cooper1985
% v = 0;  %no wind for a shivering bee
%T_gC = 11;                  %ground surface temp in C, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
T_aK = T_aC+273.15;     %air temp in K
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
kappa = 0.024;   %ish for 10-15C     %thermal conductivity of air (fill in an equation )
nu = 2.791*10^(-7)*T_aK^(0.7355)/1.225;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%% Bee Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
indicator = 3; %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)
%A_b = 2.46*10^(-5)+4.15*10^(-5)+9*10^(-5);     %sum of h, th, ab reported in Cooper1985
A_th = 9.218*10^(-5)  ; %thorax surface area in m^2, from Church1960
A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
M_b = 0.149;   %mass of the bee in g, Joos1991
M_th = 0.057;  %mass of thorax in g, Joos1991
M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 4.184*0.8;  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
%v = 0;
v = 3.1;
%alpha_si = 0.25*cscd(beta);     %shape factor for incoming solar radiation (Cooper1985), cscd is csc in degrees
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
 C_l = 2.43*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
 n = 1.98;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
E = 0.63;    %Brown2004 activation energy
delta_T_h = 3;

I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
%I_flying = 0.2097035;     %Heinrich1975, converted to W
I_flying = 0.06229515;     %Kammer1974, converted to W
%I_flying = 0.0025;     %temporarily playing around with value of I_flying to see effects

%masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
%masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)

%RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)

%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
%r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
r = 0.008;    %playing around with value to get desired behaviour

%%%%%% Solar Radiation %%%%%%%
% Does not depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation thorax
S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

%%%%%% Thermal Radiation %%%%%%%
% Term 2 does depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant

R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee thorax (from sky + earth)
R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee thorax
R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head

%%%%%% Evaporation %%%%%%%%%%%
%Probably needs to depend on T_h when used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ev = 0/(M_th*c);    %not included for now
%Ev_h = 0/(M_th*c);    %not included for now

%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_l = 0.0749;   %(fill in reference) 
% n = 0.78;       %(fill in reference)
h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;  %=1.2382 at 2021/11/09 defaults
% h = 1.372341;   %calculated according to fromula in heat transfer book and data in Church1960 (R code)
%h = 10;
%C1 = (h*A_th)/(M_th*c);           %thorax, will be multiplied by T_th
C1 = (h*A_th*0.9)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(indicator);
M_ref_th = 0.06;    %thorax weight of bees from Kammer1974
T_ref = RefTemps(indicator);
% i0_resting = exp(log(I_resting) - (3/4)*log(M_ref) + E/(k*(T_ref)));  %fit i_0 to the data (I_resting has units W)
% i0_flying = exp(log(I_flying) - (3/4)*log(M_ref) + E/(k*(T_ref))); %fit i_0 to the data (I_flying has units W)
% norm_constants = [i0_resting, i0_flying, i0_flying];   %resting/thermoregulating/flying = 1,2,3

%indicator = randsample(3,length(StartOrdered),true); %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)
norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3
i_0 = norm_constants(indicator);

%I = (i_0*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)*(1/(M_th*c)); %uses i0 fit, depends on T_air, ref temp/mass in i0 fit
%I = i_0 *  ((M_b/M_ref)^(3/4)) * exp( (-E/k)*((1/T_aK -(1/T_ref)) )) * (1/(M_th*c)); %uses i0=I, depends on T_air
I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass
%I = (i_0*(M_b)^(3/4))*(1/(M_th*c)); %uses i0 fit, depends on T_th, ref temp/mass in i0 fit
%I = (i_0*((M_th/M_ref_th)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass


%%%%%%%%%% Heat transfer to rest of body %%%%%
%Probably needs to depend on T_h when used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frac_transfer = 0.16;    %fraction of metabolic heat generated that can be transferred to the rest of the body

%T = (frac_transfer*I)/(M_th*c);    %heat transferred

Ab = r*(1/(M_th*c));      %heatsink and Tairs version 
%Ab = r*(1-(T_aK/323.15))*(1/(M_th*c));   %bothways version
%%%%%%%%%%%%%%% Calculate just the heatflux %%%%%%%%%%%%%%%%%
% T_thK = 39+273.15;
% Q_default = (S + R1 - R2*T_thK^4 - C1*T_thK - C2 + I)*M_th*c + (S_h + R1_h - R2_h.*(T_thK-delta_T_h).^4 - C1_h.*(T_thK-delta_T_h) - C2_h)*M_th*c; 


%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed

maxy=70+273.15;
%u0=[vf1m(jj) 0]';
%tspan=[0 13]';
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));
%[t,u]= ode45(@(t,u) ok(t,u,p2,m,maxu,maxT), tspan, u0,Opt);

%tiledlayout(2,2);  %to plot multiple panes, run this once before plotting
%any
%nexttile;  %to plot multiple panes, this has to go before each plot
%tspan = [min(StartOrdered),max(StartOrdered)];  %the time over which I have data
tspan = [0,1000];  %the time over which I have data
%[t,y] = ode45(@(t,y) heatflux(t,y,St,S,R1t,R1,R2,C1,C2t,C2,It,I),tspan,y0);
y0 = 39+273.15;   %initial temperature of the bee's thorax in K 
[t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0, Opt); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,Ab,T_aK),tspan,y0); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, Opt); %for use with constant environmental conditions
plot(t,y-273.15,'k');   %plot in celsius  'r' makes a red line; 'b' = blue; 'k' = black
%title('k=5');
subtitle(['P=',num2str(P),'W/m^2',', T_a=',num2str(T_aC),'C'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
yline(30, 'b');  %30 is the min thorax temp for flight, Heinrich1983
yline(50, 'r');   %lethal thorax temp Heinrich1976?ish

% min(y)-273.15
% max(y)-273.15

% n_iterations = length(y);
% last_its = n_iterations-20;
% mean(y(last_its:n_iterations))-273.15
% 


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

