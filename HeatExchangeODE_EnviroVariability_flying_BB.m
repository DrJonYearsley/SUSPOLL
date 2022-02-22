%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ArrianDataSet = readtable('HiveActivityWeatherdataset.csv','ReadVariableNames',true);
% WeatherDateIndices = find(ismember(ArrianDataSet.species, 'H') & ArrianDataSet.site==2 & ArrianDataSet.date=='28/05/2019'); 
%use weather data for site 1 on April 26th 2019 (or whenever)
% ShortDataSet = ArrianDataSet(WeatherDateIndices,:);  %: indicates all colums

%time of each observation in minutes from midnight
% StartH_string = cellfun(@(x) x(1:2),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting hour
% StartM_string = cellfun(@(x) x(5:6),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting minutes
% StartH = cellfun(@str2double,StartH_string);
% StartM = cellfun(@str2double,StartM_string);
% Start = 60*StartH+StartM;         %convert hours and minutes into minutes

%Keep just the entries with unique times
% [StartOrdered,StartIndex] = unique(Start); 
% ShortDataSetOrdered = ShortDataSet(StartIndex,:); 

Enviro_Values_flying = readtable('ParameterSample_Enviro_1000.csv','ReadVariableNames',true);
%row = sample; column = parameter (1st column is row numbers...)
n_samples = 1000;
equilibria_flying_BB_enviro = zeros(n_samples,1);
%heatflux_flying_BB_enviro = zeros(n_samples,1);

% mins = zeros(n_samples,1);
% maxs = zeros(n_samples,1);
% lasts = zeros(n_samples,1);

%%%%%%%%%% Fixed Constants %%%%%%%%%%%%%%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant

%%%%% Bee Inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%E = 0.63;    %Brown2004 activation energy
E = 0; 

indicator = 3; %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)
T_thK = 39+273.15;   %initial temperature of the bee in K
%A_b = 2.46*10^(-5)+4.15*10^(-5)+9*10^(-5);     %sum of h, th, ab reported in Cooper1985
A_th = 9.218*10^(-5)  ; %in m^2, from 
A_h = 2.46*10^(-5)  ; %HB head surface area in m^2, from Cooper1985
M_b = 0.149;   %mass of the bee in g 
M_th = 0.057;  %mass of thorax in g (mean reported in Cooper1985)
% M_h = Enviro_Values_flying.V20(sample_index);  %mass of head in g (mean reported in Cooper1985)
l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 4.184*0.8;  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
epsilon_a = 0.935;   % absorptivity of bees (Willmer1981)
v_options = [0 0 3.1];   %make the bee be out of wind when resting/shivering
v = v_options(indicator);
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
C_l = 2.429809*10^(-7);   %(fill in reference) 
n = 1.975485;       %(fill in reference)

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

%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour

% prepare the plot for ploting solution curves and histogram 
tplot= tiledlayout(1,1);   %Create a 1-by-1 tiled chart layout tplot.
ax1 = axes(tplot);   %Create an axes object ax1 by calling the axes function and specifying tplot as the parent object.

for sample_index=1:n_samples
%%%%%%%%%%% Environmental Inputs %%%%%%%%%%%%%%%%
%T_aC = 16.2549375; %morning air temp from my BB data
%P = 599.996875; %morning solar radiation from my BB data
P = Enviro_Values_flying.V1(sample_index); %mean solar irradiance from Arrian's data
%T_aC = Enviro_Values_flying.V2(sample_index);    %mean air temp in C from Arrian's data
%beta = (45.47+56.06)/2; %average sun angle in morning from my BB data
%beta = Enviro_Values_flying.V4(sample_index); %angle of the sun on average (50.93 for Dublin, day 5/5/19 between 11:53 and 1:40, taken from https://www.suncalc.org/#/53.3481,-6.2483,11/2019.05.05/11:19/1/3)
% beta = 48.52; %angle of the sun on average (??? for Dublin, day 5/5/19 between 12 and 1:40, taken from https://www.suncalc.org/#/53.3481,-6.2483,11/2019.05.05/11:19/1/3)
% T_aC = 13;  %mean air temp in C to match Cooper1985 plot
% P = 850;  %mean solar radiation in W/m^2 to match Cooper1985 plot
%T_gC = 11;                  %ground surface temp in C, very vague estimate %from https://www.met.ie/forecasts/farming/agricultural-data-report for May
T_aK = Enviro_Values_flying.V2(sample_index);     %air temp in K
T_gK = Enviro_Values_flying.V3(sample_index);        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
kappa = Enviro_Values_flying.V6(sample_index);   %ish for 10-15C     %thermal conductivity of air (fill in an equation )
nu = Enviro_Values_flying.V7(sample_index);   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = Enviro_Values_flying.V5(sample_index);   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)
% kappa = 0.024;   %ish for 10-15C     %thermal conductivity of air (fill in an equation )
% nu = 2.791*10^(-7)*T_aK^(0.7355)/1.225;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
% a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)


%%%%%% Solar Radiation %%%%%%%
% Does not depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation
S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

%%%%%% Thermal Radiation %%%%%%%
% Term 2 does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee (from sky + earth
R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee
R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head

%%%%%% Evaporation %%%%%%%%%%%
%Probably needs to depend on T_th when used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ev = 0/(M_th*c);    %not included for now
% Ev_h = 0/(M_th*c);    %not included for now

%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;

C1 = (h*A_th*0.9)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(indicator);
T_ref = RefTemps(indicator);
norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3
i_0 = norm_constants(indicator);

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
maxy=70+273.15;
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode45 stop solving at T_th=maxy

pcolors=['b' 'y' 'r'];

tspan = 0:1000;  %the time over which I have data
%[t,y] = ode45(@(t,y) heatflux(t,y,St,S,R1t,R1,R2,C1,C2t,C2,It,I),tspan,y0);
y0 = 39+273.15;   %initial temperature of the bee in K

[t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0,Opt); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_alwaysactive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK), tspan, y0); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_c_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,Ab,T_aK),tspan,y0,Opt); %for use with constant environmental conditions
%[t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, Opt); %for use with constant environmental conditions
plot(t,y-273.15,'color',[0,0,0]+.5);   %plot in celsius 'k:' creates a black dotted line
hold on;

if length(y)<1001
    y((length(y)+1):1001)=70+273.15;   %if it got cut off with T_th>70, fill in the rest of y
end

equilibria_flying_BB_enviro(sample_index) = median(y(900:1001))-273.15;     %median instead of mean should be more stable against cycles
% mins(sample_index) = min(y)-273.15;
% maxs(sample_index) = max(y)-273.15;
% lasts(sample_index) = y(n_iterations)-273.15;
sample_index   %will show how far we got

end
writematrix(equilibria_flying_BB_enviro,'equilibria_flying_BB_enviro.csv');

%add the labels to the solution curve plots
title('Flying Bumblebee');
subtitle(['P=',num2str(P),'W/m^2'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(50);   %lethal thorax temp Heinrich1976?ish
yline(43,'k:');
ylim([30,70])    %set the y-axis limits to match the rotated histogram



%%%%% This is the actual version to use! %%%%%%
ax2 = axes(tplot);  %Create an axes object ax2 by calling the axes function and specifying t as the parent object.
h = histogram(equilibria_flying_BB_enviro,'orientation','horizontal');  %establish the histogram plot
set(get(h,'Parent'),'xdir','r')   %put it on the other axis
ax2.XAxisLocation = 'top'; %Move the x-axis to the top,  
ax2.YAxisLocation = 'right';   %and move the y-axis to the right.
ylim([30,70])    %set the y-axis limits for the rotated histogram
set(ax2(1),'XTickLabel','');   %but don't actually show the axis labels
set(ax2(1),'YTickLabel','');
xlim([0,1000])  %reset the x-limits (for the histogram only), for visual effect
ax2.Color = 'none';  %Set the color of the axes object to 'none' so that the underlying plot is visible.
ax1.Box = 'off';    %Turn off the plot boxes to prevent the box edges from obscuring the x- and y-axes. 
ax2.Box = 'off';


% figure
% histogram(equilibria_flying_BB_enviro);
% xlim([30,70])
%xlabel('Equlibrium Thorax Temperature (C)') ;