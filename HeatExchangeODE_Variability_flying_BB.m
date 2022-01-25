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

Parameter_Values = readtable('ParameterSample_BB_1000.csv','ReadVariableNames',true);
%row = sample; column = parameter (1st column is row numbers...)
n_samples = 2000;
BB_Thorax_Equilibria_Variability = zeros(n_samples,1);
%heatflux_flying_BB = zeros(n_samples,1);

% mins = zeros(n_samples,1);
% maxs = zeros(n_samples,1);
% lasts = zeros(n_samples,1);

%%%%%%%%%% Constant Values %%%%%%%%%%%
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
k = 8.617333262145*10^(-5);   %Bolzmann's constant
masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
indicator = 3;  %is the bee resting/shivering/flying
cutoffTemp = 100;  %upper limit temp to cut numerical solving at

%prepare the plot for ploting solution curves and histogram 
tplot= tiledlayout(1,1);   %Create a 1-by-1 tiled chart layout tplot.
ax1 = axes(tplot);   %Create an axes object ax1 by calling the axes function and specifying tplot as the parent object.

for sample_index=1:n_samples
%%%%% General Bee Info %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_T_h = Parameter_Values.V1(sample_index);
I_resting = Parameter_Values.V2(sample_index);  %note to check T_ref if using Heinrich vs Kammer data
I_flying = Parameter_Values.V3(sample_index);  %note to check T_ref if using Heinrich vs Kammer data
M_b = Parameter_Values.V4(sample_index);   %mass of the bee in g (the default used by Cooper1985 is 100mg)
E = Parameter_Values.V5(sample_index);    %Brown2004 activation energy
M_th = Parameter_Values.V6(sample_index);  %mass of thorax in g (mean reported in Cooper1985)
c = Parameter_Values.V7(sample_index);  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
r = Parameter_Values.V8(sample_index);
T_mK = Parameter_Values.V9(sample_index);  %midpoint of abdomen function
alpha_si = Parameter_Values.V10(sample_index);     %shape factor for incoming solar radiation (Cooper1985)
epsilon_a = Parameter_Values.V11(sample_index);   % absorptivity of bees (Willmer1981)
A_th = Parameter_Values.V12(sample_index)  ; %in m^2, from Cooper1985
A_h = Parameter_Values.V13(sample_index)  ; %head surface area in m^2, from Cooper1985
alpha_so = Parameter_Values.V14(sample_index);     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = Parameter_Values.V15(sample_index);     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
a = Parameter_Values.V16(sample_index);   %called f in paper, fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)
P = Parameter_Values.V17(sample_index); %mean solar irradiance 
T_gC = Parameter_Values.V18(sample_index);                  %ground surface temp in C
%T_aC = mean(ArrianDataSet.meantempstation);    %mean air temp in C from Arrian's data
T_aC = Parameter_Values.V19(sample_index);    %mean air temp in C from Arrian's data
epsilon_e = Parameter_Values.V20(sample_index);       %(fill in the reference for this!)
f_s = Parameter_Values.V21(sample_index);   %fraction of core temp at surface
C_l = Parameter_Values.V22(sample_index);   %(fill in reference) 
n = Parameter_Values.V23(sample_index);       %(fill in reference)
l_th = Parameter_Values.V24(sample_index);   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
v_flying = Parameter_Values.V25(sample_index); %flight speed or wind speed (3.1m/s from Cooper1985)
%A_b = 2.46*10^(-5)+4.15*10^(-5)+9*10^(-5);     %sum of h, th, ab reported in Cooper1985
%alpha_si = Parameter_Values.V1(sample_index)*cscd(beta);     %shape factor for incoming solar radiation (Cooper1985)
%M_h = Parameter_Values.V20(sample_index);  %mass of head in g (mean reported in Cooper1985)

v_options = [0,0,v_flying];
v = v_options(indicator);

%%%%% General Environment Info %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%T_gC = 11;                  %ground surface temp in C, very vague estimate %from https://www.met.ie/forecasts/farming/agricultural-data-report for May
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
%T_aC = 16.2549375; %morning air temp from my BB data
%P = 599.996875; %morning solar radiation from my BB data
%P = 332.3878; %mean solar irradiance from all of Arrian's data
%beta = (45.47+56.06)/2; %average sun angle in morning from my BB data
%beta = 50.98; %angle of the sun on average (50.93 for Dublin, day 5/5/19 between 11:53 and 1:40, taken from https://www.suncalc.org/#/53.3481,-6.2483,11/2019.05.05/11:19/1/3)
% T_aC = 13;  %mean air temp in C to match Cooper1985 plot
% P = 850;  %mean solar radiation in W/m^2 to match Cooper1985 plot
%indicator = 3; %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)
%T_thK = 39+273.15;   %initial temperature of the bee in K
T_aK = T_aC+273.15;     %air temp in K
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
%a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)
%nu = 2.791*10^(-7)*T_aK^(0.7355)/1.2256;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
rho = 1.2256;  %in dry air air at sea level at 15C
nu = mu/rho; %kinematic viscosity %The Shock Absorber Handbook
kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));



%%%%%% Solar Radiation %%%%%%%
% Does not depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation thorax
S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

%%%%%% Thermal Radiation %%%%%%%
% Term 2 does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee thorax (from sky + earth)
R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee thorax
R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head

% %%%%%% Evaporation %%%%%%%%%%%
% %Probably needs to depend on T_th when used
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ev = 0/(M_th*c);    %not included for now
% Ev_h = 0/(M_th*c);    %not included for now

%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;
C1 = (h*A_th*f_s)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
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

%I = (i_0*(M_b^(3/4))*exp((-E/(k*T_aK)))*M_b)*(1/(M_th*c)); %uses i0 fit, depends on T_air, ref temp/mass in i0 fit
%I = i_0 *  ((M_b/M_ref)^(3/4)) * exp( (-E/k)*((1/T_aK -(1/T_ref)) )) * (1/(M_th*c)); %uses i0=I, depends on T_air
I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass
%I = (i_0*(M_b)^(3/4))*(1/(M_th*c)); %uses i0 fit, depends on T_th, ref temp/mass in i0 fit
%I = (i_0*((M_th/M_ref_th)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass


%%%%%%%%%% Transfer to rest of body %%%%%%%%%%%%%%%
Ab = r*(1/(M_th*c));   %for heatsink version and both ways version
%Ab = r*(1/(M_th*c))*9;      %heatsink and Tairs version %suppose for a minute, T_th-T_ab=9 always
%Ab = r*(1-(T_aK/323.15))*(1/(M_th*c)); %for Tairs version 



% %%%%%%%%%%%%% Calculate just the Heat Flux %%%%%%%%%%%%%%%%%%%%%
% heatflux_flying_BB(sample_index) = (S + R1 - R2*T_thK^4 - C1*T_thK - C2 + I)*M_th*c + (S_h + R1_h - R2_h.*(T_thK-delta_T_h).^4 - C1_h.*(T_thK-delta_T_h) - C2_h)*M_th*c; 
% %have to multiply by M_th*c to undo the division in the previous calculations
% writematrix(heatflux_flying_BB,'heatflux_flying_BB.csv');



%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed

%the times at which the data for the time dependent funcions was measured
%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed
maxy=cutoffTemp+273.15;
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode45 stop solving at T_th=maxy


pcolors=['b' 'y' 'r'];
%tiledlayout(2,2);  %to plot multiple panes, run this once before plotting
%any
%nexttile;  %to plot multiple panes, this has to go before each plot
%tspan = [0,1000,1:1000]; 
tspan = 0:2000;  %the time over which I have data
y0 = Parameter_Values.V26(sample_index)+273.15;   %initial temperature of the bee's head in K (observed by me! ish)
%y0=30+273.14;

[t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK), tspan, y0, Opt); %for use with constant environmental conditions
plot(t,y-273.15,'color',[0,0,0]+0.5);   %plot in celsius  'r' makes a red line; 'b' = blue; 'k' = black
hold on;

if length(y)<length(tspan)
    y((length(y)+1):length(tspan))=cutoffTemp+273.15;   %if it got cut off with T_th>50, fill in the rest of y
end

if length(y)<length(tspan)
    y((length(y)+1):length(tspan))=cutoffTemp+273.15;   %if it got cut off with T_th>50, fill in the rest of y
end

% n_iterations = length(y);
% last_its = n_iterations-20;
if max(y) == cutoffTemp+273.15   %if temp reached cutoff temp, use that
    BB_Thorax_Equilibria_Variability(sample_index) = cutoffTemp;
else  %otherwise, take the mean/median of the end
BB_Thorax_Equilibria_Variability(sample_index) = median(y(1000:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
end

% mins(sample_index) = min(y)-273.15;
% maxs(sample_index) = max(y)-273.15;
% lasts(sample_index) = y(n_iterations)-273.15;
sample_index
end
if indicator==3
    writematrix(BB_Thorax_Equilibria_Variability,'BB_Thorax_Equilibria_Variability_flying.csv');
end
if indicator==1
    writematrix(BB_Thorax_Equilibria_Variability,'BB_Thorax_Equilibria_Variability_resting.csv');
end

%histogram(equilibria_flying_BB-273.15,'orientation','horizontal'); %I like
%how this is transparent

% [counts,bins] = hist(equilibria_flying_BB-273.15); %# get counts and bin locations
% barh(bins,counts,'hist')   %this does put the histogram sideways, but on the left axis
% h=barh(bins,counts); %# include previous two lines from above
% set(get(h,'Parent'),'xdir','r')  %This turns the whole plot around, not just the histogram

%add the labels to the solution curve plots
title('Flying Bumblebee');
subtitle(['P=',num2str(P),'W/m^2'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(45);   %lethal thorax temp Heinrich1976?ish
ylim([0,cutoffTemp])    %set the y-axis limits to match the rotated histogram



%%%%% This is the actual version to use! %%%%%%
ax2 = axes(tplot);  %Create an axes object ax2 by calling the axes function and specifying t as the parent object.
h = histogram(BB_Thorax_Equilibria_Variability,'orientation','horizontal');  %establish the histogram plot
set(get(h,'Parent'),'xdir','r')   %put it on the other axis
ax2.XAxisLocation = 'top'; %Move the x-axis to the top,  
ax2.YAxisLocation = 'right';   %and move the y-axis to the right.
ylim([0,cutoffTemp])    %set the y-axis limits for the rotated histogram
set(ax2(1),'XTickLabel','');   %but don't actually show the axis labels
set(ax2(1),'YTickLabel','');
xlim([0,2000])  %reset the x-limits (for the histogram only), for visual effect
ax2.Color = 'none';  %Set the color of the axes object to 'none' so that the underlying plot is visible.
ax1.Box = 'off';    %Turn off the plot boxes to prevent the box edges from obscuring the x- and y-axes. 
ax2.Box = 'off';