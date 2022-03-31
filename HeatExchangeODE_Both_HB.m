%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ArrianDataSet = readtable('HiveActivityWeatherdataset.csv','ReadVariableNames',true);
WeatherDateIndices = find(ismember(ArrianDataSet.species, 'H') & ArrianDataSet.site==2 & ArrianDataSet.date=='28/05/2019'); 
%use weather data for site 1 on April 26th 2019 (or whenever)
ShortDataSet = ArrianDataSet(WeatherDateIndices,:);  %: indicates all colums

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameters that might be changed, for easy access  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indicator = 3; %indicator for when the bee is resting/thermoregulating/flying (fill in function of time later)
E = 0.63;    %Brown2004 activation energy
%E = 0;
T_aC = 9;
%T_aC = mean(ArrianDataSet.meantempstation);    %mean air temp in C from Arrian's data
P = mean(ArrianDataSet.meansolarstation); %mean solar irradiance from Arrian's data
I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
I_flying = 0.4*(80/1000); %Nachtigall1989, W/g -> W, 80mg reference mass%I_flying = 0.0025;     %temporarily playing around with value of I_flying to see effects
v_options = [0 0 3.1];   %make the bee be out of wind when resting/shivering, default
%v = 0;
%v = 3.1;
v = v_options(indicator);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %time of each observation in minutes from midnight
% StartH_string = cellfun(@(x) x(1:2),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting hour
% StartM_string = cellfun(@(x) x(5:6),ShortDataSet.timestart,'UniformOutput',false);    %extract the starting minutes
% StartH = cellfun(@str2double,StartH_string);
% StartM = cellfun(@str2double,StartM_string);
% Start = 60*StartH+StartM;         %convert hours and minutes into minutes
% 
% %Keep just the entries with unique times
% [StartOrdered,StartIndex] = unique(Start); 
% ShortDataSetOrdered = ShortDataSet(StartIndex,:); 

%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
Pr = 1.013*10^5; %atmospheric pressure in N/m^2
R_specific = 287.058;  %J/kg/K for dry air
%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%
%T_aC = ShortDataSetOrdered.meantempstation;    %air temp in C from Arrian's data
%P = ShortDataSetOrdered.meansolarstation; %solar irradiance from Arrian's data
%T_aC = 5;    %to try different air temps
%beta = 48.52; %angle of the sun on average (??? for Dublin, day 5/5/19 between 12 and 1:40, taken from https://www.suncalc.org/#/53.3481,-6.2483,11/2019.05.05/11:19/1/3)
% T_aC = 13;  %mean air temp in C to match Cooper1985 plot
% P = 850;  %mean solar radiation in W/m^2 to match Cooper1985 plot
%v = 3.1; %flight speed or wind speed (3.1m/s from Cooper1985)
T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
% v = 0;  %no wind for a shivering bee
%T_gC = 11;                  %ground surface temp in C, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
T_aK = T_aC+273.15;     %air temp in K
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
%kappa = 0.024;   %ish for 10-15C     %thermal conductivity of air (fill in an equation )
%nu = 2.791*10^(-7)*T_aK^(0.7355)/1.225;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
rho = Pr/(R_specific*T_aK);  %in dry air air at sea level at 15C - find formula for humid air?
nu = mu/rho; %kinematic viscosity %The Shock Absorber Handbook
kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%% Bee Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
A_th = 4.5*10^(-5)  ; %thorax surface area in m^2, from ???
A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 
M_b = 0.100;   %mass of the bee in g, Joos1991
M_th = 0.0407;  %mass of thorax in g, Joos1991
%M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
l_th = 0.004;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
epsilon_a = 0.91;   % absorptivity of bees (Willmer1981, ,te)
%alpha_si = 0.25*cscd(beta);     %shape factor for incoming solar radiation (Cooper1985), cscd is csc in degrees
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
T_mK = 44+273.15;     %median temp for  evaporative cooling
f = 0.9;  %fraction of internal temp at surface
a_n = ((0.23+0.19+0.44+0.6)/(4*1000))*pi*(4.18/1000);  %surface area of BB tongue, mm->m
 
masses = [0.08 0.08 0.08];   %reference weight for Rothe/Nachtigal (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp 



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


%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;
C1 = (h*A_th*f)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%% Evaporation %%%%%%%%%%%
%Depends on T_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fg = 2.33*10^6;  %latent heat of vaporization of water, J/kg 
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
rh =  mean(ArrianDataSet.meanrhstation)/100;   %relative humidity - actually calculate this from Arrian's data
Pr = 1.013*10^5; %atmospheric pressure in N/m^2
R_0 = 0.001;  %radius of nectar droplet, 1mm in m
D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s
MM_air = 0.0289652;   %molar mass of dry air in g/mol
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol

p_sat_air = A*exp(B/T_aK);  %partial pressure of saturated vapor (humid air at saturation)
p_air = rh*p_sat_air;   %partial pressure of humid air at temperature T_aK
rho_humid = ((Pr-p_air)*MM_air + p_air*MM_vapor)/(R_specific*T_aK);   %density of humid air

X_air = p_air/Pr;  %mole fraction of ambient air
Y_air = 1/( 1 + ((1-X_air)/X_air)*(MM_air/MM_vapor) );

m_evap_coef = 2*pi*R_0*rho_humid*D_A;
Ev1 = h_fg*m_evap_coef*log(1-Y_air);   %constant for a fixed T_air
Ev2 = -h_fg*m_evap_coef;   %*log(1-Y_sfc)

%Ev = Ev1-Ev2*log(1-Y_sfc)

%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(indicator);
T_ref = RefTemps(indicator);

norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3
i_0 = norm_constants(indicator);

I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass




%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed

maxy=100+273.15;
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));
pcolors=['b' 'y' 'r'];

%tspan = 0:2000;  %the time over which I have data
tspan = 0:17;  %debugging for solving up to time at which imaginary solution happens
y0 = 39+273.15;   %initial temperature of the bee's thorax in K 
[t,y] = ode23(@(t,y) heatfluxhead_Tth_active_HB(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ev1,Ev2,rh,A,B,MM_air,MM_vapor,Pr,T_mK), tspan, y0, Opt); %for use with constant environmental conditions
plot(t,y-273.15,'color',pcolors(indicator));   %plot in celsius  'r' makes a red line; 'b' = blue; 'k' = black
yline(30, 'k');  %30 is the min thorax temp for flight, Heinrich1983
yline(45, 'k');   %lethal thorax temp Heinrich1976?ish

if length(y)<length(tspan)
    y((length(y)+1):length(tspan))=maxy;   %if it got cut off with T_th>50, fill in the rest of y
end
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

