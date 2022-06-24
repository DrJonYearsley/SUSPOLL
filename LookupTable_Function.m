function TemperatureTable = LookupTable_Function(Switches,Parameter_Values,tspan)

%%%%%%%%%%% Storage vectors %%%%%%%%%%%%%%%%
Temps_median = zeros(51,10); %rows for temperature, colums for resting/shivering/flying
Temps_median(:,10) = (0:50).'; %fourth column is the temperature


%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
Bumblebee = Switches(1);    %only one of Bumblebee and Honeybee can be true
Honeybee = Switches(2);

CoolingOff = Switches(3);   %any combination of cooling can be true
CoolingSwitch = Switches(4);  
CoolingOn = Switches(5); 

Resting = Switches(6);      %any combination of metabolic states can be true 
Shivering = Switches(7); 
Flying = Switches(8);

%Fitting = Switches(9);     %only true if varying mass


%%%%%% Constant values %%%%%%%%
k = Parameter_Values(1);   %Bolzmann's constant
delta = Parameter_Values(2);   %fill in the name of this constant
sigma = Parameter_Values(3);   %fill in the name of this constant
h_fg = Parameter_Values(4);  %latent heat of vaporization of water at 50C, J/kg  (it actually depends on temp)
A = Parameter_Values(5);  %clausius-clapyron constant A for water in N/m^2
B = Parameter_Values(6);  %clausius-clapyron constant B for water in K
rh =  Parameter_Values(7);   %mean relative humidity from Arrian's data
MM_air = Parameter_Values(8);   %molar mass of dry air in kg/mol
MM_vapor = Parameter_Values(9);  %molar mass of water vapor in kg/mol

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%

P = Parameter_Values(10); %mean solar irradiance from all of Arrian's data
%T_gC = Parameters(11);                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
T_gK = Parameter_Values(12);        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
Pr = Parameter_Values(13); %atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = Parameter_Values(14);  %J/kg/K for dry air
a = Parameter_Values(15);   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%% Bee Parameters same for BB and HB %%%%%%%%%%%%%%%%
s = Parameter_Values(16);  %fraction of internal temp at surface - BB fitted value
c = Parameter_Values(17);  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
C_l = Parameter_Values(18);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = Parameter_Values(19);       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = Parameter_Values(20);
alpha_si = Parameter_Values(21);     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = Parameter_Values(22);     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = Parameter_Values(23);     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
maxy=Parameter_Values(24);    %don't solve above 60C because it's not biologically relevant
r = Parameter_Values(25);  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
R_0 = Parameter_Values(26);  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue
   

if Bumblebee==true
   %%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%

    A_th = Parameter_Values(27)  ; %thorax surface area in m^2, from Church1960
    A_h = Parameter_Values(28)  ; %head surface area in m^2, own data
    M_b = Parameter_Values(29);   %mass of the bee in g, Joos1991, default
    M_th = Parameter_Values(30);  %mass of thorax in g, Joos1991
    l_th = Parameter_Values(31);   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = Parameter_Values(32);   % absorptivity of bees (Willmer1981, ,te)
    v_options = Parameter_Values(33:35);   %make the bee be out of wind when resting/shivering, default
    epsilon_e = Parameter_Values(36);       %(fill in the reference for this!)
    T_mK = Parameter_Values(37);     %median temp for  abdomen cooling
    %I_resting = Parameters(38);     %Kammer1974, table 1, for 25C, converted to W
    %I_flying = Parameters(39);   %experiment with the value
    masses = Parameter_Values(40:42);   %reference weight for Kammer (flying)
    RefTemps = Parameter_Values(43:45);   %Reference temp is 25C for Kammer (resting & flying)
    norm_constants = Parameter_Values(46:48);   %resting/shivering/flying = 1,2,3
    E_opts = Parameter_Values(49:51); %fitted values

    y0 = Parameter_Values(52);   %initial temperature of the bee's thorax in K 
    %LethalTemp = Parameters(53);  %lethal thorax/air temp in C
    %CoolingTemp = Parameters(54);   %thorax temp where cooling begins
    %FlyingTemp = Parameters(55);      %thorax temp where flight can begin
end

if Honeybee==true
    A_th = Parameter_Values(27)  ; %thorax surface area in m^2, from Church1960
    A_h = Parameter_Values(28)  ; %head surface area in m^2, own data
    M_b = Parameter_Values(29);   %mass of the bee in g, Joos1991, default
    M_th = Parameter_Values(30);  %mass of thorax in g, Joos1991
    l_th = Parameter_Values(31);   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = Parameter_Values(32);   % absorptivity of bees (Willmer1981, ,te)
    v_options = Parameter_Values(33:35);   %make the bee be out of wind when resting/shivering, default
    epsilon_e = Parameter_Values(36);       %(fill in the reference for this!)
    T_mK = Parameter_Values(37);     %median temp for  abdomen cooling
    %I_resting = Parameters(38);     %Kammer1974, table 1, for 25C, converted to W
    %I_flying = Parameters(39);   %experiment with the value
    masses = Parameter_Values(40:42);   %reference weight for Kammer (flying)
    RefTemps = Parameter_Values(43:45);   %Reference temp is 25C for Kammer (resting & flying)
    norm_constants = Parameter_Values(46:48);   %resting/shivering/flying = 1,2,3
    E_opts = Parameter_Values(49:51); %fitted values

    y0 = Parameter_Values(52);   %initial temperature of the bee's thorax in K 
    %LethalTemp = Parameters(53);  %lethal thorax/air temp in C
    %CoolingTemp = Parameters(54);   %thorax temp where cooling begins
    %FlyingTemp = Parameters(55);      %thorax temp where flight can begin
end

for i = 1:51  %go through temps 0 to 50 
    for j = nonzeros([Resting Shivering Flying].*[1 2 3])'  %go through resting/shivering/flying as indicated at start
    
    %%%%%%%%%%% Environmental Parameters that depend on air temp %%%%%%%%%%%%%%%%
    T_aC = i-1;    %varying air temp
    T_aK = T_aC+273.15;     %air temp in K
    %nu = 2.791*10^(-7)*T_aK^(0.7355)/Pr;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
    mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
    rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
    nu = mu/rho; %kinematic viscosity for humid air %The Shock Absorber Handbook
    kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));


    %%%%%%%%% Bee parameters that depend on metabolic state j
    v = v_options(j);
    D_A = v*(2*R_0);    %diffusion coefficient equivalent to advection
    M_ref = masses(j);
    T_ref = RefTemps(j);
    i_0 = norm_constants(j);
    E = E_opts(j);


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
    C1 = (h*A_th*s)/(M_th*c);           %thorax, will be multiplied by T_th, *s is for thorax surface temperature
    C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
    C1_h = (h*A_h)/(M_th*c);           %head
    C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



    %%%%%%%%%%% Metabolic %%%%%%%%%
    %Does  depend on T_th 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass

    
     
    
%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatfluxhd.m, which may need to be
%renamed
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode23 stop solving at T_th=maxy


%solve and plot solution curves separately for passive or active model 

if CoolingOff == true %physiological model with abdomen off    
    saving_index = j;   %index for storage vector is 1, 2, or 3
    Ab = 0;
    Ev1 = 0;
    Ev2 = 0;   %no abdomen or evaporative cooling
    CoolingSwitch_indicator = 0;  %no switching
[t,y] = ode23(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0, Opt); %for use with constant environmental conditions
    y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary
    if max(y) == maxy   %if temp reached maximum, use that
        Temps_median(i,saving_index) = maxy-273.15;
    else  %otherwise, take the mean/median of the end
        Temps_median(i,saving_index) = median(y(1500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
    end
    [i j saving_index]; %print out to show where we are
end

%clear y t saving_index

if CoolingSwitch == true %behavioural model with abdomen switching
    saving_index = j+3;   %index for storage vector
    Cooling_vals = Cooling_Flux(Bumblebee,r,M_th,c,A,B,T_aK,rh,Pr,MM_air,MM_vapor,R_specific,R_0,D_A,h_fg,v);
    Ab = Cooling_vals(1);
    Ev1 = Cooling_vals(2);
    Ev2 = Cooling_vals(3);
    CoolingSwitch_indicator = 1;
    [t,y] = ode23(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0,Opt); %for use with constant environmental conditions
    y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary
    if max(y) == maxy   %if temp reached maximum, use that
        Temps_median(i,saving_index) = maxy-273.15;
    else  %otherwise, take the mean/median of the end
        Temps_median(i,saving_index) = median(y(1500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
    end
    [i j saving_index]; %print out to show where we are

end

%clear y t saving_index

if CoolingOn == true %behavioural model with abdomen coolingalways on    
    saving_index = j+6;   %index for storage vector
    Cooling_vals = Cooling_Flux(Bumblebee,r,M_th,c,A,B,T_aK,rh,Pr,MM_air,MM_vapor,R_specific,R_0,D_A,h_fg,v);
    Ab = Cooling_vals(1);
    Ev1 = Cooling_vals(2);
    Ev2 = Cooling_vals(3);
    CoolingSwitch_indicator = 0;  %no switch for cooling
    
    [t,y] = ode23(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0,Opt); %for use with constant environmental conditions
    y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary    
    if max(y) == maxy   %if temp reached maximum, use that
        Temps_median(i,saving_index) = maxy-273.15;
    else  %otherwise, take the mean/median of the end
        Temps_median(i,saving_index) = median(y(1500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
    end
    [i j saving_index]; %print out to show where we are
end
    end
end
%writematrix(Temps_median,'LookupTable_Bumblebee_default.csv');


TemperatureTable = Temps_median;