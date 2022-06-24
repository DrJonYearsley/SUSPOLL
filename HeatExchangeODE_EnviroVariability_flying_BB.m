%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
Bumblebee = false;    %only one of Bumblebee and Honeybee can be true
Honeybee = true;
CoolingOff = false;   %only one of cooling can be true
CoolingSwitch = true;  
CoolingOn = false; 
Resting = false;      %only one of metabolic states can be true 
Shivering = false; 
Flying = true;
indicator = 3;  %is the bee resting/shivering/flying

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Code Structure Stuff %%%%%%%%
n_samples = 1000;
cutoffTemp = 100;  %upper limit temp to cut numerical solving at
maxy=cutoffTemp+273.15;
pcolors=['b' 'y' 'r'];
tspan = 0:2000;  %the time over which I have data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Storage Vectors %%%%%%%%%%%%
Thorax_Equilibria_Variability = zeros(n_samples,1);



%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
h_fg = 2.3819*10^6;  %latent heat of vaporization of water at 50C, J/kg  (it actually depends on temp)
%https://www.engineeringtoolbox.com/water-properties-d_1573.html
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
%D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s - %calculated version is below
MM_air = 0.0289652;   %molar mass of dry air in kg/mol
%https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol
R_specific = 287.058;  %J/kg/K for dry air
Pr = 1.013*10^5; %atmospheric pressure in N/m^2 (value used in Sidebotham)





%%%%%%%%%%% Bee Parameters same for BB and HB %%%%%%%%%%%%%%%%
s = 0.9965;  %ratio calculated from Church1960 data converted to K
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
maxy=70+273.15;    %don't solve above 60C because it's not biologically relevant
tspan = 0:2000; 
r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue
   

if Bumblebee==true
   %%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%

    A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
    A_h = 3.61375*10^(-5)  ; %head surface area in m^2, own data
    M_b = 0.149;   %mass of the bee in g, Joos1991, default
    %M_b = 0.035;   %mass of the bee in g, Joos1991, minimum
    %M_b = 0.351;   %mass of bee in g, Joos1991, maximum
    M_th = 0.057;  %mass of thorax in g, Joos1991
    %M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
    l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0 0.1 4.1];   %make the bee be out of wind when resting/shivering, default
    %v_options = [0 0 5.5];   %make the bee be out of wind when resting/shivering, maximum (Osborne2013)
    %v_options = [0 0 1];   %make the bee be out of wind when resting/shivering, minimum (Osborne2013)
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 42+273.15;     %median temp for  abdomen cooling
    I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
    I_flying = 0.06229515;     %Kammer1974, converted to W
    %I_flying = 0.004020538;     %fitted value
    %I_flying = 0.2097035;     %Heinrich1975, converted to W
    %I_flying = 0.03;   %experiment with the value
    %masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
    %masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
    masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    %RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
    norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3
    E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
    %E_opts = [0.63 0.63 0]; %fitted values

    y0 = 30+273.15;   %initial temperature of the bee's thorax in K 
    LethalTemp = 45;  %lethal thorax/air temp in C
    CoolingTemp = 42;   %thorax temp where cooling begins
    FlyingTemp = 30;      %thorax temp where flight can begin
end

if Honeybee==true
    A_th = 4.5*10^(-5)  ; %thorax surface area in m^2, from ???
    A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 
    M_b = 0.100;   %mass of the bee in g, Joos1991
    M_th = 0.0407;  %mass of thorax in g, Joos1991
    l_th = 0.004;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.91;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0.1 0.1 5.6];   %make the bee be out of wind when resting/shivering, default
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 47.9+273.15;     %median temp for  evaporative cooling

    I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
    I_flying = 0.4*(80/1000); %Nachtigall1989, W/g -> W, 80mg reference mass%    
    %I_flying = 0.005;     %%fitted value
    %I_flying = 0.004;     %experiment with value
    masses = [0.08 0.08 0.08];   %reference weight for Rothe/Nachtigal (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp 
    norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3
    E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
    %E_opts = [0.63 0.63 0]; %fitted values

    y0 = 39+273.15;   %initial temperature of the bee's thorax in K (fill in ref)
    LethalTemp = 52;  %lethal thorax/air temp in C
    CoolingTemp = 47.9;   %thorax temp where cooling begins
    FlyingTemp = 35;      %thorax temp where flight can begin


end

%prepare the plot for ploting solution curves and histogram 
tplot= tiledlayout(1,1);   %Create a 1-by-1 tiled chart layout tplot.
ax1 = axes(tplot);   %Create an axes object ax1 by calling the axes function and specifying tplot as the parent object.

for sample_index=1:n_samples
    
    Parameter_Values = readtable('ParameterSample_Enviro_1000.csv','ReadVariableNames',true);

    P = Parameter_Values.V1(sample_index);
    T_aK = Parameter_Values.V2(sample_index);     %air temp in K
    T_gK = Parameter_Values.V3(sample_index);
    rh = Parameter_Values.V4(sample_index);
    a = Parameter_Values.V5(sample_index);
    %nu = 2.791*10^(-7)*T_aK^(0.7355)/Pr;   %kinematic viscosity of air, equation, divided by air pressure - 'standard sea level' condition from Wikipedoa
    mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
    rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
    nu = mu/rho; %kinematic viscosity for humid air %The Shock Absorber Handbook
    kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));


    %%%%%%%%% Bee parameters that depend on metabolic state j
    v = v_options(indicator);
    D_A = v*(2*R_0);    %diffusion coefficient equivalent to advection
    M_ref = masses(indicator);
    T_ref = RefTemps(indicator);
    i_0 = norm_constants(indicator);
    E = E_opts(indicator);


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


%for now, only set for using cooling switching 
Cooling_vals = Cooling_Flux(Bumblebee,r,M_th,c,A,B,T_aK,rh,Pr,MM_air,MM_vapor,R_specific,R_0,D_A,h_fg);
Ab = Cooling_vals(1);
Ev1 = Cooling_vals(2);
Ev2 = Cooling_vals(3);
CoolingSwitch_indicator = 1;
[t,y] = ode23(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0,Opt); %for use with constant environmental conditions
plot(t,y-273.15,'color',[0,0,0]+0.5);   %plot in celsius  'r' makes a red line; 'b' = blue; 'k' = black
hold on;

y = calc_length(tspan,y,maxy);  %check if bee maxed out temp and fill in if necessary

if max(y) == cutoffTemp+273.15   %if temp reached cutoff temp, use that
    Thorax_Equilibria_Variability(sample_index) = cutoffTemp;
else  %otherwise, take the mean/median of the end
    Thorax_Equilibria_Variability(sample_index) = median(y(1500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
end

sample_index  %print to show where we are

end

if Honeybee==true
    if indicator==3
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_flying_1000_HB.csv');
    end
    if indicator==2
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_shivering_1000_HB.csv');
    end
    if indicator==1
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_resting_1000_HB.csv');
    end
end

if Bumblebee==true
    if indicator==3
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_flying_1000_BB.csv');
    end
    if indicator==2
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_shivering_1000_BB.csv');
    end
    if indicator==1
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_EnviroVariability_resting_1000_BB.csv');
    end
end

%add the labels to the solution curve plots
title('Flying Bee');
subtitle(['P=',num2str(P),'W/m^2'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
yline(30);  %30 is the min thorax temp for flight, Heinrich1983
yline(45);   %lethal thorax temp Heinrich1976?ish
ylim([0,cutoffTemp])    %set the y-axis limits to match the rotated histogram



%%%%% Plot a nice sideways histogram on top of it! %%%%%%
ax2 = axes(tplot);  %Create an axes object ax2 by calling the axes function and specifying t as the parent object.
h = histogram(Thorax_Equilibria_Variability,'orientation','horizontal');  %establish the histogram plot
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