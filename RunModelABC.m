function KeyTemps = RunModelABC(i0_guess, r_guess)
Temps = zeros(51,4); %rows for temperature, colums for resting/shivering/flying
Temps(:,1) = (0:50).'; %first column is the temperature

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
E = 0.63;    %Brown2004 activation energy
%E = 0; 

A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
M_b = 0.149;   %mass of the bee in g, Joos1991
M_th = 0.057;  %mass of thorax in g, Joos1991
M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
c = 3.349;  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
v_options = [0 0 4.1];   %make the bee be out of wind when resting/shivering
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = 0.97;       %(fill in the reference for this!)
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
T_mK = 42+273.15;     %median temp for  abdomen cooling
f = 0.9;  %fraction of internal temp at surface

masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
%i0 will be provided by fitting algorithm

%r will be provided by fitting algorithm

for i = 1:51  %go through temps 0 to 50 
    for j = 1:3    %go through the states (shivering-active; flying-passive; flying-active)
% j = 3;     %temporarily do just flying       

metabolic_state=[2 3 3]; %go through each metabolic state
model_state=[2 1 2]; %1 for passive, 2 for active

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
C1 = (h*A_th*f)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head



%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(metabolic_indicator);
T_ref = RefTemps(metabolic_indicator);


I = (i0_guess*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass

%%%%%%%%%% Transfer to rest of body %%%%%%%%%%%%%%%
Ab = r_guess*(1/(M_th*c));   %for heatsink version and both ways version


%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed
maxy=70+273.15;
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy));   %using this will make ode45 stop solving at T_th=maxy


%pcolors=['b' 'y' 'r'];
tspan = 0:1000;  %the time over which I have data
y0 = 39+273.15;   %initial temperature of the bee's head in K (observed by me! ish)

%solve and plot solution curves separately for passive or active model 

if model_indicator == 1 %passive model
%    figure(1)
    [t,y] = ode45(@(t,y) heatfluxhead_Tth_passive(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k), tspan, y0, Opt); %for use with constant environmental conditions
%     plot(t,y-273.15,'color',pcolors(metabolic_indicator));   %plot in celsius 'k:' creates a black dotted line
else  %active model
%    figure(2)
    [t,y] = ode45(@(t,y) heatfluxhead_Tth_active(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,Ab,T_aK,T_mK), tspan, y0,Opt); %for use with constant environmental conditions
%     plot(t,y-273.15,'color',pcolors(metabolic_indicator));   %plot in celsius 'k:' creates a black dotted line
end

%hold on; 

if length(y)<1001
    y((length(y)+1):1001)=70+273.15;   %if it got cut off with T_th>70, fill in the rest of y
end

Temps(i,j+1) = median(y(850:1001))-273.15;     %median instead of mean should be more stable against cycles
%[i j]

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify when shivering bee goes above 30C
shiver30 = find(Temps(:,2) >= 30);  %get the array indces where shivering bee temp>=30
if isempty(shiver30) %if the bee never gets above 30, set temp to 50 as max
    shiver30Temp = 50;
else
    shiver30Temp = Temps(shiver30(1),1);   %otherwise, get the air temp for the first index in shiver30
end

%identify when passive and active flying bees diverge
paDifference = Temps(:,3)-Temps(:,4);   %difference between active and passive thorax temps
paDiverge = find(paDifference > 0.5);  %they are always a little bit different, so call it diverging when the difference is >0.25
if isempty(paDiverge)
    paDivergeTemp = 50;  %if they never diverge, set to 50 as maximum
else
    paDivergeTemp = Temps(paDiverge(1),1);  %if they do, get the air temp for the first index in paDiverge
end

%identify when active flying bee goes above 42C
fly42 = find(Temps(:,4) >= 42);  %get the array indces where shivering bee temp>=42
if isempty(fly42)
    fly42Temp = 50;  %if it never goes above 42C, set temp to 50 as max
else
    fly42Temp = Temps(fly42(1),1);   %otherwise, get the air temp for the first index in fly42
end


KeyTemps = [shiver30Temp paDivergeTemp fly42Temp];



% %add the labels to the solution curve plots
% figure(1)
% title('Passive Model');
% subtitle(['P=',num2str(P),'W/m^2'])
% xlabel('Time (s)') ;
% ylabel('Thorax Temperature (C)') ;
% yline(30);  %30 is the min thorax temp for flight, Heinrich1983
% yline(50);   %lethal thorax temp Heinrich1976?ish
% yline(43,'k:');
% legend('resting','shivering','flying','Location','southeast')
% 
% figure(2)
% title('Active Model');
% subtitle(['P=',num2str(P),'W/m^2'])
% xlabel('Time (s)') ;
% ylabel('Thorax Temperature (C)') ;
% yline(30);  %30 is the min thorax temp for flight, Heinrich1983
% yline(50);   %lethal thorax temp Heinrich1976?ish
% yline(43,'k:');
% legend('resting','shivering','flying','Location','southeast')
% 
% 
% %plot thorax temp as a function of air temp
% figure(3) %create a new plot window for the lookup table plot
% hold on
% air=0:50;
% thorax=air;
% plot(air,Temps(:,2),'y-')  %active, shivering
% plot(air,Temps(:,3),'r')  %passive, flying
% plot(air,Temps(:,4),'r-.')  %active, flying
% plot(air,thorax,'k:')
% %title('Thorax temp')
% ylim([10,45])
% xlim([0,55])
% xlabel('Air Temperature (C)')
% ylabel('Equilibrium Thorax Temperature (C)')
% yline(30);  %30 is the min thorax temp for flight, Heinrich1983
% yline(45);   %lethal thorax temp Heinrich1976?ish
% II=area([0,55],[45 45],42,'EdgeColor', 'none', 'FaceColor', 'r');       
% alpha(0.1)
% legend('active-shivering','passive-flying','active-flying','T_{th} = T_{air}','Location','southeast');
% hold off;



