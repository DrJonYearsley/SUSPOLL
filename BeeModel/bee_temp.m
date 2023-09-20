% Calculate the temperature of a bee's thorax based on rates of heat
% loss and gain
%
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set scenario parameters
species='honeybee';       % Either 'honeybee' or 'bumblebee'
metabolism='resting';     % Either 'resting' or 'active'
cooling = "none";         % Either "none", "head" or "abdomen"
behaviour = "flying";    % Either "flying", "shivering", "onflower" or "resting"
% Note flying bees must be metabolically active
crop = "oilseed";         % Either "oilseed" or "fieldbeans"



% Initial conditions
Tth_initial = 30+273.15;   % Statrting thorax temperature (Kelvin)

% Time span to integrate temperature over (s)
tspan = 0:2000; 

% Call function to define parameters based upon scenario
params = set_parameters(species, metabolism, cooling, behaviour, crop);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++







%% Perform additional calculations to add some parameters
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% Constants for heat convection calculations ++++++++++++++++++++++

% dynamic viscosity of air (needs checking)
mu = (1.458e-6*params.T_air^1.5)/(params.T_air+110.4); 

% in humid air air at sea level 
params.rho = params.Pr/(params.R_specific * params.T_air);        
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water

%kinematic viscosity for humid air %The Shock Absorber Handbook
params.nu = mu / params.rho;

% % Calculation for thermal conductvity of air
params.kappa =  (0.02646*params.T_air^1.5)/(params.T_air+245.4*10^(-12/params.T_air));

% Constants for evaporation calculations +++++++++++++++++++++++++++

% Assume diffusion rate equal to advection rate (Peclet number = 1)
params.D_A = params.v * (2*params.R_0);   

%partial pressure of humid air at temperature T_air (K)
p_air = params.rh * params.A*exp(params.B/params.T_air);  

% density of humid air
rho_humid = ((params.Pr-p_air)*params.MM_air + p_air*params.MM_vapor) / ...
    (params.R_specific * params.T_air);   

% Molar fraction and mass fraction of air
X_air = p_air/params.Pr;  % molar fraction of ambient air
Y_air = 1/( 1 + ((1-X_air)/X_air)*(params.MM_air/params.MM_vapor) );

m_evap_coef = 2*pi* params.R_0 * rho_humid * params.D_A;

params.evap_coef1 = params.h_fg * m_evap_coef * log(1-Y_air);   %constant for a fixed T_air
params.evap_coef2 = -params.h_fg * m_evap_coef;   %*log(1-Y_sfc)

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






%% Calculate rate of change of Tth
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

opts = odeset('Events',@(t,y)myEvent(t,y,maxy),'RelTol',1e-12,'AbsTol',1e-14);

dTth_dt(0,23+273.15, params)

% [t,Tth] = ode45(@(t,Tth) dTth_dt(t, Tth, params), tspan, Tth_initial, opts);
