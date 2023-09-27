function out = set_parameters(species, cooling, behaviour, crop)
% Set the parameters for the bee thermodynamic model
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check the input arguments
if species~= "honeybee" && species~="bumblebee"
    error("Species must be either bumblebee or honeybee")
end

if cooling~= "none" && cooling~="head" && cooling~="abdomen"
    error("Cooling must be either none, head or abdomen")
end


% Define a structure that keeps all the physical constants for the model
phys = struct('k',8.617333262145*10^(-5),...   % Boltzmann's constant
    'delta', 5.31*10^(-13), ...   % Empirical parameter
    'sigma', 5.67*10^(-8), ...    % Stefan-Boltzmann constant
    'e',1.60217663e-19,...        % Charge of an electron (Coulombs)
    'R_specific',287.058,...      % specific gas constant J/kg/K for dry air
    'A', 9.1496e10, ...           % clausius-clapyron constant A for water in N/m^2 (Sidebottom 2015)
    'B', -5.1152e3, ...           % clausius-clapyron constant B for water in K  (Sidebottom 2015)
    'kappa',0.02534, ...          % thermal conductivity of air
    'h_fg', 2.3819e6, ...           % latent heat of vaporization of water, J/kg 
    'MM_air', 0.0289652, ...      % molar mass of dry air in kg/mol
    'MM_vapor', 0.018016);        % molar mass of water vapor in kg/mol


% This will be calculated using Peclet number assumption Pc=1 (meaning rate
% of diffusion=rate of advection)
%    'D_A', 2.06e-5, ...           % diffusion coefficient of air into itself in m^2/s



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for the environment
env = struct('P', 332.3878, ... % mean solar irradiance from all of Arrian's data
    'T_air', 19 + 273.15, ...   % Air temperature (K)
    'T_g', 17.1 + 273.15, ...   % Ground temperature (K)
    'Pr', 1.013e5,...           % atmospheric pressure in N/m^2 (value used in Sidebotham)
    'f', 0.25, ...              % fraction of solar radiation from sun reflected back by earth (albedo/ground reflectance)
    'rh', 0.6908);              % mean relative humidity from Arrian's data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for all bees
bee.s = 0.9965;             % fraction of internal temp at surface
                            % calculated from Church1960 data converted to K
bee.c = 3.349;              % specific heat capacity of thorax
                            % (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
bee.Cd = 2.429809*10^(-7);  % fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
bee.n = 1.975485;           % fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
bee.delta_Th = 2.9;         % temp difference between head and thorax (K)
bee.alpha_p = 0.25;         % shape factor for incoming solar radiation (Cooper1985)
bee.alpha_np = 0.5;         % fraction of bee's surface (shape factor) that is irradiated with
                            % outgoing solar radiation or thermal radiation (Cooper1985)
bee.Tth_max=70+273.15;      % max thorax temp (K). Stop calculations if temp exceeds this
bee.D_fl = 0.95;            % relative decrease in E and i0 from flying to resting (but metabolically active_

% Indicator parameter to turn cooling behaviour on or off
if cooling=="none"
    bee.active_head = 0;
    bee.active_abdomen = 0;
elseif cooling=="head"
    bee.active_head = 1;
    bee.active_abdomen = 0;
elseif cooling=="abdomen"
    bee.active_head = 0;
    bee.active_abdomen = 1;
end



%%%%%% Forage Patch values (needed in Set_Rates)
bee.distance_to_patch = 100;   %distance to patch in m
bee.t_bout = 3600;             % (s) 60 min bout default for now

% t_forage = t_bout - 2*t_0;   %time available for foraging
% flowervisits_max = floor((t_forage-t_handling)/(t_handling + t_flybetween));   %number of flowers possible to visit
% t_bout_actual = 2*t_0+flowervisits_max*(t_handling+t_flybetween)+t_handling;
% %using t_forage-t_handling because flower visiting has to start and end
% %with a flower visit



%E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
%E_opts = [0.63 0 0]; 
%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
%r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour
%R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue




   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for a bumblebee
bumblebee.A_th = 9.3896e-5;             % thorax surface area in m^2, from Church1960
bumblebee.A_h = 3.61375e-5;             % head surface area in m^2, from Cooper1985 - will need to update this to BB


if behaviour=="flying"
    bumblebee.i0 =  0.06404973;        % reference active metabolic rate (fitted) (J/s)
    bumblebee.E = 0.009176471*phys.e;  % activation energy (fitted) (J)
    bumblebee.v = 4.1;                 % flight speed (m / s) (could be wind speed too)
elseif behaviour=="shivering"
    bumblebee.i0 =  0.06404973;        % reference active metabolic rate (fitted) (J/s)
    bumblebee.E = 0.009176471*phys.e;  % activation energy (fitted) (J)
    bumblebee.v = 0;                   % bee is out of wind
elseif behaviour=="onflower"
    bumblebee.i0 =  0.06404973 * bee.D_fl + 0.001349728 *(1-bee.D_fl); % (J/s)
    bumblebee.E = (0.63 * bee.D_fl + 0.009176471 * (1-bee.D_fl))*phys.e;         % activation energy (fitted) (J)
    bumblebee.v = 0;                 % flight speed (m / s) (could be wind speed too)
elseif behaviour=="resting"
    bumblebee.i0 =  0.001349728;       % reference resting metabolic rate (J/s) Kammer & Heinrich 1974 (table 1)
    bumblebee.E = 0.63*phys.e;         % activation energy (fitted) (J)
    bumblebee.v = 0;                   % bee is out of wind
else
    error("behaviour not recognised for a metabolically active bumblebee")
end


bumblebee.M_b = 0.149;              % mass of the bee in g, Joos1991, default
bumblebee.M_th = 0.057;            %mass of thorax in g, Joos1991
bumblebee.M_ref = 0.177;           % Reference mass of bee Kamer (g)
    %bubblebee.M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
bumblebee.d_th = 0.005467;         % thorax diameter (m)
                                   % (avg thorax diam, from Mitchell1976/Church1960)
bumblebee.epsilon_a = 0.935;       % absorptivity of bees (Willmer1981, ,te)
bumblebee.epsilon_e = 0.97;       % emmisivity of bee
bumblebee.T_mK = 42+273.15;       % median temp for  abdomen cooling (K)
bumblebee.T_ref = 25+273.15;  % Reference temp is 25C for Kammer (resting & flying)
bumblebee.r0 = 0.0367/9;      % rate of heat transfer to abdomen J/s/K
                              % Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
bumblebee.R_0 = 0.00;         % radius of nectar droplet (m)
bumblebee.LethalTemp = 45;    %lethal thorax/air temp (C)
bumblebee.CoolingTemp = 42;   %thorax temp where cooling begins (C)
bumblebee.FlyingTemp = 30;    %thorax temp where flight can begin (C)


bumblebee.Tth_initial = 30+273.15;      %initial temperature of the bee's thorax in K 



if crop=="oilseed"
    bumblebee.t_handline = 5.03;   % Flower handling time (s?)
    bumblebee.t_flybetween = 2.23; % Flower handling time (s?)
elseif crop=="fieldbeans"
    bumblebee.t_handline = 3.4;    % Flower handling time (s?)
    bumblebee.t_flybetween = 2.23; % Flower handling time (s?)
else
    error("Crop not recognised")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for a honeybee


honeybee.A_th = 4.5e-5;            % thorax surface area in m^2, from ???
honeybee.A_h = 2.46e-5;            % head surface area in m^2, from Cooper1985 

if metabolism=="active"
    honeybee.i0 = 0.0335;          % reference active metabolic rate (fitted) (J/s)
else
    honeybee.i0 = 0.000425;        % reference resting metabolic rate (J/s)
end
honeybee.E = 0.008*phys.e;         % activation energy (fitted) (J)


if behaviour=="flying"
    honeybee.i0 =   0.034027;        % reference active metabolic rate (fitted) (J/s)
    honeybee.E = 0.008*phys.e;         % activation energy (fitted) (J)
    honeybee.v = 5.6;                 % flight speed (m / s) (could be wind speed too)
elseif behaviour=="shivering"
    honeybee.i0 =   0.034027;        % reference active metabolic rate (fitted) (J/s)
    honeybee.E = 0.008*phys.e;         % activation energy (fitted) (J)
    honeybee.v = 0;                   % bee is out of wind
elseif behaviour=="onflower"
    honeybee.i0 =   0.034027 * bee.D_fl + 5.65e-3*(80/1000)* (1-bee.D_fl); % (J/s)
    honeybee.E = (0.008 * bee.D_fl + 0.63* (1-bee.D_fl))*phys.e;         % activation energy (fitted) (J)
    honeybee.v = 0;                 % flight speed (m / s) (could be wind speed too)
elseif behaviour=="resting"
    % i0 from Rothe1989, mW/g -> W, 80mg reference mass
    honeybee.i0 = 5.65e-3*(80/1000);  % reference resting metabolic rate (J/s)
    honeybee.E = 0.63*phys.e;         % activation energy (fitted) (J)
    honeybee.v = 0;                   % bee is out of wind
else
    error("behaviour not recognised for a metabolically active bumblebee")
end



honeybee.M_b = 0.100;              % mass of the bee in g, Joos1991
honeybee.M_th = 0.0407;            % mass of thorax in g, Joos1991
honeybee.M_ref = 0.08;             % Reference mass of bee Joos1991 (g)

honeybee.d_th = 0.004;             % diameter of thorax (m) (avg thorax diam, from Mitchell1976/Church1960)
honeybee.epsilon_a = 0.91;         % absorptivity of bees (Willmer1981, ,te)
honeybee.epsilon_e = 0.97;         % emmisivity of a bee thorax
honeybee.T_mK = 47.9+273.15;       % median temp for  evaporative cooling
honeybee.T_ref = 25+273.15;        % Reference temp for metabolic rate calculations
honeybee.r0 = 0;                   % rate of heat transfer to abdomen J/s/◦C
honeybee.R_0 = 0.000616/2;         % radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue
honeybee.LethalTemp = 52;          % lethal thorax/air temp in C
honeybee.CoolingTemp = 47.9;       % thorax temp where cooling begins
honeybee.FlyingTemp = 35;          % thorax temp where flight can begin

if crop=="oilseed"
    honeybee.t_handline = 4.0;     % Flower handling time (s?)
    honeybee.t_flybetween = 2.23;  % Flower handling time (s?)
elseif crop=="fieldbeans"
    honeybee.t_handline = 11.9;    % Flower handling time (s?)
    honeybee.t_flybetween = 2.23;  % Flower handling time (s?)
else
    error("Crop not recognised")
end


%% Produce a final structure of all parameters for the bee and the environment
mergestructs = @(x,y) cell2struct([struct2cell(x); struct2cell(y)],[fieldnames(x); fieldnames(y)],1);

if species=="honeybee"
    bee_final = mergestructs(bee, honeybee);
else
    bee_final = mergestructs(bee, bumblebee);
end

% COmbine physical constants and environmental parameters
const_final = mergestructs(phys, env);

% Final parameters to return
out = mergestructs(const_final, bee_final);

   
