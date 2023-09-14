% Calculate the temperature of a bee's thorax based on rates of heat
% loss and gain
%
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Define a structure that keeps all the physical constants for the model
phys = struct('k',8.617333262145*10^(-5),...   % Bolzmann's constant
    'delta', 5.31*10^(-13), ...   % Empirical parameter
    'sigma', 5.67*10^(-8), ...    % Stefan-Boltzmann constant
    'A', 9.1496*10^10, ...        % clausius-clapyron constant A for water in N/m^2
    'B', -5.1152*10^3, ...        % clausius-clapyron constant B for water in K
    'MM_air', 0.0289652, ...      % molar mass of dry air in g/mol
    'MM_vapor', 0.018016);        % molar mass of water vapor in kg/mol


%D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s - %calculated version is below
%h_fg = 2.33*10^6;  %latent heat of vaporization of water, J/kg 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for the environment
env = struct('P', 332.3878, ... % mean solar irradiance from all of Arrian's data
    'T_gK',17.1 + 273.15, ...   % ground surface temp in K
    'T_air', 25 + 273.15, ...   % Air temperature (K)
    'Pr', 1.013e5,...           % atmospheric pressure in N/m^2 (value used in Sidebotham)
    'R_specific',287.058,...    % specific gas constant J/kg/K for dry air
    'a', 0.25, ...              % fraction of solar radiation from sun reflected back by earth (albedo)
    'rh', 0.6908);              % mean relative humidity from Arrian's data

% Calculate some additional environmental parameters
env.mu = (1.458e-6*env.T_air^1.5)/(env.T_air+110.4); %dynamic viscosity of air
env.rho = env.Pr/(phys.R_specific*env.T_air);  %in humid air air at sea level 
%kappa = 0.02534;   %15C dry air at sea level    %thermal conductivity of air (fill in an equation )
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for all bees
bee.s = 0.9965;             % fraction of internal temp at surface
                            % calculated from Church1960 data converted to K
bee.c = 3.349;              % specific heat capacity of thorax
                            % (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
bee.C_l = 2.429809*10^(-7); % fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
bee.n = 1.975485;           %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
bee.delta_Th = 2.9;         % temp difference between head and thorax (K)
bee.alpha_si = 0.25;        % shape factor for incoming solar radiation (Cooper1985)
bee.alpha_so = 0.5;         % fraction of bee's surface that is irradiated with
                            % outgoing solar radiation (Cooper1985)
bee.alpha_th = 0.5;         % fraction of surface of bee that is irradiated with
                            % thermal radiation (Cooper1985)
bee.maxy=70+273.15;         % max thorax temp (K). Stop calculations if temp exceeds this
bee.active = 0;             % 1 = active cooling behaviour is happening


tspan = 0:2000; 

%E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
%E_opts = [0.63 0 0]; 
%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
%r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour
%R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue




   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for a bumblebee
bubblebee.A_th = 9.3896*10^(-5)  ; % thorax surface area in m^2, from Church1960
bubblebee.A_h = 3.61375*10^(-5)  ; % head surface area in m^2, from Cooper1985 - will need to update this to BB
bumblebee.i0_active = 0.06335;     % reference active metabolic rate (fitted) (J/s)
bumblebee.i0_rest = 0.00135;       % reference resting metabolic rate (J/s) Kammer & Heinrich 1974
bumblebee.E = 0.02*1.602e-19;      % activation energy (fitted) (J)
bubblebee.M_b = 0.149;             %mass of the bee in g, Joos1991, default
%bubblebee.M_b = 0.035;   %mass of the bee in g, Joos1991, minimum
%bubblebee.M_b = 0.351;   %mass of bee in g, Joos1991, maximum
bubblebee.M_th = 0.057;            %mass of thorax in g, Joos1991
bumblebee.M_ref = 0.177;           % Reference mass of bee Kamer (g)
    %bubblebee.M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
bubblebee.l_th = 0.005467;         % thorax diameter (m)
                                   % (avg thorax diam, from Mitchell1976/Church1960)
bubblebee.epsilon_a = 0.935;       % absorptivity of bees (Willmer1981, ,te)
bubblebee.v_options = [0 0.1 4.1]; % make the bee be out of wind when resting/shivering, default
    %v_options = [0 0 5.5];   %make the bee be out of wind when resting/shivering, maximum (Osborne2013)
    %v_options = [0 0 1];     %make the bee be out of wind when resting/shivering, minimum (Osborne2013)
bubblebee.epsilon_e = 0.97;       % emmisivity of bee
bubblebee.T_mK = 42+273.15;       % median temp for  abdomen cooling (K)

%I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
    %I_flying = 0.06229515;     %Kammer1974, converted to W
    %I_flying = 0.0018375;     %fitted value
    %I_flying = 0.2097035;     %Heinrich1975, converted to W
    %I_flying = 0.03;   %experiment with the value
    %masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
    %masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
bubblebee.masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    %RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
bubblebee.T_ref = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
    %norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3
bumblebee.r0 = 0.004;         % rate of heat transfer to abdomen J/s/◦C
bubblebee.LethalTemp = 45;    %lethal thorax/air temp (C)
bubblebee.CoolingTemp = 42;   %thorax temp where cooling begins (C)
bubblebee.FlyingTemp = 30;    %thorax temp where flight can begin (C)


y0 = 30+273.15;   %initial temperature of the bee's thorax in K 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define structure with parameters for a honeybee


honeybee.A_th = 4.5*10^(-5)  ;     % thorax surface area in m^2, from ???
honeybee.A_h = 2.46*10^(-5)  ;     % head surface area in m^2, from Cooper1985 
honeybee.i0_active = 0.0335;       % reference active metabolic rate (fitted) (J/s)
honeybee.i0_rest = 0.000425;       % reference resting metabolic rate (J/s) 
honeybee.E = 0.008*1.602e-19;      % activation energy (fitted) (J)

honeybee.M_b = 0.100;              % mass of the bee in g, Joos1991
honeybee.M_th = 0.0407;            % mass of thorax in g, Joos1991
honeybee.M_ref = 0.08;             % Reference mass of bee Joos1991 (g)

honeybee.l_th = 0.004;             % diameter of thorax (m) (avg thorax diam, from Mitchell1976/Church1960)
honeybee.epsilon_a = 0.91;         % absorptivity of bees (Willmer1981, ,te)
honeybee.v_options = [0.1 0.1 5.6];% make the bee be out of wind when resting/shivering, default
    % v_options = [0.1 0.1 10];   %make the bee be out of wind when resting/shivering, maximum (fill in ref)
    %v_options = [0.1 0.1 1];   %make the bee be out of wind when resting/shivering, minimum (fill in ref)
honeybee.epsilon_e = 0.97;         % emmisivity of a bee thorax
honeybee.T_mK = 47.9+273.15;     %median temp for  evaporative cooling

    %I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
    %I_flying = 0.4*(80/1000); %Nachtigall1989, W/g -> W, 80mg reference mass%    
    %I_resting = 0.001349728;     %%experiment with value
    %I_flying = 0.004;     %experiment with value
honeybee.masses = [0.08 0.08 0.08];%reference weight for Rothe/Nachtigal (flying)
honeybee.T_ref = [25+273.15, 25+273.15, 25+273.15];   %Reference temp 
honeybee.r0 = 0;                   % rate of heat transfer to abdomen J/s/◦C
    %norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3

honeybee.LethalTemp = 52;          % lethal thorax/air temp in C
honeybee.CoolingTemp = 47.9;       % thorax temp where cooling begins
honeybee.FlyingTemp = 35;           % thorax temp where flight can begin



y0 = 39+273.15;   %initial temperature of the bee's thorax in K (fill in ref)






   
