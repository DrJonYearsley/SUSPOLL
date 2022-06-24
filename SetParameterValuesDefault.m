pcolors=['b' 'y' 'r'];

%%%%%% Constant values %%%%%%%%
Params_Const
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
h_fg = 2.33*10^6;  %latent heat of vaporization of water, J/kg 
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
rh =  0.6908;   %mean relative humidity from Arrian's data
%D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s - %calculated version is below
MM_air = 0.0289652;   %molar mass of dry air in g/mol
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%

P = 332.3878; %mean solar irradiance from all of Arrian's data
T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
Pr = 1.013*10^5; %atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = 287.058;  %J/kg/K for dry air
%kappa = 0.02534;   %15C dry air at sea level    %thermal conductivity of air (fill in an equation )
%https://www.tec-science.com/mechanics/gases-and-liquids/viscosity-of-liquids-and-gases/#Formulas_for_calculating_the_viscosity_of_air_and_water
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%% Bee Parameters same for BB and HB %%%%%%%%%%%%%%%%
%s = 0.9;  %fraction of internal temp at surface - default value calculated from Church1960 data in C
s = 0.9141875;  %fraction of internal temp at surface - BB fitted value
%s = 1; %temporarily make surface and internal temp equal 
%s = 0.9965;  %ratio calculated from Church1960 data converted to K
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
maxy=50+273.15;    %don't solve above 50C because it's not biologically relevant
tspan = 0:2000; 
E_opts = [0.63 0.63 0.63];    %Brown2004 activation energy
%E_opts = [0.63 0 0]; 
%r = 0.0367/20;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_air = 20C
r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
%r = 0.0367;   %from Heinrich1976, in J/s
%r = 0.03;    %playing around with value to get desired behaviour
R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue
   

if Bumblebee==true
   %%%%%%%%%%%%%%% Bee Parameters %%%%%%%%%%%%%%%

    A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
    A_h = 3.61375*10^(-5)  ; %head surface area in m^2, from Cooper1985 - will need to update this to BB
    M_b = 0.149;   %mass of the bee in g, Joos1991, default
    %M_b = 0.035;   %mass of the bee in g, Joos1991, minimum
    %M_b = 0.351;   %mass of bee in g, Joos1991, maximum
    M_th = 0.057;  %mass of thorax in g, Joos1991
    M_h = 0.039;  %mass of head in g, Joos1991 (body-thorax-abdomen; need to account for wings & legs)
    l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0 0.1 4.1];   %make the bee be out of wind when resting/shivering, default
    %v_options = [0 0 5.5];   %make the bee be out of wind when resting/shivering, maximum (Osborne2013)
    %v_options = [0 0 1];   %make the bee be out of wind when resting/shivering, minimum (Osborne2013)
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 42+273.15;     %median temp for  abdomen cooling
    I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
    %I_flying = 0.06229515;     %Kammer1974, converted to W
    %I_flying = 0.0018375;     %fitted value
    %I_flying = 0.2097035;     %Heinrich1975, converted to W
    I_flying = 0.03;   %experiment with the value
    %masses = [M_b M_b M_b];   %reference weight for Kammer only data is just M_b for now
    %masses = [M_b (0.25+0.60)/2 (0.25+0.60)/2];   %reference weight for Heinrich (flying)
    masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    %RefTemps = [25+273.15, 19.55556+273.15, 19.55556+273.15];   %Reference temp is 25C for Kammer (resting), 35-44C for Heinrich (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
    norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3

    y0 = 30+273.15;   %initial temperature of the bee's thorax in K 

end

if Honeybee==true
    A_th = 4.5*10^(-5)  ; %thorax surface area in m^2, from ???
    A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 
    M_b = 0.100;   %mass of the bee in g, Joos1991
    M_th = 0.0407;  %mass of thorax in g, Joos1991
    l_th = 0.004;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.91;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0.1 0.1 3.1];   %make the bee be out of wind when resting/shivering, default
    % v_options = [0.1 0.1 10];   %make the bee be out of wind when resting/shivering, maximum (fill in ref)
    %v_options = [0.1 0.1 1];   %make the bee be out of wind when resting/shivering, minimum (fill in ref)
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 44+273.15;     %median temp for  evaporative cooling

    I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
    %I_flying = 0.4*(80/1000); %Nachtigall1989, W/g -> W, 80mg reference mass%    
    %I_resting = 0.001349728;     %%experiment with value
    I_flying = 0.004;     %experiment with value
    masses = [0.08 0.08 0.08];   %reference weight for Rothe/Nachtigal (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp 
    norm_constants = [I_resting, I_flying, I_flying];   %resting/shivering/flying = 1,2,3

    y0 = 39+273.15;   %initial temperature of the bee's thorax in K (fill in ref)
end