function out = dTth_dt(~, Tth, params)
% Function to calculate rate of change of thorax temp by accounting
%  for rates of heat loss and gain
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define parameters required for metabolic term and heat transfer in the head


%% Define functions for each heat transfer term

% Thorax metabolic heat gain (anonymous function)
Q_I = params.i0 * (params.M_b/params.M_ref)^3/4 * ...
    exp(-params.E/params.k *(1./Tth - 1/params.T_ref));

% Heat lost to abdomen (behaviour for bumblebees)
Q_ab =  params.active_abdomen * params.r0 * (Tth - params.T_air);


% Heat lost from evaporation of nectar  (active behaviour for honeybees)
if params.active_head
    X_sfc = (params.rh*params.A*exp(params.B/(Tth-params.delta_Th)))/params.Pr;
    Y_sfc = 1/( 1 + ((1-X_sfc)/X_sfc)*(params.MM_air/params.MM_vapor) );
else
    Y_sfc = 0;
end

Q_N = params.active_head * (params.evap_coef1 + params.evap_coef2 * log(1-Y_sfc));


% Following are defined as anonymous functions because they are used for
% thorax heat transfer and head heat transfer

% Heat gain from longwave radiation (i.e. atmosphere)
Q_L = @(T, P) P.alpha_np * P.A_th * P.epsilon_a * ...
    (P.sigma * P.T_g^4 + P.delta*P.T_air^6);

% Heat gain from shortwave radiation (i.e. direct and reflected sunlight)
Q_S = @(T, P) P.A_th * params.epsilon_a * P.P * (P.alpha_p + P.alpha_np*P.f);

% Radiative heat loss
Q_R = @(T, P) P.epsilon_e * P.A_th * P.sigma * T^4;

% Conductive heat loss
Q_C = @(T, P) (P.kappa*P.Cd/P.d_th) * (P.d_th * P.v / P.nu)^P.n * ...
    P.A_th * (P.s * T - P.T_air);



% Heat transfer parameters for the head
param_H = params;
% Update surface area for the head
param_H.A_th = params.A_h;

% Head rate of heat gain (assuming head is in thermal equilibrium with the thorax)
Q_H = @(T, P) Q_S(T, P) + Q_L(T, P) - Q_R(T, P) -  Q_C(T, P) - Q_N;



%% Rate of change of thorax temperature
out = 1/ params.M_th / params.c * (Q_S(Tth, params) + Q_L(Tth, params) + ...
    Q_I + Q_H(Tth-params.delta_Th, param_H) - Q_R(Tth, params) - ...
    - Q_C(Tth, params) - Q_ab);



