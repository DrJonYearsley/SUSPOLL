function out = dTth_dt(~, Tth, params)

% Function to calculate rate of change of thorax temp by accounting
%  for rates of heat loss and gain
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define parameters required for each term
param_I = params;
% Update reference metabolic rate for resting
param_I.i0 = params.i0_rest;


param_H = params;
% Update surface area for the head
param_H.A_th = params.A_h;

% Thorax metabolic heat gain (anonymous function)
Q_I = @(Tth, P) {P.i0 (P.M_b/P.M_ref)^3/4 * exp(-P.E/P.k *(1/Tth - 1/P.T_ref))};

% Heat lost to abdomen (behaviour for bumblebees)
Q_ab =  @(Tth, P) {P.active * P.r0 * (Tth - P.T_air)};

% Heat gain from longwave radiation (i.e. atmosphere)
Q_L = @(Tth, P) {P.alpha_np * P.A_th * P.epsilon_a * ...
    (P.sigma * P.Tg^4 + P.delta*P.Tair^6)};

% Radiative heat loss
Q_R = @(Tth, P) {P.epsilon_e * P.A_th * P.sigma * Tth^4};


% Head rate of heat gain (assuming head is in thermal equilibrium with the thorax)
Q_H = @(T_head, P) {Q_S(T_head, P) + Q_L(T_head, P) - ...
    Q_R(T_head, P) -  Q_C(T_head, P) - Q_N(T_head, P)};



% Rate of change of thorax temperature
out = 1/ params.M_th / params.c * (Q_S(Tth, params) + Q_L(Tth, params) + ...
    Q_I(Tth, params) + Q_H(Tth-bee1.delta_Th, param_H) - ...
    Q_R(Tth, params) - Q_C(Tth, params) - Q_ab(Tth, params));












%% Conductive heat loss
function qc = Q_C(Tth, P)

qc = 



% Heat gain from shortwave radiation (i.e. solar radiation, direct and reflected)
Q_S = bee2.A_th * bee2.epsilon_a * env.P * ( + bee1.alpha_so * bee1.alpha_si)


% Heat lost from evaporation of nectar  (behaviour for honeybees)
Q_N() =
