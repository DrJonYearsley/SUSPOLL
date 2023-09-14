function out = dTth_dt(~, Tth, params)

% Function to calculate rate of change of thorax temp by accounting
%  for rates of heat loss and gain
%
% Jon Yearsley and Sarah MacQueen
% Sept 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define parameters required for each term
param_I = struct('k',params.k, ...
    'E', params.E, ...
    'M_b',params.M_b, ...
    'M_ref',params.M_ref, ...
    'T_ref',params.T_ref, ...
    'i0', params.i0_rest);

param_S = ;
param_L = ;
param_H = struct('k',params.k, ...
    'A_h', bee.A_h, ...
    );
param_R = ;
param_C = ;
param_ab = ;


% Thorax metabolic heat gain (anonymous function)
Q_I = @(Tth, P) {P.i0 (P.M_b/P.M_ref)^3/4 * exp(-P.E/P.k *(1/Tth - 1/P.T_ref))};

% Heat lost to abdomen (behaviour for bumblebees)
Q_ab =  @(Tth, P) {P.active * P.r0 * (Tth - P.T_air)};


% Rate of change of thorax temperature
dTth_dt = 1/ params.M_th / params.c * (Q_S(Tth, params) + Q_L(Tth, params) + ...
    Q_I(Tth, params) + Q_H(Tth-bee1.delta_Th, param_H) - ...
    Q_R(Tth, params) - Q_C(Tth, params) - Q_ab(Tth, params));




%% Head rate of heat gain (assuming head is in thermal equilibrium with the thorax)
function qh = Q_H(T_head)


qh = Q_S(T_head, param_S) + Q_L(T_head, param_L) - ...
    Q_R(T_head, param_R) -  Q_C(T_head, param_C) - ...
    Q_N(T_head, param_N);







%% Conductive heat loss
function qc = Q_C(Tth, P)

qc = 

%% Radiative heat loss
function qr = Q_R(Tth, P)

qr = param_R.epsilon_e * param_R.A_th * param_R.sigma * Tth^4;

%% Heat gain from longwave radiation (i.e. atmosphere)
Q_L = param_L.alpha_np * param_L.A_th * param_L.epsilon_a * ...
    (param_L.sigma * param_L.Tg^4 + param_L.delta*param_L.Tair^6);


% Heat gain from shortwave radiation (i.e. solar radiation, direct and reflected)
Q_S = bee2.A_th * bee2.epsilon_a * env.P * ( + bee1.alpha_so * bee1.alpha_si)


% Heat lost from evaporation of nectar  (behaviour for honeybees)
Q_N() =
