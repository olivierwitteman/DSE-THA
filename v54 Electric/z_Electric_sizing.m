clc

% Really quick and course first iteration on battery weight for a fully
% electric aircraft

E_sp_f = 43; % MJ/kg
W_f = 150; % fuel weight
eta_f_tot = 0.3;


E_req = W_f * eta_f_tot * E_sp_f; % MJ required (without losses)
disp(['Required energy (without losses) [MJ]: ', num2str(E_req)])

E_sp_batt = 0.87*1.5; % MJ/kg
eta_batt_tot = 0.9;

W_batt = E_req / (eta_batt_tot * E_sp_batt);

disp(['Battery mass without iterating for snowball effect [kg]: ', num2str(round(W_batt, 0))])