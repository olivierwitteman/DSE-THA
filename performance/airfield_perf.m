%%% AIRFIELD TO PERFORMANCE

%% Assumptions

%% TAKE OFF
%% INPUTS
W_to=1667*9.81          ;                                       %%INPUT 
S=12.1                  ;                                       %%INPUT 
rho=1.225               ;                                       %%INPUT    
Cl_max=2.               ;                                       %%INPUT 
Cd_to=0.505             ;                                       %%INPUT
g=9.81                  ;                                       %%INPUT    
T_av_to=1               ;     341kW                             %%INPUT                                           %%INPUT 
Dg_av_to=1              ;                                       %%INPUT 
gamma_to=1              ;                                       %%INPUT                  
H_scr=15.24             ;%[M]  %input from guide                %%INPUT

%% Calculations
%Ground distance
V_lof = 1.05*V_min

V_min = sqrt(2*W_to/(S*rho*Cl_max))

D_av_to = 0.5 * rho * S * Cd_to * (V_lof/sqrt(2))^2              

a_av = (g/W_to)*(T_av_to-D_av_to-Dg_av_to)

X_gd_to = V_lof^2/(2*a_av)

% Transition phase
X_tp_to = V_lof^2*sin(gamma_to)/(0.15*g)

% Climb distance
X_cd_to = H_scr-(1-cos((gamma_to)*(V_lof^2/(0.15*g))))/tan(gamma_to)

% TO total distance
X_TO = X_gd_to + X_tp_to + X_cd_to



%% LANDING
%% INPUTS

W_la = 1                ;                                       %%INPUT
T_av_la = 1             ;                                       %%INPUT
D_av_la = 1             ;                                       %%INPUT
Dg_av_la = 1            ;                                       %%INPUT
gamma_la = 1            ;                                       %%INPUT
Trans_time = 2          ; %[sec] time to transition             %%INPUT
mu = 1                  ;                                       %%INPUT
delta_n = 0.15*9.81     ; %[] 0.15*g for prop                   %%INPUT
T_av_rev = 1            ;                                       %%INPUT
L_av_la = 1             ;                                       %%INPUT

%% Calculations

V_ap = 1.3*V_min*Cl_max*(delta_n*g)

% Airborne distance

R = 1.3^2*W_la*2/(S*rho*Cl_max*delta_n*g)

X_ab_la = R * sin(gamma_la) + (H_scr-(1-cos(gamma_la))*R/tan(gamma_la))

% Transition phase
X_tp_la = Trans_time * V_ap

% Braking phase

X_br_la = W_la^2*2*1.3^2/(2*g*S*rho*Cl_max*(T_av_rev+D_av_la + mu*(W_la-L_av_la)))
 
X_LA = X_br_la + X_tp_la + X_ab_la






