clc; clear all
% <------ INPUT    Means you can change/adjust it
%Bonjour
%% Read reference aircraft for getting starting MTOW form average weights
filename = 'Reference_Olivier.xlsx';    % reference file with aircraft
sheet = 1;
MTOW_weights_set = 'C17:C39'; % where does matlab have to look in the excel
MTOW_weights = xlsread(filename,sheet,MTOW_weights_set);
average_MTOW = mean(MTOW_weights);
%hallohallo

% uheafiudhiue
%% General input and input for Fuel fractions
% Input MTOW, A, e, V_cruise, V_stall
%
% Output: Cd0, LD_cruise (also loiter if needed),
% Wl_Wto (for wingloading),
% W4W5 (for wing loading),
% OEW, W_fuel_used (cruise or all????)

MTOW = 1835; % or choose average MTOW % <------ INPUT
A = 7.5;                              % <------ INPUT
e = 0.75;                             % <------ INPUT
V_cruise = 180;  % kts                % <------ INPUT
V_stall = 61;    % kts                % <------ INPUT


[summary, Wl_Wto, Cd0, LD_cruise, W4W5, m_cruise] = Fuel_Frac(MTOW, A, e, V_cruise, V_stall); % Reference aircraft and fuel fractions
%                                       !!!!!!!!

%% Input for wingloading and payload range diagram
h = 2300;                   % <------ INPUT

Cl_max = 2.2;               % <------ INPUT
cl_to = 1.9                 % <------ INPUT
c = 5       % 1.2*V_stall*grad(0.083) = 3.1 minimum . % <------ INPUT
V_land = 32 % ms From requirements?     % <------ INPUT

OEW = summary(1,2);                 % Input from fuel fractions
W_fuel_used = summary(3,2);         % Input from fuel fractions
Wl_Wto;                             % Input from fuel fractions
cd0_clean = Cd0                     % Input from fuel fractions
A = double(summary(4,2));           % Input from fuel fractions
e_clean = double(summary(5,2));     % Input from fuel fractions
V_stall = double(summary(6, 2));    % Input from fuel fractions
V_cruise = double(summary(7, 2));   % Input from fuel fractions
m_cruise = double(summary(8, 2));   % Input from fuel fractions  !!!!!
W_fuel_total = double(summary(9, 2))% Input from fuel fractions  !!!!!


Wing_Loading_Func(h,A,e_clean,cd0_clean, Cl_max,cl_to ,c, Wl_Wto, V_land, V_stall, V_cruise)

W4W5;                               % Input from fuel fractions
payload = 363;                      % from requirements
percent_emptiness_payload = 0.5;    % Input : how empty the payload is

PR_func(MTOW, OEW, W_fuel_used, payload, LD_cruise, W4W5, percent_emptiness_payload);
% maybe change W_fuel_used to W_fuel_total



% choose design point so we can add S and P
prompt_WS = 'What wing loading did you choose: ';
WS = double(input(prompt_WS))
S = MTOW*9.81/WS; % mˆ2


prompt_WP = 'What WP did you choose: ';
WP = double(input(prompt_WP))
P = MTOW*9.81/WP;   % Watts

%% Progress summary
summary_end = ["MTOW: ", MTOW;
     "OEW: ", OEW;
     "W_fuel: ", W_fuel_used;
     "Payload: ", payload;
     "A: ", A;
     "e: ", e;
     "LD_cruise: ", LD_cruise;
     "Cl_max: ", Cl_max;
     "Cl_to: ", cl_to;
     "V_cruise: ", V_cruise;
     "V_stall: ", V_stall;
     "V_land: ", V_land;
     "Percent emptiness of payload: ", percent_emptiness_payload;
     "h: ", h;
     "Power: ", P;
    ];


%% wing planform BASED ON WING AREA FROM WINGLOADING DIAGRAM
[summary_wing] = wing_planform_design(V_cruise, A, S, m_cruise, h);

summary_wing = [summary_wing; ["Wing Area", S]]

MAC = double(summary_wing(9, 2))

%% eng dimensions
N = 1
[D_p, w_ee, l_ee, h_ee] = engine_dim_func(P, N);

%% CG VERY ROUGH ESTIMATION
fus_length = 5.2;   % <------ INPUT
X_oew = 0.5;        % <------ INPUT Assume position of the OEW cg
X_payload = 0.35;   % <------ INPUT Assume position of the Payload(including passengers) cg
xc_oewcg = 0.3;     % <------ INPUT
xc_wcg = 0.4;       % <------ INPUT

wing_x = 0.4;       % <------ INPUT Assume position of the wing cg from the nose
empen_x = 0.8;      % <------ INPUT Assume position of the empennage cg
fus_x = 0.5;        % <------ INPUT Assume position of the fuselage cg
nacell_x = 0.4;     % <------ INPUT Assume position of the nacelle cg = same for engines
                    % fixed equipment is the same position as the fuselage
                    % cg

% [x_lemac, most_aft_cg, most_forward_cg] = CG_calc_func(MAC, payload, fus_length, W_fuel_total, MTOW, OEW, X_oew, X_payload, xc_oewcg, xc_wcg)
[x_lemac, most_aft_cg, most_forward_cg] = CG_calc_func(MAC, payload, fus_length, W_fuel_total,...
    double(MTOW), double(OEW), X_oew, X_payload, xc_oewcg, xc_wcg, wing_x, empen_x, fus_x, nacell_x)


% maybe add the input in the wingloading function
