clc;
clear variables;
% <------ INPUT    Means you can change/adjust it

%% Read reference aircraft for getting starting MTOW form average weights (Test4Terence)
filename = 'Reference_Olivier.xlsx';    % reference file with aircraft
sheet = 1;
MTOW_weights_set = 'C17:C40'; % where does matlab have to look in the excel
MTOW_weights = xlsread(filename,sheet,MTOW_weights_set);
average_MTOW = mean(MTOW_weights);
MTOW_median = median(MTOW_weights)



%% General input and input for Fuel fractions
% Input MTOW, A, e, V_cruise, V_stall
%
% Output: Cd0, LD_cruise (also loiter if needed),
% Wl_Wto (for wingloading),
% W4W5 (for wing loading),
% OEW, W_fuel_used (cruise or all????)

%MTOW = average_MTOW; % or choose average MTOW % <------ INPUT
MTOW=5298*9.80665




A = 8.5;                              % <------ INPUT
e = 0.78;                             % <------ INPUT

V_cruise = 180;  % kts                % <------ INPUT
V_stall = 61;    % kts                % <------ INPUT


[summary, Wl_Wto, Cd0, LD_cruise, W4W5, m_cruise, V_stall] = Fuel_Frac(MTOW, A, e, V_cruise, V_stall); % Reference aircraft and fuel fractions
%                                       !!!!!!!!

%% Input for wingloading and payload range diagram
h = 2400;                   % <------ INPUT
g = 9.80665;
rho = 1.225;

CL_to = 1.9; % <------ INPUT
CL_max = 2.1; % <------ INPUT
c = 5;       % 1.2*V_stall*grad(0.083) = 3.1 minimum . % <------ INPUT
V_land = 1.2*32; % ms From requirements?     % <------ INPUT

OEW = 4935*9.80665;                   % Input from fuel fractions
W_fuel_used = summary(3,2);           % Input from fuel fractions
Wl_Wto;                               % Input from fuel fractions
cd0_clean = Cd0;                      % Input from fuel fractions
A = double(summary(4,2));             % Input from fuel fractions
e_clean = double(summary(5,2));       % Input from fuel fractions
V_stall = double(summary(6, 2));      % Input from fuel fractions
V_cruise = double(summary(7, 2));     % Input from fuel fractions
m_cruise = double(summary(8, 2));     % Input from fuel fractions  !!!!!
W_fuel_total = double(summary(9, 2)); % Input from fuel fractions  !!!!!

Wing_Loading_Func(h,A,e_clean,cd0_clean, CL_max,CL_to ,c, Wl_Wto, V_land, V_stall, V_cruise)

W4W5;                               % Input from fuel fractions
payload = 363;                      % from requirements
percent_emptiness_payload = 00.5;    % Input : how empty the payload is


PR_func(MTOW, OEW, W_fuel_used, payload, LD_cruise, W4W5, percent_emptiness_payload);
% maybe change W_fuel_used to W_fuel_total



% choose design point so we can add S and P
prompt_WS = 'What wing loading did you choose: ';
WS = double(input(prompt_WS))
S = MTOW*9.81/WS; % mˆ2


prompt_WP = 'What WP did you choose: ';
WP = double(input(prompt_WP))
P = MTOW*9.81/WP   % Watts

%% Progress summary
summary_end = ["MTOW: ", MTOW;
     "OEW: ", OEW;
     "W_fuel: ", W_fuel_used;
     "Payload: ", payload;
     "A: ", A;
     "e: ", e;
     "LD_cruise: ", LD_cruise;
     "Cl_max: ", CL_max;
     "CL_to: ", CL_to;
     "V_cruise: ", V_cruise;
     "V_stall: ", V_stall;
     "V_land: ", V_land;
     "Percent emptiness of payload: ", percent_emptiness_payload;
     "h: ", h;
     "Power: ", P;
    ];


%% wing planform BASED ON WING AREA FROM WINGLOADING DIAGRAM
[summary_wing] = wing_planform_design(V_cruise, A, S, m_cruise, h); % m_cruise

summary_wing = [summary_wing; ["Wing Area", S]];



%% eng dimensions
% N = 1
% 
% [D_p, w_ee, l_ee, h_ee] = engine_dim_func(P, N);

%% CG VERY ROUGH ESTIMATION
% % % % fus_length = 6.6;   % <------ INPUT
% % % % X_oew = 0.40;        % <------ INPUT Assume pos5ition of the OEW cg
% % % % X_payload = 0.35;   % <------ INPUT Assume position of the Payload(including passengers) cg
% % % % xc_oewcg = 0.3;     % <------ INPUT
% % % % xc_wcg = 0.4;       % <------ INPUT5
% % % % 
% % % % wing_x = 0.55;       % <------ INPUT Assume position of the wing cg from the nose
% % % % empen_x = 0.8;      % <------ INPUT Assume position of the empennage cg
% % % % fus_x = 0.5;        % <------ INPUT Assume position of the fuselage cg
% % % % nacell_x = 0.01;     % <------ INPUT Assume position of the nacelle cg = same for engines
% % % %                     % fixed equipment is the same position as the fuselage
% % % %                     % cg


tr = double(summary_wing(3, 2));
MAC = double(summary_wing(11, 2));
b = double(summary_wing(2,2));
tc = double(summary_wing(14,2));
cl_cruise = double(summary_wing(8,2));
cr = double(summary_wing(9,2));
ct = double(summary_wing(10,2));
M_cruise = double(summary_wing(1, 2));

sweep_LE = double(summary_wing(6, 2));
sweep_TE = double(summary_wing(7, 2));
sweep_2c = double(summary_wing(5, 2));
sweep_4c = double(summary_wing(4, 2));





%% CG VERY ROUGH ESTIMATION
prompt_prop_pos = 'Where to put the enignes? wing or fuselage . 1/2: ';
prop_pos = double(input(prompt_prop_pos));

if prop_pos == 2 % fuselage
    fus_length = 8.;   % <------ INPUT
    X_oew = 0.40;       % <------ INPUT Assume position of the OEW cg
    X_payload = 0.667;    % <------ INPUT Assume position of the Payload(including passengers) cg
    xc_oewcg = 0.3;     % <------ INPUT
    xc_wcg = 0.4;       % <------ INPUT

    wing_x = 0.8;       % <------ INPUT Assume position of the wing cg from the nose
    empen_x = 0.9;      % <------ INPUT Assume position of the empennage cg
    fus_x = 0.6;        % <------ INPUT Assume position of the fuselage cg
    nacell_x = 0.14*0;    % <------ INPUT Assume position of the nacelle cg = same for engines
                        % fixed equipment is the same position as the fuselage
    propul_x = 0.14     % ADDED WHEN CHANGING STUFF
    N = 1
end


if prop_pos == 1            % wing
    fus_length = 6.6;   % <------ INPUT
    X_oew = 0.550;      % <------ INPUT Assume position of the OEW cg
    X_payload = 0.5;    % <------ INPUT Assume position of the Payload(including passengers) cg
    xc_oewcg = 0.25;    % <------ INPUT
    xc_wcg = 0.4;       % <------ INPUT

    wing_x = 0.8;      % <------ INPUT Assume position of the wing cg from the nose
    empen_x = 0.9;     % <------ INPUT Assume position of the empennage cg
    fus_x = 0.75;      % <------ INPUT Assume position of the fuselage cg
%     nacell_x = -0.2;    % <------ INPUT Assume position of the nacelle cg = same for engines
%                         % fixed equipment is the same position as the fuselage
    propul_x = -0.1     % <--- INPUT in front of the wing
    N = 2;
end


% [x_lemac, most_aft_cg, most_forward_cg] = CG_calc_func(MAC, payload, fus_length, W_fuel_total, MTOW, OEW, X_oew, X_payload, xc_oewcg, xc_wcg)
[x_lemac, most_aft_cg, most_forward_cg] = CG_calc_func(MAC, payload, fus_length, W_fuel_total,...
    double(MTOW), double(OEW), X_oew, X_payload, xc_oewcg, xc_wcg, wing_x, empen_x, fus_x, propul_x, prop_pos);
                                                                                          %!!!!!!!!
                                                                                          
most_aft_cg = most_aft_cg;              % <---- change IF it is too bullshit
most_forward_cg = most_forward_cg;      % <---- CHANGE iF it is too bullshit
%% eng dimensions
% N = 1;
[D_p, w_ee, l_ee, h_ee] = engine_dim_func(P, N);
%% Horizontal and vertical control surface areas

SF_S = 1.0;
[S_h, S_v] = control_surf_func(MAC, S, b, fus_length, empen_x, most_aft_cg);
S_h = SF_S*S_h;
S_v = SF_S*S_v;

A_h = 5.6;    % <----- INPUT   [3, 5] slide 68 lecture 7 ADSEE 1
A_v = 1.84;  % <----- INPUT   [1, 2] slide 68 lecture 7 ADSEE 1

b_h = sqrt(A_h * S_h);
b_v = sqrt(A_v * S_v);

CL_to_end = MTOW*g/(0.5*rho*S*(V_stall*1.2)^2);
CL_max_end = CL_to * 1.1^2;              

save('variables_ADSEE_I.mat', 'A', 'MTOW', 'OEW', 'S', 'V_cruise', 'W4W5',...
    'W_fuel_used', 'tr', 'sweep_LE', 'sweep_TE', "sweep_2c",'sweep_4c',...
    "h", "S_h", "S_v", "b", "N", "W_fuel_total", "CL_max", "CL_to", "e",...
    "tc", "cl_cruise", "cr", "ct", "b", "V_stall",...
    "h", "S_h", "S_v", "b", "N", "W_fuel_total", "M_cruise", "MAC",...
    "A_h", "A_v", "b_h", "b_v")

