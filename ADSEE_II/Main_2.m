clc;
clear variables;
%% ADSEE II
% Inputs

vars = load('../ADSEE_I/variables_ADSEE_I.mat');

A = double(vars.A);
MTOW = double(vars.MTOW);
OEW = double(vars.OEW)
S_ref = double(vars.S);
v = double(vars.V_cruise);
W4W5 = double(vars.W4W5);
W_f = double(vars.W_fuel_total);
taper_ratio = double(vars.tr);
sweep_c4 = double(vars.sweep_4c);
sweep_c2 = double(vars.sweep_2c);
sweep_LE = double(vars.sweep_LE);
sweep_TE = double(vars.sweep_TE);
b = double(vars.b);
V_stall = double(vars.V_stall);

LAMBDA = sweep_c4;    % Wingsweep at 0.25MAC
Wfiml = 0.97 * MTOW * 9.81;      % Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48
WSsc = MTOW * 0.98 * 9.81 / S_ref;      % Wing loading at the start of the cruise
WSec = W4W5 * WSsc;     % Wing loading at the end of the cruise


%%%%%%%% added martin isa
h = double(vars.h);
lambda_isa = 0.0065;
T_0_isa = 288.15;
g_isa = 9.81;
R_isa = 287.1;
P0_isa=101.325*10.^3;
T_isa = T_0_isa - lambda_isa * h;
P_isa = P0_isa * (T_isa / T_0_isa)^(g_isa / (lambda_isa * R_isa));
rho = P_isa / (R_isa * T_isa);
a = sqrt(1.4*R_isa*T_isa); % speed of sound
M = v/a;


%inputs (page 477-479 Raymer) (everything is in retard units) (lbs,
%gallons, ft^3, ft^2, inch etc.)
W_dg = MTOW * 0.9 * 2.2; % Design gross weight
N_z = 4.4; % Load factor
N_gear = 2.5; % Find Raymer!!!!
lambda = taper_ratio; % taper ratio

LAMBDA_ht = sweep_c4; % Sweep at 25% MAC
A_ht = double(vars.A_h); % Aspect ratio horizontal tailwing                   ??????
H_t_over_H_v = 0.; % = 0 for conventional tail, 1 for T-tail
LAMBDA_vt = LAMBDA_ht; % Sweep at 25% of vertical tail MAC      ??????
A_vt = double(vars.A_v); % Aspect ratio vertical tail                         ??????
lambda_vt = 1.; % taper raio vertical tail                      ??????
lambda_h = 1; %Taper ratio horizontal tail                      ??????
L_t = 4.; % Tail length, wing quarter MAC to tail quarter MAC   ??????
W_press = 0 ;%11.9+(V_pr*P_delta)^0.271; %Weight penalty due to pressurization; PROBABLY ZERO FOR OUR DESIGNS BECAUSE WE DON'T PRESSURIZE OUR CABIN
W_l = (MTOW - W_f) * 2.2; %Landing design gross weight

V_t = W_f / (0.840 * 3.79); %Total fuel volume in gallons

L_m = 15.; %Extended length of main landing gear                ??????
L_n = 15.; %Extended nose gear length (inch)                    ??????


W_en = 345. * 2.2; %Engine weight (each) in pounds              <---- INPUT
N_en = double(vars.("N")); %Number of engines\                  XXXX
V_i = V_t * 1.05; %Integral tanks volume in gallons

N_t = 1; %Number of fuel tanks                                  ??????
W_uav = 0.03 * MTOW * 2.2; %Uninstalled avionics weight in pounds
N_p = 5; %Number of personal onboard

% cl = 0.30647;
cl_cruise  = double(vars.("cl_cruise"));
clmax = double(vars.("CL_max"));                % XXXX

cambered = 0; % 1 for True, 0 for False
e = double(vars.("e"));
c = sqrt(S_ref/A);  %                           ?????? WHICH CHORD IS THIS????

S_ht = double(vars.("S_h"));
S_vt = double(vars.("S_v"));
tc_avg = double(vars.("tc")); % (t/c)_avg is the average thickness to chord
xc_max = 0.25; % (x/c)_max is the position of maximum thickness         ????????

% C_f_e = 0.0055; % light AC - single engine
% C_f_e = 0.0045; % light AC - twin engine

% k = 0.152E-5; % polished sheet metal
k = 0.634E-5; % smooth paint
% k = 0.052E-5; % smooth molded composite

L1 = 1.4; % nosecone length                               ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
L2 = 3.57; % main fuselage length                          ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
L3 = 8-L1-L2; % tailcone length                               ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate
A_cs = 2.9;
D = sqrt(A_cs/pi) % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)

P_req = deg2rad(60)/1.3; % requirement of roll rate


% % % % theta = sweep_TE*180/pi; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
c_l_alpha = 0.32; % Airfoil lift curve slope


%%%Aileron geometry input (DO NOT CHANGE)!%%%
aileron_length = [0:0.05:b/2]; % aileron length in meters
tau = 0.6 ; % Function of ratio of the aileron chord over the wing chord (aileron effectiveness) (See slide 10 of ADSEE-II lecture 4 of 2016 for the graph, or look in aircraft design by Mohammed Sadraey)
            % The aileron should be placed after the rear spar, this
            % determines the maximum chord ratio
chordratio_ail_total = 0.41;
%chordratio_ail_total = [0.075, 0.19, 0.41, 0.7];
%tau = [0.2, 0.4, 0.6, 0.8];
da_max = 30. ; %maximum aileron deflection angle in degrees (reference Mohammed Sadraey)

b2 = b/2;
b1 = b2/2;


%% ADSEE II - Lecture 1

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, v, WSsc, WSec)
%% ADSEE II - Lecture 2

[CLmax, alpha_stall] = clalpha(A, clmax, LAMBDA, CLdes, S_ref)

%% ADSEE II - Lecture 3 - Drag coefficient estimation
%% Component contributions propeller AC (fast method)
Component = ['Wing', 'Fuselage multi-engine', 'Fuselage single-engine', 'Nacelles', 'Tail (hor + ver)', 'misc'];
C_D_Cs = [0.007, 0.08, 0.11, 0.06, 0.008, 0.15];
% Change these according to component name (defined above)
A_Cs = [15, 0, 3.0, 0, 0.3];
Fast_Cd0 = ADSEE_II_Drag.fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref)

S_w = ADSEE_II_Drag.S_wet_c(S_ref, S_ht, S_vt, D, L1, L2, L3);

%% Component method
% [Fuselage, Wing, horizontal tail, vertical tail]
C_f_c_fuselage = ADSEE_II_Drag.fp_skin_friction(0.1, k, rho, v, L2, mu, a);
C_f_c_wingtail = ADSEE_II_Drag.fp_skin_friction(0.4, k, rho, v, L3, mu, a);
C_f_cs = [C_f_c_fuselage, C_f_c_wingtail, C_f_c_wingtail, C_f_c_wingtail];
option = [2, 1, 1, 1];
FF_cs = [ADSEE_II_Drag.form_factor(option(1), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(2), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(3), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(4), L2, D, tc_avg, xc_max, LAMBDA, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = ADSEE_II_Drag.S_wet_c(S_ref, S_ht, S_vt, D, L1, L2, L3);


cd0_c = ADSEE_II_Drag.tot_comp_drag0(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0);
misc = ADSEE_II_Drag.cD_misc0(0.034, A_cs, L2*D*0.1, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c, c);

total_cD0 = cd0_c + misc

% cD = total_cD0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2
cD = Fast_Cd0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2

L_D = CLdes/cD


%% ADSEE II - Lecture 4
c_r = double(vars.("cr"));
c_t = double(vars.("ct"));
sweep_LE; % sweep at leading edge in degrees (positive number)
theta = atan((c_r-c_t)/(b/2.))*180/pi;
% theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"

prompt_dclda = 'What is your lift curve slope: default is 0.32  ';
c_l_alpha = double(input(prompt_dclda));
c_l_alpha = 0.32; % Airfoil lift curve slope  <------- INPUT FROM BOOK
S_ref = S_ref; % Wing surface in square meters
c_d0 = Fast_Cd0; % 2D zero lift drag coefficient        
V = 1.2*V_stall; %speed in m/s                          % XXXXXXX
b = b; %wingspan in meters

aileron_l = aielron_22222(c_r, c_t, lambda, theta, c_l_alpha,...
    S_ref, c_d0, V, b);


% P = AileronSizing.Intergral(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v);
% [b1, Inner_Ail_Chord, Outer_Ail_Chord] = AileronSizing.Iteration(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v, P, P_req, chordratio_ail_total);
% P = AileronSizing.Intergral(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v);
% [b1, Inner_Ail_Chord, Outer_Ail_Chord] = AileronSizing.Iteration(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v, P, P_req, chordratio_ail_total);
% disp('The total aileron size is from the tip of the wing up until: in [m] from the base of the fuselage'), disp(b1);
% disp('Inner Aileron Chord:'), disp(Inner_Ail_Chord), disp('Inner Aileron Chord:'), disp(Outer_Ail_Chord);


aileron_l = aielron_22222(c_r, c_t, sweep_LE*180/pi, theta, c_l_alpha,...
    S_ref, c_d0, V, b);

% disp('The total aileron size is from the tip of the wing up until: in [m] from the base of the fuselage'), disp(b1);
% disp('Inner Aileron Chord:'), disp(Inner_Ail_Chord), disp('Inner Aileron Chord:'), disp(Outer_Ail_Chord);
%%

disp(["Final answer: ", num2str(aileron_l)])

%% ADSEE II - Lecture 6 - Class II Weights


%W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_w,A,tc_avg,lambda,LAMBDA,S_f,L_over_D,W_fw,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M)

% W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_ref,A,tc_avg,lambda,LAMBDA,S_W,L_D,W_f,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M)

%W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_w,A,tc_avg,lambda,LAMBDA,S_f,L_over_D,W_fw,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M)

W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_ref*10.7639,A,tc_avg,lambda,LAMBDA,W_f*2.2,L/D,W_f*2.2,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M)

%W_breakdown = [W_wing, W_horizontaltail, W_verticaltail, W_fuselage, W_mainlandinggear, W_noselandinggear, W_installedengines, W_fuelsystem, W_flightcontrols, W_hydraulics, W_avionics, W_electrical, W_airco_and_anti_ice, W_furnishings]/2.2;

FG_OEW_arms = [Xlemac+0.4*MAC Xlemac+0.25*MAC+L_t+0.15*MAC_lt Xlemac+0.25*MAC+L_t+0.15*MAC_ht 0.4*L Xlemac+0.3*MAC 0.15*L 0.15*L_OR_Xlemac+0.25*MAC 0 0.5*L 0.5*L 0.3*L 0.3*L Xlemac+0.4*MAC 0.4*L] %Still put in if/or statement and define variables.
W_total = sum(W_breakdown)
CG_OEW = Xlemac +0.13*MAC
CG_OEW = (W_breakdown*FG_OEW_arms)/(W_total)


