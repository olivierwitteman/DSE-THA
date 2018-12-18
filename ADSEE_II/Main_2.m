clc;
clear variables;
%% ADSEE II
% Inputs
vars = load('../ADSEE_I/variables_ADSEE_I.mat');

P = double(vars.P);
A = double(vars.A);               % <---- CHANGE FOR ELECTRIC/HYBRID
MAC = double(vars.MAC);           % <---- CHANGE FOR ELECTRIC/HYBRID
MTOW = double(vars.MTOW);         % <---- CHANGE FOR ELECTRIC/HYBRID
OEW = double(vars.OEW);           % <---- CHANGE FOR ELECTRIC/HYBRID
S_ref = double(vars.S);           % <---- CHANGE FOR ELECTRIC/HYBRID
v = double(vars.V_cruise);        % <---- CHANGE FOR ELECTRIC/HYBRID
W4W5 = double(vars.W4W5);         % <---- CHANGE FOR ELECTRIC/HYBRID
W_f = double(vars.W_fuel_total);  % <---- CHANGE FOR ELECTRIC/HYBRID
taper_ratio = double(vars.tr);    % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_c4 = double(vars.sweep_4c); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_c2 = double(vars.sweep_2c); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_LE = double(vars.sweep_LE); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_TE = double(vars.sweep_TE); % <---- CHANGE FOR ELECTRIC/HYBRID
b = double(vars.b);               % <---- CHANGE FOR ELECTRIC/HYBRID
b_h = double(vars.b_h);           % <---- CHANGE FOR ELECTRIC/HYBRID
b_v = double(vars.b_v);           % <---- CHANGE FOR ELECTRIC/HYBRID
V_stall = double(vars.V_stall);
LAMBDA = sweep_c4;    % Wingsweep at 0.25MAC
Wfiml = 0.97 * MTOW * 9.81;      % Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48 !!!!! ELECTRIC = 1
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
N_gear = 3; % Find Raymer!!!!
lambda = taper_ratio; % taper ratio

LAMBDA_ht = sweep_c4; % Sweep at 25% MAC
A_ht =double(vars.A_h); % Aspect ratio horizontal tailwing                   ??????
H_t_over_H_v = 0; % = 0 for conventional tail, 1 for T-tail
LAMBDA_vt = LAMBDA_ht; % Sweep at 25% of vertical tail MAC                    ??????
A_vt = double(vars.A_v); % Aspect ratio vertical tail                         ??????


lambda_vt = 1.; % taper raio vertical tail                                    ??????
lambda_h = 1; %Taper ratio horizontal tail                                    ??????
L_t = 3.9*3.2808; % Tail length, wing quarter MAC to tail quarter MAC in ft   ??????
N_en = double(vars.N); %nr of engines 
W_press = 0 ;%11.9+(V_pr*P_delta)^0.271; %Weight penalty due to pressurization; PROBABLY ZERO FOR OUR DESIGNS BECAUSE WE DON'T PRESSURIZE OUR CABIN
W_l = (MTOW - W_f) * 2.2; %Landing design gross weight

V_t = W_f / (0.840 * 3.79); %Total fuel volume in gallons

L_m = 15.; %Extended length of main landing gear                ??????
L_n = 15.; %Extended nose gear length (inch)                    ??????


%% Engine Weight Inputs
W_en = P/1340 * 2.2; %Total engine weight (each) in pounds              <---- INPUT

%%
V_i = V_t * 1.05; %Integral tanks volume in gallons

N_t = 1; %Number of fuel tanks                                  ??????
W_uav = 0.03 * MTOW * 2.2; %Uninstalled avionics weight in pounds
N_p = 5; %Number of personal onboard

% cl = 0.30647;
cl_cruise  = double(vars.("cl_cruise"));
clmax = double(vars.("CL_max"));               % XXXX

cambered = 0; % 1 for True, 0 for False
e = double(vars.("e"));
c = sqrt(S_ref/A);  %                           ?????? WHICH CHORD IS THIS????

S_ht = double(vars.("S_h"));
S_vt = double(vars.("S_v"));
MAC_ht = b_h/A_ht;
MAC_vt = b_v/A_vt;

tc_avg = double(vars.("tc")); % (t/c)_avg is the average thickness to chord
xc_max = 0.25; % (x/c)_max is the position of maximum thickness         ????????

% C_f_e = 0.0055; % light AC - single engine
% C_f_e = 0.0045; % light AC - twin engine

% k = 0.152E-5; % polished sheet metal
k = 0.634E-5; % smooth paint

% k = 0.052E-5; % smooth molded composite
% L1 = 1.4; % nosecone length                               ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
% L2 = 3.57; % main fuselage length                          ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
% L3 = 8-L1-L2; % tailcone length                               ??????? SHOULD BE DONE WITH DRAWINGS I GUESS?????
% L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate
% A_cs = 2.9;
% D = sqrt(A_cs/pi) % derived from frontal area (even though fuselage may not be cilindrical)


%L1_pos = 'L1:  ';
%L1 = double(input(L1_pos));

%L2_pos = 'L2:  ';
%L2 = double(input(L2_pos));

% L1_pos = 'L1:  ';   % nosecone length
% L1 = double(input(L1_pos));
% 
% L2_pos = 'L2:  ';   % main fuselage length
% L2 = double(input(L2_pos));
% 
% L3_pos = 'L3:  ';   % tailcone length
% L3 = double(input(L3_pos));

L1 = 1.7;
L2 = 2.8;
L3 = 3;

L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate


A_cs = 3;
D = sqrt(A_cs/pi); % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)

P_req = deg2rad(60)/1.3; % requirement of roll rate

% % % % theta = sweep_TE*180/pi; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"

%%%Aileron geometry input (DO NOT CHANGE)!%%%
aileron_length = [0:0.05:b/2]; % aileron length in meters
tau = 0.6 ; % Function of ratio of the aileron chord over the wing chord (aileron effectiveness) (See slide 10 of ADSEE-II lecture 4 of 2016 for the graph, or look in aircraft design by Mohammed Sadraey)
            % The aileron should be placed after the rear spar, this
            % determines the maximum chord ratio
chordratio_ail_total = 0.41;
%chordratio_ail_total = [0.075, 0.19, 0.41, 0.7];
%tau = [0.2, 0.4, 0.6, 0.8];
%da_max = 30. ; %maximum aileron deflection angle in degrees (reference Mohammed Sadraey)

b2 = b/2;
b1 = b2/2;


%% ADSEE II - Lecture 1

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, v, WSsc, WSec);
%% ADSEE II - Lecture 2

[CLmax, alpha_stall,CL_alpha_clean] = clalpha(A, clmax, CLdes, sweep_c2, M, sweep_LE, tc_avg);

%% ADSEE II - Lecture 3 - Drag coefficient estimation
%% Component contributions propeller AC (fast method)
Component = ['Wing', 'Fuselage multi-engine', 'Fuselage single-engine', 'Nacelles', 'Tail (hor + ver)', 'misc'];
C_D_Cs = [0.007, 0.08, 0.11, 0.06, 0.008, 0.15];
% Change these according to component name (defined above)
A_Cs = [15, 0, 3.0, 0, 0.3];
Fast_Cd0 = ADSEE_II_Drag.fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref);

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


total_cD0 = cd0_c + misc;

%cD = total_cD0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2
cD = Fast_Cd0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2


L_D = CLdes/cD


%% ADSEE II - Lecture 4
c_r = double(vars.("cr"));
c_t = double(vars.("ct"));
sweep_LE; % sweep at leading edge in degrees (positive number)
%theta = atan((c_r-c_t)/(b/2.))*180/pi;
% theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"

c_l_alpha = double(CL_alpha_clean)*180/pi %[rad]
% c_l_alpha = 0.32; % Airfoil lift curve slope  <------- INPUT FROM BOOK
%S_ref = S_ref; % Wing surface in square meters
c_d0 = Fast_Cd0; % 2D zero lift drag coefficient        
%V = 1.2*V_stall; %speed in m/s                          % XXXXXXX
%b = b; %wingspan in meters

aileron_l = AileronNEW(c_r, c_t, sweep_LE, sweep_TE, c_l_alpha,...
   S_ref, c_d0, V_stall, b);

% [b1, Inner_Ail_Chord, Outer_Ail_Chord] = AileronSizing.Iteration(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v, P, P_req, chordratio_ail_total);
% disp('The total aileron size is from the tip of the wing up until: in [m] from the base of the fuselage'), disp(b1);
% disp('Inner Aileron Chord:'), disp(Inner_Ail_Chord), disp('Inner Aileron Chord:'), disp(Outer_Ail_Chord);

% disp('The total aileron size is from the tip of the wing up until: in [m] from the base of the fuselage'), disp(b1);
% disp('Inner Aileron Chord:'), disp(Inner_Ail_Chord), disp('Inner Aileron Chord:'), disp(Outer_Ail_Chord);

%% ADSEE II - Lecture 6 - Class II Weights

W_battery_hybrid = 321;  %just some values for code testing
W_battery_electric = 2053; %just some values for code testing

W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_ref*10.7639,A,tc_avg,lambda,LAMBDA,W_f*2.2,L/D,W_f*2.2,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M);
W_breakdownHYB = W_breakdown + [0 0 0 0 0 0 0 W_battery_hybrid 0 0 0 0 0 0];
W_breakdownELEC = W_breakdown + [0 0 0 0 0 0 0 W_battery_electric-W_breakdown(8) 0 0 0 0 0 0];

%The variables in the matrix W_breakdown are given below
%W_breakdown = [W_wing, W_horizontaltail, W_verticaltail, W_fuselage, W_mainlandinggear, W_noselandinggear, W_installedengines, W_fuelsystem, W_flightcontrols, W_hydraulics, W_avionics, W_electrical, W_airco_and_anti_ice, W_furnishings]/2.2;

W_total = sum(W_breakdown);
W_totalHYB = sum(W_breakdownHYB);
W_totalELEC = sum(W_breakdownELEC);

config = menu('Do you want the prelimnary XLEMAC results for Design 1 (Hybrid), Design 2 (Fuel) or Design 3 (Electric)?', '1','2','3'); 
weight_pie= W_breakdownELEC;

perc_weights =  round(double((W_breakdownELEC/W_totalELEC)*100));
perc_weightsFUEL=round(double((W_breakdown/W_total)*100));
perc_weightsHYB=round(double((W_breakdownHYB/W_totalHYB)*100));

label_pie = {'Wing '+ string(perc_weights(1))+'%',...
    'Horizontal Tail '+ string(perc_weights(2))+'%',...
    'Vertical Tail '+ string(perc_weights(3))+'%',...
    'Fuselage '+ string(perc_weights(4))+'%',...
    'Main Gear '+ string(perc_weights(5))+'%',...
    'Nose Gear '+ string(perc_weights(6))+'%',...
    'Engines '+ string(perc_weights(7))+'%',... 
    'Battery Pack '+ string(perc_weights(8))+'%',...
    'Flight Controls '+ string(perc_weights(9))+'%',...
    'Hydraulics '+ string(perc_weights(10))+'%',...
    'Avionics '+ string(perc_weights(11))+'%',...
    'Electrical '+ string(perc_weights(12))+'%',...
    'Airco and anti-ice '+ string(perc_weights(13))+'%',...
    'Furnishings '+ string(perc_weights(14))+'%'};

label_pieFUEL = {'Wing '+ string(perc_weightsFUEL(1))+'%',...
    'Horizontal Tail '+ string(perc_weightsFUEL(2))+'%',...
    'Vertical Tail '+ string(perc_weightsFUEL(3))+'%',...
    'Fuselage '+ string(perc_weightsFUEL(4))+'%',...
    'Main Gear '+ string(perc_weightsFUEL(5))+'%',...
    'Nose Gear '+ string(perc_weightsFUEL(6))+'%',...
    'Engines '+ string(perc_weightsFUEL(7))+'%',... 
    'Fuel system '+ string(perc_weightsFUEL(8))+'%',...
    'Flight Controls '+ string(perc_weightsFUEL(9))+'%',...
    'Hydraulics '+ string(perc_weightsFUEL(10))+'%',...
    'Avionics '+ string(perc_weightsFUEL(11))+'%',...
    'Electrical '+ string(perc_weightsFUEL(12))+'%',...
    'Airco and anti-ice '+ string(perc_weightsFUEL(13))+'%',...
    'Furnishings '+ string(perc_weightsFUEL(14))+'%'};

label_pieHYB = {'Wing '+ string(perc_weightsHYB(1))+'%',...
    'Horizontal Tail '+ string(perc_weightsHYB(2))+'%',...
    'Vertical Tail '+ string(perc_weightsHYB(3))+'%',...
    'Fuselage '+ string(perc_weightsHYB(4))+'%',...
    'Main Gear '+ string(perc_weightsHYB(5))+'%',...
    'Nose Gear '+ string(perc_weightsHYB(6))+'%',...
    'Engines '+ string(perc_weightsHYB(7))+'%',... 
    'Battery Pack '+ string(perc_weightsHYB(8))+'%',...
    'Flight Controls '+ string(perc_weightsHYB(9))+'%',...
    'Hydraulics '+ string(perc_weightsHYB(10))+'%',...
    'Avionics '+ string(perc_weightsHYB(11))+'%',...
    'Electrical '+ string(perc_weightsHYB(12))+'%',...
    'Airco and anti-ice '+ string(perc_weightsHYB(13))+'%',...
    'Furnishings '+ string(perc_weightsHYB(14))+'%'};

figure;
pELEC=pie(W_breakdownELEC, label_pie);
set(pELEC(2:2:end),'FontSize',30)
figure;
pFUEL=pie(W_breakdown, label_pieFUEL);
set(pFUEL(2:2:end),'FontSize',30)
figure;
pHYB=pie(W_breakdownHYB, label_pieHYB);
set(pHYB(2:2:end),'FontSize',30)

%Here two different paths are taken to taken xlemac for different wing
%positions

%% CG Calculation (With Questions go to Pieter)
L = L/3.2808; %Changing L to meters for the upcoming calculation
L_t = L_t/3.2808; %Changint L_t to meters

CG_OEW_MAC = 0.13;
syms Xlemac
OEWDES1_ARMS = [Xlemac+0.4*MAC,Xlemac+0.25*MAC+L_t+0.15*MAC_ht,Xlemac+0.25*MAC+L_t+0.15*MAC_vt,0.4*L,...
    Xlemac+0.3*MAC,0.15*L, Xlemac+0.25*MAC,Xlemac+0.4*MAC,0.5*L,0.5*L,0.3*L,0.3*L,Xlemac+0.4*MAC,0.4*L].';
    %14x1 matrix with arms of aircraft design 1 (fuel + battery + wing mounted engines)
    %W_breakdown = [W_wing, W_horizontaltail, W_verticaltail, W_fuselage, W_mainlandinggear, W_noselandinggear, 
    %W_installedengines, W_fuelsystem, W_flightcontrols, W_hydraulics, W_avionics, W_electrical, W_airco_and_anti_ice,
    %W_furnishings]/2.2;

OEWDES2_ARMS = [Xlemac+0.4*MAC,Xlemac+0.25*MAC+L_t+0.15*MAC_ht,Xlemac+0.25*MAC+L_t+0.15*MAC_vt,0.4*L,...
    Xlemac+0.3*MAC,0.15*L,Xlemac+0.25*MAC,Xlemac+0.4*MAC,0.5*L,0.5*L,0.3*L,0.3*L,Xlemac+0.4*MAC,0.4*L].';
    %14x1 matrix with arms of aircraft design 2 (fuel + front fuselage mounted engine)
OEWDES3_ARMS = [Xlemac+0.4*MAC,Xlemac+0.25*MAC+L_t+0.15*MAC_ht,Xlemac+0.25*MAC+L_t+0.15*MAC_vt,0.4*L,...
    Xlemac+0.3*MAC,0.15*L,Xlemac+0.25*MAC,Xlemac+0.4*MAC,0.5*L,0.5*L,0.3*L,0.3*L,Xlemac+0.4*MAC,0.4*L].';
    if config == 1
        eqn1 = Xlemac + CG_OEW_MAC*MAC == W_breakdownHYB*OEWDES1_ARMS/W_totalHYB;
        XLEMAC = double(solve(eqn1, Xlemac))
    elseif config == 2
        eqn2 = Xlemac + CG_OEW_MAC*MAC == W_breakdown*OEWDES2_ARMS/W_total;
        XLEMAC = double(solve(eqn2, Xlemac))
    elseif config == 3
        eqn3 = Xlemac + CG_OEW_MAC*MAC == W_breakdownELEC*OEWDES2_ARMS/W_totalELEC;
        XLEMAC = double(solve(eqn3, Xlemac))
    end
 L = L*3.2808; %Changing L back to feet
 L_t=L_t*3.2808; %Changing L_t back to feet
 
 CG_OEW = double(XLEMAC + CG_OEW_MAC*MAC)
