clc;
clear variables;
%% ADSEE II
% Inputs

vars = load('../ADSEE_I/variables.mat');

A = vars.A;
MTOW = vars.MTOW;
OEW = vars.OEW;
S_ref = vars.S;
v = vars.V_cruise;
W4W5 = vars.W4W5;

LAMBDA = 0.;    %Wingsweep at 0.25MAC
Wfiml = 0.97*MTOW*9.81;      %Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48
WSsc = MTOW*0.98*9.81/S_ref;      %Wing loading at the start of the cruise
WSec = W4W5*WSsc;     %Wing loading at the end of the cruise

% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound
% v = 92.6; % [m/s], 180kts
M = v/a;

%inputs (page 477-479 Raymer) (everything is in retard units) (lbs,
%gallons, ft^3, ft^2, inch etc.)
W_dg = 1.; %Design gross weight
N_z = 1.; %Load factor
N_gear = 1;
S_w = 1.; % Wing surface
lambda = 1.; %taper ratio
S_f = 1.; %Wetted area
W_fw = 1.;

LAMBDA_ht = 1.; % Sweep at 25% MAC
A_ht = 1.; % Aspect ratio horizontal tailwing
H_t_over_H_v = 1.; % = 0 for conventional tail, 1 for 1 tail
LAMBDA_vt = 1.; % Sweep at 25% of vertical tail MAC
A_vt = 1.; % Aspect ratio vertical tail
lambda_vt = 1.; % taper raio vertical tail
lambda_h = 1; %Taper ratio horizontal tail
L_t = 1.; % Tail length, wing quarter MAC to tail quarter MAC
W_press = 0 ;%11.9+(V_pr*P_delta)^0.271; %Weight penalty due to pressurization; PROBABLY ZERO FOR OUR DESIGNS BECAUSE WE DON'T PRESSURIZE OUR CABIN
W_l = 1.; %Landing design gross weight
L_m = 1.; %Extended length of main landing gear
L_n = 1. ; %Extended nose gear length (inch)
W_en = 1.; %Engine weight (each) in pounds
N_en = 1; %Number of engines\
V_t = 1. ; %Total fuel volume in gallons
V_i = 1. ;%Integral tanks volume in gallons
N_t = 1; %Number of fuel tanks
W_uav = 1.; %Uninstalled avionics weight in pounds
N_p = 1; %Number of personal onboard

cl = 0.30647;
clmax = 2.1;

cambered = 0; % 1 for True, 0 for False
e = 0.75; % default assumed Oswald factor
c = sqrt(S_ref/A);
b = S_ref/c
S_ht = 0.15*S_ref;
S_vt = 0.1*S_ref;
tc_avg = 0.15; % (t/c)_avg is the average thickness to chord
xc_max = 0.25; % (x/c)_max is the position of maximum thickness

C_f_e = 0.0055; % light AC - single engine
% C_f_e = 0.0045; % light AC - twin engine

S_W = 4 * S_ref; % assumed wetted area

% k = 0.152E-5; % polished sheet metal
k = 0.634E-5; % smooth paint
% k = 0.052E-5; % smooth molded composite

L1 = 1; % nosecone length
L2 = 4; % main fuselage length
L3 = 2; % tailcone length
L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate
A_cs = 3;
D = sqrt(A_cs/pi); % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)


P_req = degtorad(60)/1.3; % requirement of roll rate

%Input here your wing  parameters
c_r = 1.67; % root chord
c_t = 0.67; %tip chord

theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
c_l_alpha = 0.32; % Airfoil lift curve slope

%%%
%b1 = [0:0.5:(b/2-aileron_length)]; %   the length in meters where the aileron starts measured from the wing root
%b2 = b1+aileron_length ; % end aileron '
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
A_Cs = [15, 0, 3.0, 0.1, 0.3];
Fast_Cd0 = ADSEE_II_Drag.fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref)

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
c_r = 1.67; % root chord
c_t = 0.67; %tip chord
lambda = 0.; % sweep at leading edge in degrees (positive number)
theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
c_l_alpha = 0.32; % Airfoil lift curve slope
S_ref = 12.3; % Wing surface in square meters
c_d0 = 0.02; % 2D zero lift drag coefficient
V = 190.; %speed in m/s
b = 10.51; %wingspan in meters


aileron_l = aielron_22222(c_r, c_t, lambda, theta, c_l_alpha,...
    S_ref, c_d0, V, b);


disp("Final answer")
disp(aileron_l)
%% ADSEE II - Lecture 6 - Drag coefficient estimation

%W_breakdown = C2W.calculation(W_dg,N_z,N_gear,S_w,A,tc_avg,lambda,LAMBDA,S_f,L_over_D,W_fw,v,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,b,W_uav,N_p,M)
