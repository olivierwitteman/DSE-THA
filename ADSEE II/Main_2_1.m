clc;
%% ADSEE II
% Inputs

LAMBDA = 0.;    %Wingsweep at 0.25MAC
Wfiml = 2;      %Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48
Vinf = 102;     %Cruise Speed
WSsc = 70;      %Wing loading at the start of the cruise
WSec = 50;      %Wing loading at the end of the cruise

% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound
v = 92.6; % [m/s], 180kts



% cl = 0.30647;
clmax = 2.1;

cambered = 0; % 1 for True, 0 for False
A = 8;
e = 0.75; % default assumed Oswald factor
S_ref = 15;
c = sqrt(S_ref/A);
b = S_ref/c;
S_h = 0.15*S_ref;
S_v = 0.1*S_ref;
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
A_cs = 3;
D = sqrt(A_cs/pi); % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)

%% ADSEE II - Lecture 1

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, Vinf, WSsc, WSec)
%% ADSEE II - Lecture 2

[CLmax, alpha_stall] = clalpha(A, clmax, LAMBDA, CLdes, S_ref)

%% ADSEE II - Lecture 3 - Drag coefficient estimation
%% Component contributions propeller AC (fast method)
Component = ['Wing', 'Fuselage multi-engine', 'Fuselage single-engine', 'Nacelles', 'Tail (hor + ver)', 'misc'];
C_D_Cs = [0.007, 0.08, 0.11, 0.06, 0.008, 0.15];
% Change these according to component name (defined above)
A_Cs = [15, 0, 3.0, 0.1, 0.3];
% Fast_Cd0 = ADSEE_II_Drag.fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref)

%% Component method
% [Fuselage, Wing, horizontal tail, vertical tail]
C_f_c_fuselage = ADSEE_II_Drag.fp_skin_friction(0.1, k, rho, v, L2, mu, a);
C_f_c_wingtail = ADSEE_II_Drag.fp_skin_friction(0.4, k, rho, v, L3, mu, a);
C_f_cs = [C_f_c_fuselage, C_f_c_wingtail, C_f_c_wingtail, C_f_c_wingtail];
option = [2, 1, 1, 1];
FF_cs = [ADSEE_II_Drag.form_factor(option(1), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(2), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(3), L2, D, tc_avg, xc_max, LAMBDA, v, a), ADSEE_II_Drag.form_factor(option(4), L2, D, tc_avg, xc_max, LAMBDA, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = ADSEE_II_Drag.S_wet_c(S_ref, S_h, S_v, D, L1, L2, L3);


cd0_c = ADSEE_II_Drag.tot_comp_drag0(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0);
misc = ADSEE_II_Drag.cD_misc0(0.034, A_cs, L2*D*0.1, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c, c);

total_cD0 = cd0_c + misc

cD = total_cD0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2

L_D = CLdes/cD

%% ADSEE II - Lecture 4
P_req = degtorad(60)/1.3;%requirement of roll rate

%Input here your wing  parameters
c_r = 1.67; % root chord
c_t = 0.67; %tip chord
lambda = 0.; % sweep at leading edge in degrees (positive number)
theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
c_l_alpha = 0.32; % Airfoil lift curve slope
S_ref = 12.3; % Wing surface in square meters
c_d0 = 0.02; % 2D zero lift drag coefficient
V = 190.; %speed in m/s
b = 10.51; %wingspan in meters
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
P = AileronSizing.Intergral(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, c_d0, c_r, da_max, V);
[b1, Inner_Ail_Chord, Outer_Ail_Chord] = AileronSizing.Iteration(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, c_d0, c_r, da_max, V, P, P_req, chordratio_ail_total)


