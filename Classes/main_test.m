%% ADSEE II - Lecture 3 - Drag coefficient estimation
clc

%% Inputs
% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound
v = 92.6; % [m/s], 180kts

cl = 0.30647;
cl = 2.;

cambered = 0; % 1 for True, 0 for False
A = 8;
sweep = 0;
e = 0.75; % default assumed Oswald factor
S_ref = 15;
c = sqrt(S_ref/A);
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

%%

dc = Drag_class();

%% Component method
% [Fuselage, Wing, horizontal tail, vertical tail]
C_f_c_fuselage = dc.fp_skin_friction(0.1, k, rho, v, L2, mu, a);
C_f_c_wingtail = dc.fp_skin_friction(0.4, k, rho, v, L3, mu, a);
C_f_cs = [C_f_c_fuselage, C_f_c_wingtail, C_f_c_wingtail, C_f_c_wingtail];
option = [2, 1, 1, 1];
FF_cs = [dc.form_factor(option(1), L2, D, tc_avg, xc_max, sweep, v, a), dc.form_factor(option(2), L2, D, tc_avg, xc_max, sweep, v, a), dc.form_factor(option(3), L2, D, tc_avg, xc_max, sweep, v, a), dc.form_factor(option(4), L2, D, tc_avg, xc_max, sweep, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = dc.S_wet_c(S_ref, S_h, S_v, D, L1, L2, L3);


cd0_c = dc.tot_comp_drag0(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0);
misc = dc.cD_misc0(0.034, A_cs, L2*D*0.1, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c, c);

total_cD0 = cd0_c + misc

cD = total_cD0 + dc.k_f(A, sweep, cl) * cl^2

