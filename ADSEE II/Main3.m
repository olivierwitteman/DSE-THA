clc; clear all
%%%%%%%%%%%%%
%Input file
%Airfoil selection
%INPUTS
LAMBDA =  0.1;           %Wingsweep at 0.25MAC
rho=      1.225;          %rho at cruise alt
Wfiml=    2;          %Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48
Vinf =    102;          %Cruise Speed
WSsc =     70;         %Wing loading at the start of the cruise
WSec =     50;         %Wing loading at the end of the cruise

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, Vinf, WSsc, WSec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General inputs
A=0.002; %aspect ratio 
sweepc2 =0.25; %1/2 cord sweep
M=0; %zero because of uncompressability
Beta=sqrt(1-M^2);% 92.6m/s
eta=0.95;
%%Airfoil ans wing
sweepLE=0.25; %Calculate
alpha0L=0; %Given by the airfoil
%CLdes=0.65; %Given
taper=0.6; %Given
Sharpfactor=0.02; %from the airfoil.
%%statistical
%general
C1=0.02; %From Table (slide:15)
C2=0.02;  %From table slide: 15 raymer.

CLmax_base=0.2;
%% Use the datcom method to get the CL slope from statistics
Datcomtop=2*pi*A; 
Datcombottom=sqrt(4+(A*Beta/eta)^2*(1+((tan(sweepc2))^2/Beta^2)))+2;
CLalpha=Datcomtop/Datcombottom;  

%next is the trim angle.this is the angle the aircraft needs to fly at to
%fly at CLdes

alphatrim=CLdes/CLalpha+alpha0L;

%Next is the CLmax, 2 methods are being introduced. The Datcom method and a
%general  one. Datcom is prefferred. CLmax

Datcom_choose_method=((C1+1)*cos(sweepLE));

delta_CLmax=0.5; %Term for M>0.2. Not for take-off and landing.
delta_alpha_CLmax=0.02; %From table slide 19 (Raymer)
CLM_clmax=0.05; %slide 17 use right table 
CLmax_base=0.2; %table from Raymer
alpha_CLmax=0.04;

[CLmax, alpha_stall] = clalpha(A, Datcom_choose_method,CLM_clmax, delta_CLmax ,CLmax_base, alpha_CLmax, delta_alpha_CLmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ADSEE II - Lecture 3 - Drag coefficient estimation
%clc;

%% Inputs
% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound
v = 92.6; % [m/s], 180kts

cl = 0.30647;


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
FF_cs = [ADSEE_II_Drag.form_factor(option(1), L2, D, tc_avg, xc_max, sweep, v, a), ADSEE_II_Drag.form_factor(option(2), L2, D, tc_avg, xc_max, sweep, v, a), ADSEE_II_Drag.form_factor(option(3), L2, D, tc_avg, xc_max, sweep, v, a), ADSEE_II_Drag.form_factor(option(4), L2, D, tc_avg, xc_max, sweep, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = ADSEE_II_Drag.S_wet_c(S_ref, S_h, S_v, D, L1, L2, L3);


cd0_c = ADSEE_II_Drag.tot_comp_drag0(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0);
misc = ADSEE_II_Drag.cD_misc0(0.034, A_cs, L2*D*0.1, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c, c);

total_cD0 = cd0_c + misc

cD = total_cD0 + ADSEE_II_Drag.k_f(A, sweep, cl) * cl^2


