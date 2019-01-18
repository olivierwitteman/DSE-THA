%S_v  vertical tail surface
%b_v  vertical tail span
%A_v vertical aspect ratio
%lambda_v Taper ratio vertical tail
%LAMDA_qc quarter chord sweep of the vertical tail
%trailingedgesweep_v sweep of the vertical tail at the trailing edge
%delta_engine rudder deflection needed for 1 critical engine
%delta_crosswind rudder deflection needed for landing with crosswind
%delta_spin rudder deflection needed to take account for spin (probably
%incorrect and just ignore this value, delete it if you want)
%c_v_root the root chord of the vertical tail
%c_v_tip the tip chord of the vertical tail
%c_v_mac the MAC of the vertical tail
%V_stall stall speed (mind the captital V)
%l_v distance from nose to vertical tail
%l_f Fuselage length
%S_W surface are of the wing
%v_cruise cruise speed (mind the lowercase v)
%fuselage_width fuselage width
%fuselage_height fuselage height
%T_cruise temperature at cruise altitude
%c_l_alpha the change in lift coefficient with respect to the angle of
%attack
%b wing span (of the main wing)
%x_m length from nose to c.g.
%rho_cruise air density at cruise altitude
%C_n_beta The change in moment N coefficient with respect to the sideslip

function [S_v, b_v, A_v, lambda_v,LAMBDA_qc, trailingedgesweep_v,...
    delta_engine, delta_crosswind, delta_spin, c_v_root, c_v_tip, c_v_mac]...
    = rudder(V_stall, l_v, l_f, S_W, v_cruise,fuselage_width, fuselage_height,...
    T_cruise, c_l_alpha, b, x_m, rho_cruise, C_n_beta)
%V_stall = 31.4;
rho = 1.225;
%l_v = 8.2; %distance to vertical tail
%S_W = 12;
%l_f = 8.2 ;%fuselage length
%v_cruise = 90;
kinematic_viscosity = 14.61; % 14.61
d_f = (fuselage_width + fuselage_height)/2; % 1.416 en 1.65 fuselage diameter
A_v = 1.5; %Sadraey, raymer
a = sqrt(1.4*T_cruise*287); %290
M = v_cruise/a;
beta = sqrt(1-M^2);
%c_l_alpha = 0.32;
        %c_l_alpha = 1.05/(beta)*(0.82*c_l_alpha1);
k = c_l_alpha/(2*pi/beta);
        %b_v =1.6854;
%b = 11.2891;
b_w = b;
V_MC = 1.2 * V_stall;

T_TO = 180000/V_MC ;%N
n_E = 2 ;%Number of engines
y_E = b/2 ;%The distance between the failed engine and the plane of symmetry
N_E = T_TO/n_E*y_E; %Moment caused by active engines
N_D = 0.75*N_E;     %Drag caused by failed engines  25
delta_F = degtorad(25); %25 or less
LAMBDA_qc = degtorad(35); %sweep (Sadreay p.321) (most aircraft DO have sweep)
lambda_v = 0.5; %Sadraey, Raymer
trailingedgesweep_v = atan(tan(LAMBDA_qc)-4/A_v*((1-0.25)*(1-lambda_v)/(1+lambda_v)));
LAMBDA_hc = atan(tan(LAMBDA_qc)-4/A_v*((0.5-0.25)*(1-lambda_v)/(1+lambda_v))); %sweep @ half chord
K_LAMBDA = (1-0.08*cos(LAMBDA_qc)^2)*cos(LAMBDA_qc)^(3/4) ;%Empirical correction factor (DATCOM 1978)
c_L_delta_th = 5 ;%c_r/c_v~cessna 182 = 0.42 &S_r/S_v = 0.38 and NACA 63-015A airfoil (t/c=0.15) Look in graph 
...Look at Fig 11.16 if rudder deflection is 25 degrees or smaller, if higher look in Fig 11.17 
...@http://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_11_EmpennageSizing.pdf
K = 0.7;%Dependent on the flap deflection, look in graph 11.17
FUNC = 0.85;
%x_m = 0.35*l_f; %length from nose to c.g.)
k_N = 0.01*(0.27*x_m/l_f-0.168*log(l_f/d_f)+0.416)-0.0005;
Re = v_cruise*l_f/kinematic_viscosity;
k_ri = 0.46*log10(Re/10^6)+1;
C_NbetaF = -360/(2*pi)*k_N*k_ri*l_f^2*d_f/(S_W*b);
C_L_alpha = 2*pi*A_v/(2+sqrt(A_v^2*beta^2/k^2*(1+tan(LAMBDA_hc)^2/beta^2)+4));
%C_L_alpha = 2*pi*A_v/(2+sqrt(A_v^2*(1+tan(LAMBDA_hc)^2-M^2)+4));
C_YbetaV = -C_L_alpha;


%Control requirement
N_V = 0.5*rho*V_MC^2*delta_F*FUNC*c_L_delta_th*K*K_LAMBDA*l_v;
S_V1 = (N_E+N_D)/N_V;

%Stability requirement
C_Nbeta = 0.057; %Or larger
%C_Nbeta = C_NbetaF - C_YbetaV * S_V2 * l_v/(S_W*b_w);
S_V2 = abs(S_W*(((C_Nbeta - C_NbetaF)/(-C_YbetaV))*(b_w/l_v)));


%THIS SIZING HAS TO FULFILL THE 'SPIN REQUIREMENT' (I.E. ABLE TO RECOVER
%FROM A WHOLE SPIN OR A 3 SECOND LONG SPIN (WHICHEVER TAKES LONGER) BUT NO
%CALCULATION HOW TO SIZE THE VERTICAL TAIL WING FOR SPIN CAN BE FOUND.
S_ref = S_W;
S_v = max(S_V1,S_V2);
    %l_vt = 8.; %Distance between aircraft c.g. and aerodynamic centre of the tail
    %l_v = l_vt; %distance between wing aerodynamic centre and tail aerodynamic centre
%K_f1 = 0.75; %Parameter that takes into account the fuselage directional stability (Sadraey p.322)
C_L_alpha_v = 0.4997; %vertical lift curve slope %Lift coefficient units per RADIANS!  %Airfoil 63-015A
eta_v = .95; %Dynamic pressure ratio
%V_v = 0.12; %Often between 0.02 and 0.12 (Sadraey p.320) (Assume Fokker F-27)
V_v = l_v*S_v/(b*S_ref);
    %l_h =8.; %Value from Sumant
    %l_vt = l_v;
%A_v = 1.84; %Sadraey, Raymer

%b = 11.2;
%S_v = b*S_ref*V_v/l_v;
%rho = 1.225;
%V = v_cruise;
%rho_cruise = 0.966632;
    %C_L_v = 0.;
    %C_n = 0.183; % flight dynamics program (martin)
%C_n_beta = -3.283e-3;%flight dynamics program (martin)
    %C_n_p = -2.151e-4;%flight dynamics program (martin)
    %C_n_r = -5.16e-5;%flight dynamics program (martin)
    %vertical_tail_sidewash = .1;
    %C_n_beta_v = K_f1*C_L_alpha_v*(1-vertical_tail_sidewash)*eta_v*V_v/l_v*l_vt; %This has to be positive to be stable


%V_v = l_v*S_v/(b*S_ref); %should be between 0.02 and 0.12
%L_v = 0.5*rho*V^2*S_v*C_L_v; %Lift of the vertical tail
%N_cg = L_v*l_v; %should equal to zero


b_v = sqrt(A_v*S_v); %A_v = Aspect ratio vertical tail
c_v_mac = S_v/b_v; %vertical tail chord at MAC
c_v_root = c_v_mac/((1+lambda_v+lambda_v^2)/(1+lambda_v))*(3/2);
c_v_tip = lambda_v*c_v_root;

    %c_r_over_c_v = 0.4;%cessna 182
tau_r = 0.55; %Rudder effectiveness (rudder equals (max?) 35% of tailplane (slides))
C_n_dR = -c_L_delta_th*V_v*eta_v*tau_r*b_v;
C_y_dR = C_L_alpha_v*eta_v*tau_r*S_v/S_W;

%Engine faillure req
T_max = T_TO/n_E;% Thrust CFR 23.147  
    %C_L_v = C_L_alpha * delta_F;
    %q_cruise = 0.5*rho_cruise*v_cruise^2;
    %L_V_max =  S_v*q_cruise*C_L_v;
    %N_A = q_cruise * S_ref * C_n * b;
delta_engine = T_max/(-rho_cruise*(1.4*V_stall)^2*S_ref*C_n_dR);
delta_engine = radtodeg(delta_engine);

%Crosswind req.
V_landing = 1.1*V_stall;
S_S = 11.7; %Value from Stan
C_D_y = 0.6; %Typically between 0.5 and 0.8
V_w = max(0.2 * V_stall, 20*0.51444444444444444444);%assume cross wind of 20 knots
V_T = sqrt(V_landing^2+V_w^2);
sideslip_crosswind = atan(V_w/V_landing);
F_w = 0.5*rho*V_w^2*S_S*C_D_y;
d_c = 1.8;% distance between c.g. and centre of the projected side area of the aircraft
syms sigma1 delta_crosswind
eqns = [0.5*rho*V_T^2*S_W*b*(C_n_beta*(sideslip_crosswind-sigma1)+C_n_dR*delta_crosswind)+F_w*d_c*cos(sigma1)==0, F_w==0.5*rho*V_T^2*S_W*(C_YbetaV*(sideslip_crosswind-sigma1)+C_y_dR*delta_crosswind)];
solvation = solve(eqns, [sigma1 delta_crosswind]);
delta_crosswind = radtodeg(solvation.delta_crosswind);

%Spin req.
R_SR = degtorad(240/3); %CFR
alpha_spin = degtorad(45); %typical value (sadraey p.700)
I_XXB = 1523.706;
I_ZZB = 2033.9;
I_XZB = 232.4;
I = [cos(alpha_spin)^2 sin(alpha_spin)^2 -sin(2*alpha_spin);
sin(alpha_spin)^2 cos(alpha_spin)^2 sin(2*alpha_spin);
0.5*sin(2*alpha_spin) -0.5*sin(2*alpha_spin) cos(2*alpha_spin)];
IB = [I_XXB; I_ZZB; I_XZB];
IB_solved = linsolve(I,IB);
N_SR = ((IB_solved(1)*IB_solved(2)-IB_solved(3)^2)/IB_solved(1))*R_SR;
S_Ve = S_v - (0.2*b_v*c_v_mac); %with 0.2 the percentage of area that isnt effective
V_Ve =  l_v*S_Ve/(b*S_ref);
C_n_dR_e = -C_L_alpha_v*V_Ve*eta_v*tau_r;
delta_spin = radtodeg(2*N_SR/(rho_cruise*V_stall^2*S_ref*b*C_n_dR_e));

end