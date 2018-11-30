%inputs (page 477-479 Raymer) (everything is in retard units) (lbs,
%gallons, ft^3, ft^2, inch etc.)
W_dg = 1.; %Design gross weight
N_z = 1.; %Load factor
N_gear = 1;
S_w = 1.; % Wing surface
A = 1.; %Aspect ratio
t_over_c = 1.; %Wing thickness ratio
lambda = 1.; %taper ratio
LAMBDA = 1.; %Sweep at 25% MAC
S_csw = 1.;
K_door = 1.;
K_lg = 1.;
L = 1.;
S_f = 1.; %Wetted area
K_ws = 1.;
L_over_D = 1.;
W_fw = 1.;
V_cruise = 1.;
rho = 1.;
q = 0.5*rho*V_cruise^2;
S_ht = 1.;% Horizontal tailwing surface
LAMBDA_ht = 1.; % Sweep at 25% MAC
A_ht = 1.; % Aspect ratio horizontal tailwing
lamba_h = 1.; %Taper ratio tail
H_t_over_H_v = 1.; % = 0 for conventional tail, 1 for 1 tail
S_vt = 1.; % Surface vertical tail
LAMBDA_vt = 1.; % Sweep at 25% of vertical tail MAC
A_vt = 1.; % Aspect ratio vertical tail
lambda_vt = 1.; % taper raio vertical tail
lambda_h = 1; %Taper ratio horizontal tail 
L_t = 1.; % Tail length, wing quarter MAC to tail quarter MAC
W_press = 0 ;%11.9+(V_pr*P_delta)^0.271; %Weight penalty due to pressurization; PROBABLY ZERO FOR OUR DESIGNS BECAUSE WE DON'T PRESSURIZE OUR CABIN
N_l = 1.5 * N_gear; %Ultimate alnding load factor
W_l = 1.; %Landing design gross weight
L_m = 1.; %Extended length of main landing gear
L_n = 1. ; %Extended nose gear length (inch)
W_en = 1.; %Engine weight (each) in pounds
N_en = 1; %Number of engines\
V_t = 1. ; %Total fuel volume in gallons
V_i = 1. ;%Integral tanks volume in gallons
N_t = 1; %Number of fuel tanks
L = 1. ; %Fuselage structural length in ft
B_w = 1.; %Wingspan
W_uav = 1.; %Uninstalled avionics weight in pounds
N_p = 1; %Number of personal onboard
M = 1.; %Mach number
W_wing = 0.036*S_w^0.758*W_fw^0.0035*(A/((cos(LAMBDA))^2.))^0.6*q^0.006*lambda^0.04*(100*t_over_c/cos(LAMBDA))^-0.3*(N_z*W_dg)^0.49
W_horizontaltail = 0.016*(N_z*W_dg)^0.414*q^0.168*S_ht^0.896*(100*t_over_c/cos(LAMBDA))^-0.12*(A_ht/(cos(LAMBDA_ht))^2)^0.043*lambda_h^-0.02
W_verticaltail = 0.073*(1+0.2*H_t_over_H_v)*(N_z*W_dg)^0.376*q^0.122*S_vt^0.873*(100*t_over_c/cos(LAMBDA_vt))^-0.49*(A_vt/(cos(LAMBDA_vt))^2)^0.357*lambda_vt^0.039
W_fuselage =0.052*S_f^1.086*(N_z*W_dg)^0.177*L_t^-0.051*L_over_D^-0.072*q^0.241+W_press
W_mainlandinggear = 0.095*(N_l*W_l)^0.768*(L_m/12)^0.409
W_noselandinggear = 0.125*(N_l*W_l)^0.566*(L_n/12)^0.845
W_installedengines = 2.575*W_en^0.922*N_en
W_fuelsystem = 2.49*V_t^0.726*(1/(1+V_i/V_t))^0.363*N_t^0.242*N_en^0.157
W_flightcontrols = 0.053*L^1.536*B_w^0.371*(N_z*W_dg*10^(-4))^0.80
W_hydraulics = 0.001*W_dg
W_avionics = 2.117*W_uav^0.933
W_electrical = 12.57*(W_fuelsystem+W_avionics)
W_airco_and_anti_ice = 0.265*W_dg^0.52*N_p^0.68*W_avionics^0.17*M^0.08
W_furnishings = 0.0582*W_dg - 65