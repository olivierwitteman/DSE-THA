%% ADSEE II - Lecture 3 - Drag coefficient estimation
clc

%% Inputs
% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound
v = 92.6; % [m/s], 180kts

cambered = 0; % 1 for True, 0 for False
A = 8;
sweep = 0;
e = 0.75; % default assumed Oswald factor
S_ref = 15;
S_h = 0.15*S_ref;
S_v = 0.1*S_ref;
tc_avg = 0.15; % (t/c)_avg is the average thickness to chord
xc_max = 0.25; % (x/c)_max is the position of maximum thickness

C_f_e = 0.0055; % light AC - single engine
% C_f_e = 0.0045; % light AC - twin engine

S_W = 4 * S; % assumed wetted area

% k = 0.152E-5; % polished sheet metal
k = 0.634E-5; % smooth paint
% k = 0.052E-5; % smooth molded composite

L1 = 1; % nosecone length
L2 = 4; % main fuselage length
L3 = 2; % tailcone length
D = 1; % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)

%% Component contributions propeller AC (fast method)
Component = ['Wing', 'Fuselage multi-engine', 'Fuselage single-engine', 'Nacelles', 'Tail (hor + ver)', 'misc'];
C_D_Cs = [0.007, 0.08, 0.11, 0.06, 0.008, 0.15];
% Change these according to component name (defined above)
A_Cs = [15, 0, 3.0, 0.1, 0.3];


Fast_Cd0 = fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref)

%% Component method
% [Fuselage, Wing, horizontal tail, vertical tail]
C_f_c_fuselage = fp_skin_friction(0.1, k, rho, v, L2, mu, a);
C_f_c_wingtail = fp_skin_friction(0.4, k, rho, v, L3, mu, a);
C_f_cs = [C_f_c_fuselage, C_f_c_wingtail, C_f_c_wingtail, C_f_c_wingtail];
option = [2, 1, 1, 1];
FF_cs = [form_factor(option(1), L2, D, tc_avg, xc_max, sweep, v, a), form_factor(option(2), L2, D, tc_avg, xc_max, sweep, v, a), form_factor(option(3), L2, D, tc_avg, xc_max, sweep, v, a), form_factor(option(4), L2, D, tc_avg, xc_max, sweep, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = S_wet_c(S_ref, S_h, S_v, D, L1, L2, L3);


tot_comp_drag(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0)

%% Functions

function misc = cD_misc(A_cs_max, A_base, v, a, S_A, S_r)
% Fuselage upsweep (A_max is max fuselage cross sectional area, u is upsweep in rad)
dcD = 3.83*A_cs_max*u^2.5;
% Fuselage base drag
dcD = [dcD, (0.139 + 0.419*(v/a - 0.161)^2)*A_base];
% Landing gear (nose)
m = S_A/S_r;
n = S_A/A_cs_max;
dcD = [dcD, n*0.3];

% Landing gear (main)
c = 0.05328; % open wheel wells
% c = 0.04955; % closed wheel wells
dcD = [dcD, n*c*exp(5.615*m)];


misc = sum(dcD);

function cD0 = fast_sum_C_D_0(C_D_Cs, A_Cs, S)
iterated = 0;
for i = 1:length(C_D_Cs)-1
    iterated = [iterated, dCDzero(A_Cs(i), C_D_Cs(i), S)];
end
cD0 = (1 + C_D_Cs(end)) * sum(iterated);
end

function cD_c = comp_drag(C_f_c, FF_c, IF_c, S_c)
cD_c = C_f_c * FF_c * IF_c * S_c;
end

function ff = form_factor(type, L2, D, tc_avg, xc_m, sweep_m, v, a)
f = L2/D;
if type == 1
    % Wing, tail, strut and pylon
    ff = (1 + 0.6/xc_m * tc_avg + 100*tc_avg^4)*(1.34*(v/a)^0.18 * cos(sweep_m)^0.28);

elseif type == 2
    % Fuselage and smooth canopy
    ff = 1 + 60/f^3 + f/400;

elseif type == 3
    % Nacelle and smooth external store
    ff = 1 + 0.35/f;
else
    print('Not type given, ff=0')
    ff = 0; 
end
end

function cf = fp_skin_friction(laminar, k, rho, v, l, mu, a)
Re_sub = min(rho*v*l/mu, 38.21*(l/k)^1.053);
lam = 1.328/sqrt(Re_sub);
turb = 0.455/(log10(Re_sub)^2.58*(1+0.144*(v/a)^2)^0.65);

cf = laminar*lam + (1-laminar)*turb;
end

function s_c = S_wet_c(S_w_exp, S_HT_exp, S_VT_exp, D, L1, L2, L3)
s_c = [S_F_wet(D, L1, L2, L3), S_w_wet(S_w_exp), S_HT_wet(S_HT_exp), S_VT_wet(S_VT_exp)];
end

function cD0 = tot_comp_drag(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, CD0_misc)
cd0 = 0;
for j = 1:length(C_f_cs)
    cd0 = [cd0, comp_drag(C_f_cs(j), FF_cs(j), IF_cs(j), S_cs(j))];
end
cD0 = 1/S_ref * sum(cd0) + CD0_misc;
end

function dC_D_0 = dCDzero(A_c, C_D_c, S)
dC_D_0 = 1/S * C_D_c * A_c;
end

function s = S_w_wet(S_w_exp)
s = 1.07*2*S_w_exp;
end

function s = S_HT_wet(S_HT_exp)
s = 1.05*2*S_HT_exp;
end

function s = S_VT_wet(S_VT_exp)
s = 1.05*2*S_VT_exp;
end

function s = S_F_wet(D, L1, L2, L3)
s = pi*D/4 * ((1/(3*L1^2) * (4*L1^2 + 0.25*D^2)^1.5 - 1/8 * D^3) - D + 4*L2 + 2*sqrt(L3^2 + D^2 /4));
end



%% To be implemented
% % Standard C_D calculation
% if cambered == True
%     C_D = C_D_min + K * (C_L - C_L_minD);
% else
%     C_D = C_D_0 + K * C_L * C_L;
% end

% C_D_0 = C_f_e * S_W/S;

% % Calculated inputs:
% K = 1/(pi * A * e);