function [summary_wing] = wing_planform_design(V_cruise, A, S, h, MTOW)

thick = 6;  % used for graphs

% ISA calculation
lambda = 0.0065; 
T_0 = 288.15; 
g = 9.81;
R = 287.1;
P0=101.325*10.^3;
T = T_0 - lambda * h;
P = P0 * (T / T_0).^(g / (lambda * R));
rho_isa = P / (R * T);

gamma = 1.4;

a_cruise = sqrt(gamma*R*T);
M_cruise = V_cruise./a_cruise;

M_plus = 0.935;
M_dd = M_plus + 0.935;

% Sweep
if M_cruise < 0.7    
    sweep_c4 = acos(1);
end
if M_cruise > 0.7
    sweep_c4 = acos(0.75.*M_plus/M_dd);
end

% Taper
taper = 0.2*(2 - sweep_c4);  % from slide 16

b = sqrt(S*A);

cr = 2*S/((1+taper).*b);
ct = taper*cr;

% half chord sweep
sweep_c2 = atan(tan(sweep_c4) - 4/A*1/4*((1-taper)/(1+taper)));
sweep_c2_deg = sweep_c2*180/pi;


% thickeness to chord
cl_cruise = MTOW*9.81*2*1./(rho_isa*V_cruise*V_cruise*S);
tc = 0.18;

% MAC
MAC = cr * 2/3 * ((1 + taper + taper.^(2))/(1+taper));
MAC_pos = b./2 *(cr - MAC)./(cr - ct);

% dihedral
dihedral = 3 - sweep_c4/10;

dsweep_dc = sweep_c4 - sweep_c2;
sweep_LE = dsweep_dc + sweep_c4;
sweep_TE = sweep_c4 - 3 * dsweep_dc;

summary_wing = ["M_cruise",M_cruise;
    "Span", b;
    "Taper", taper;
    "Sweep c_4 (rad): ", sweep_c4;
    "Sweep c_2 (rad): ", sweep_c2;
    "Sweep LE (rad): ", sweep_LE;
    "Sweep TE (rad): ", sweep_TE;
    "Cl cruise: ", cl_cruise;
    "c root", cr;
    "c tip", ct;
    "MAC: ", MAC;
    "MAC position", MAC_pos;
    "Dihedral", dihedral;
    "t/c", tc;
    ];

% f4 = figure;
% plot([0, b/2], [5, 5], "b-",'LineWidth',thick)  % span
% hold on
% plot([b/2, b/2], [5-cr, 5], "b-",'LineWidth',thick)  % cr
% plot([0, 0], [5-ct, 5], "b-",'LineWidth',thick)  % ct
% plot([0, b/2] , [5-ct, 5-cr], "b-",'LineWidth',thick)  % TE

end


