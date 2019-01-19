clc

% CZ_q
Cmq = 0.0005;       % < -------- INPUT
MAC = 1.14;         % < -------- INPUT
lh = 3.9;           % < -------- INPUT

CZ_qh = Cmq * MAC / lh * 1 / 1.1;
CZ_q = 2 * CZ_qh

% CZ_alpha and CZ_alphad
CL_alpha = 0.005;   % < -------- INPUT
CL_alphad = 0.0004; % < -------- INPUT

CZ_alpha = - CL_alpha; % page 170
CZ_alphad = CL_alphad;

% PHUGOID
i_p = 3*180/pi; % incidence angle
alpha_0 = 0.2
CX_u = 0; % -3*Cd*(1 - 1) cuz Cd is equal to Tc
CZ_u = - 2 * CL + CD * (-(alpha_0 + i_p) + 3*1);
CZ_0 = - CL - CD*(alpha_0 + i_p)
CM_u = - 3 * Cd * zv;s
