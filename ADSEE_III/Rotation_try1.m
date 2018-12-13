clc;
vars = load('../ADSEE_I/variables_ADSEE_I.mat');
MAC = double(vars.MAC)
CL_max = 1.9;
CL_h = -0.6215;
Cm_ac = -0.3278;
V_s1 = 31.2; % m/s
V_r = 1.05*V_s1; % CHECK REGULATIONS
x_g = x_lemac + 0.4*MAC; % 40% of mac
zh = 1.08;  % ????
zt = 0.68; % from landing gear file
CL_R = 1.2; % assumed



lh = 3.8; % from previous file
x_h = 6.8;
Vh_Vr = 0.85;

nh = (x_h - x_g)/lh*(Vh_Vr)^(2)

CLa_h = 4.2;  % NO IDEA which1 it should be
theta = 0.05; % ASSUMED

nq = 1 + CLa_h/CL_h * theta * (x_h - x_g)/V_r

Sh_S_rot = (CL_max/(nh*nq*CL_h)*(Cm_ac/CL_max - (V_s1/V_r)^(2)*(x_g - zt)/MAC) + ...
    CL_r/CL_h*(x_g/MAC - 0.25))*MAC/lh




