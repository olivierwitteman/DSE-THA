function [S_h, S_v] = control_surf_func(MAC, S, b, l_fus, empen_x, x_aft_cg)
%% Procedure explained from ADSEE I lecture 7 slides 54
% MAC = 1.5;              % <------ INPUT
% 
% S = 15.3;               % <------ INPUT
% 
% b = 10;                 % <------ INPUT
% 
% l_fus = 6.6;            % <------ INPUT
% 
% empen_x = 0.8;          % <------ INPUT
empen_x = empen_x * l_fus
% 
% x_aft_cg = 2.7;         % <------ INPUT


V_h = [0.92, 0.60, 0.71, 0.61, 0.48, 0.7,...
    0.61, 0.76, 0.49, 0.83, 0.63];                    % from slides (source Roskam)
V_h = mean(V_h);


V_v = [0.0466, 0.038, 0.047, 0.037, ...
    0.026, 0.038, 0.037, 0.024, 0.039, 0.086, 0.062]; % from slides (source Roskam)
V_v = mean(V_v);

X_h = empen_x + 0.05*l_fus;
X_v = X_h;


S_h = (V_h * S * MAC)/(X_h - x_aft_cg)
S_v = (V_v * S * b)/(X_v - x_aft_cg)
end






