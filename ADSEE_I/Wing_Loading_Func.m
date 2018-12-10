function [] = Wing_Loading_Func(h, A ,e_clean,cd0_clean, Cl_max, cl_to ,c, Wl_Wto, V_land, V_stall, V_cruise)
clc;
% The goal is to pick a point at the top right corner. This will give
% information on what S, A, Cl and W/P you need.
thick = 5;
rho_0 = 1.225;

%% sizing for stall
% Cl_max so for all flaps and HLDs
WS_stall = 0.5*rho_0*V_stall*V_stall*Cl_max;     

f2 = figure; %%%%%%%%%%%%%
plot([WS_stall, WS_stall], [0, 2], "magenta",'LineWidth',thick);
hold on
% Higher Cl_max for landing shifts the magenta line to the right


%% Sizing for Take-Off


s_to = 762; % m from project guide requirements

% How to find the TOP:
% read off your value for 2500 ft S_to and then multiply it by 0.2512!
top = 250*0.2512 % checked with formula on slide 41

sigma = 1; % rho/rho0 slide 27 . Assume from Sizing mission - sea level
WS_x = [1:1:1800];


for i = cl_to
    PW_TO = top./WS_x*i*sigma;  % slide 36
    ylim([0, 1.5])
    plot(PW_TO, "red",'LineWidth',thick);
end
% The lowest Cl_to value represents the lowest line on the graph!

%% Sizing for landing
s_land = 762;   % from project requirements or CFR reuq?


for i = Cl_max % CHANGE CL_LAND_GROUP# ACCORDING TO YOUR GROUP
    WS_land = i*rho_0*s_land/(0.5915)/(2*Wl_Wto);   % CHANGE f_group# to your group
    plot([WS_land, WS_land], [0, 2], "green",'LineWidth',thick);
end
% the lowest cl_land is the lowest green vertical line

%% Cruise Speed Sizing
lambda = 0.0065; 
T_0 = 288.15; 
g = 9.81;
R = 287.1;
P0=101.325*10.^3;
T = T_0 - lambda * h;
P = P0 * (T / T_0)^(g / (lambda * R));
rho_isa = P / (R * T)

power_setting = 0.9; % assumed slide 66. % <-------INPUT
per_mtow = 1; % as asked in sizing mission


% V_cruise = 92.6; % m/s from requirements project guide

eta_p = 0.8; % Assumption from slides

for i = A
    WPto = power_setting./per_mtow.*eta_p.*(rho_isa./rho_0).^(3/4)*...
        ((cd0_clean.*0.5.*rho_isa.*V_cruise.^(3)./(per_mtow.*WS_x) + ...
        per_mtow.*WS_x./(pi.*i.*e_clean.*0.5.*rho_isa.*V_cruise)).^(-1));
    ylim([0, 0.9])
    plot(WPto, "black",'LineWidth',thick);
end
% The lowest A represents the lowest black line in the graph

%% Sizing for climb rate performance %%%% XXX
% e_clean, cd0_clean

for i = A
    
    WP_climb = eta_p./(sqrt(WS_x).*sqrt(2./rho_isa)./...
        (1.345.*(i.*e_clean).^(0.75)./(cd0_clean.^(1./4))));
    plot(WP_climb, "blue",'LineWidth',thick);
end
% The lower the A the lower is the blue line in the grph

%% Sizing for climb gradient
% e_clean, cd0_clean

% https://www.law.cornell.edu/cfr/text/14/23.2120 website for climb grad
% 8.3% is also given in slide 93
grad = 0.083; % <------- INPUT


for i = A
    cd = 4*cd0_clean;   % slide 77
    cl_grad = sqrt(3*cd0_clean*pi*i*e_clean);   % slide 77

    WP_grad = eta_p./(sqrt(WS_x).*(grad+cd./cl_grad).*sqrt(2./(rho_isa.*cl_grad)));
    plot(WP_grad, "cyan",'LineWidth',thick);
end
% The lower the A the lower is the cyan line for gradient

%% Reference AC

scatter([9.81*112, 9.81*125], [9.81/130, 9.81/150]);

%% Legend code, if you change your Cls, Cds or A just change the values here as well
set(gca, "FontSize", 20)
ylim([0, 0.2])
xlabel("W/S [N/m^2]",'FontSize', 20)
ylabel("W/P [N/W]",'FontSize', 20)


legend({'Stall land. config.'...
    ,'Cl_{to} = ' + string(cl_to)...
    ,'Cl_{la} = ' + string(Cl_max)...
    ,'V_{cruise} at A = ' + string(A)...
    ,'c at A = ' + string(A)...
    ,'c/V at A = ' + string(A)...
},'Location','northeast', 'FontSize', 16)
end


