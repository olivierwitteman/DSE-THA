                                                           
plot([WS_stall, WS_stall], [0, 2], "magenta")
hold on
% Higher Cl_max for landing shifts the magenta line to the right


%% Sizing for Take-Off
cl_to = [1.9, 2.0, 2.1]/1.21; % <-------INPUT cl for takeoff choices

s_to = 762; % m from project guide requirements


for i = cl_to
    PW_TO = top./WS_x*i*sigma;  % slide 36
    ylim([0, 1.5])
    plot(PW_TO, "red")
end
% The lowest Cl_to value represents the lowest line on the graph!

%% Sizing for landing
s_land = 762;   % from project requirements
V_land = 31.38; % This is in accordance to the regulation from CFR 23,


% f = Wl/Wto based on each group weight estimation
f_group1 = 0.8704 % <-------INPUT
f_group2 = 0.87   % <-------INPUT
f_group3 = 0.93322    % <-------INPUT
Cl_land_group3 = [2.0, 2.1, 2.2]   % <-------INPUT copied from slide 11 - single piston
Cl_land_group2 = [1.6, 1.8, 2.0, 2.2]   % <-------INPUT copied from slide 11 - single piston
Cl_land_group1 = [1.6, 1.9, 2.2, 2.5]   % <-------INPUT copied from slide 11 - twin engine


for i = Cl_land_group3 % CHANGE CL_LAND_GROUP# ACCORDING TO YOUR GROUP
    WS_land = i*rho_0*s_land/(0.5915)/(2*f_group3);   % CHANGE f_group# to your group
    plot([WS_land, WS_land], [0, 2], "green")
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


V_cruise = 92.6; % m/s from requirements project guide
eta_p = 0.8; % Assumption from slides

for i = A
    WPto = power_setting./per_mtow.*eta_p.*(rho_isa./rho_0).^(3/4)*...
        ((cd0_clean.*0.5.*rho_isa.*V_cruise.^(3)./(per_mtow.*WS_x) + ...
        per_mtow.*WS_x./(pi.*i.*e_clean.*0.5.*rho_isa.*V_cruise)).^(-1));
    ylim([0, 0.9])
    plot(WPto, "black")
end
% The lowest A represents the lowest black line in the graph

%% Sizing for climb rate performance %%%% XXX

for i = A
    
    WP_climb = eta_p./(sqrt(WS_x).*sqrt(2./rho_isa)./...
        (1.345.*(i.*e_clean).^(0.75)./(cd0_clean.^(1./4))));
    plot(WP_climb, "blue")
end
% The lower the A the lower is the blue line in the grph

%% Sizing for climb gradient


% https://www.law.cornell.edu/cfr/text/14/23.2120 website for climb grad
% 8.3% is also given in slide 93
grad = 0.083; 


for i = A
    cd = 4*cd0_clean;   % slide 77
    cl_grad = sqrt(3*cd0_clean*pi*i*e_clean);   % slide 77

    WP_grad = eta_p./(sqrt(WS_x).*(grad+cd./cl_grad).*sqrt(2./(rho_isa.*cl_grad)));
    plot(WP_grad, "cyan")
end
% The lower the A the lower is the cyan line for gradient

%% Reference AC

scatter([9.81*112, 9.81*125], [9.81/130, 9.81/150]);

% design_point = 1190, 0.052119

scatter(1190, 0.052119)



%% Legend code, if you change your Cls, Cds or A just change the values here as well
set(gca, "FontSize", 40)
ylim([0, 0.2])
xlabel("W/S [N/m^2]",'FontSize', 20)
ylabel("W/P [N/W]",'FontSize', 20)
legend({'Stall land. config.'...
    ,'Cl_{to} = 1.9','Cl_{to} = 2.0', 'Cl_{to} = 2.1'...
    ,'Cl_{la} = 2.0','Cl_{la} = 2.1','Cl_{la} = 2.2'...
    ,'V_{cruise} at A = 8'...
    ,'c at A = 8'...
    ,'c/V at A = 8'...
    ,'Reference Aircraft', 'Design point'},'Location','northeast', 'FontSize', 16)

