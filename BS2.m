rho = 1.225;

% Assumptions :
%               - Single engine piston props
%               - Cl assumptions

Cl_max_clean = (1.3 + 1.9)/2;   % Slide 11
Cl_max_TO = (1.3 + 1.9)/2;  % Slide 11
Cl_max_LA = (1.6 + 2.3)/2;  % Slide 11

% Determine Drag polar
% determine delta cd and oswald factor in different configurations

dCd0_to_flaps = 0.015;    % slide 13
de_to_flaps = 0.05;       % slide 13
dCd0_LA_flaps = 0.065;    % slide 13
de_LA_flaps = 0.1;        % slide 13
dC0_undercarriage = 0.02; % slide 13
de_undercarriage = 0;     % slide 13

% drag coefficients (propeller)
cd0_clean = 0.0280;  % slide 15
e_clean = 0.78;

cd0_to_lg_down = 0.0380;
e_to_lg_down = 0.83;
Cl_max_lg_down = 2.0;

cd0_la_lg_down = 0.0730;
e_la_lg_down = 0.88;
Cl_max_lg_down = 2.4;


% sizing for stall
V_stall = 31.38; % m/s assumed from slide 18        CHECK FAA regualtion!!!
WS_stall = 0.5*rho*V_stall*V_stall*Cl_max_clean  % N/mˆ2
plot([WS_stall, WS_stall], [0, 2])
hold on

%% Sizing for Take-Off
% Cl taken 1.6, 1.9, 2.0

cl_to = [1.6, 1.9, 2.0]/1.21;  
s_to = 700; % m Assumed from slide 41 CHECK HOW TO DO IT CORRECTLY!!!!

%top = [0.0577, 8,6726, -700];
%roots(top)
top = 58.18; % assumed from slide 41 CHANGE LATER!!!!

sigma = 1; % rho/rho0 slide 27 . Assume from Sizing mission - sea level
WS_x = [1:1:1800];


for i = cl_to
    PW_TO = top./WS_x*i*sigma;  % slide 36
    ylim([0, 1.5])
    plot(PW_TO, "red")
end
%V_to = 1.1*V_stall

%% Sizing for landing
% Cl considered 1.6, 1.8, 2.0, 2.2
V_land = 33.4; % m/s assumed from website TO BE REVISITED!!!!
s_land = 0.5915*V_land.^2; % assumed slide 52!!!!


f = 0.874 % based on group 2 weight estimation
Cl_land = [1.6, 1.8, 2.0, 2.2]

for i = Cl_land
    WS_land = i*rho*s_land/(0.5915)/(2*f);
    plot([WS_land, WS_land], [0, 2], "green")
end


%% Cruise Speed Sizing
h = 2000; % m Assumed TO BE REVISITED !!!!
lambda = 0.0065; 
T_0 = 288.15; 
g = 9.81;
R = 287.1;
P0=101.325*10^3;
T = T_0 - lambda * h;
P = P0 * (T / T_0)^(g / (lambda * R));
rho_isa = P / (R * T)

power_setting = 0.9; % assumed slide 66 OUR DECISION/REGULATION!!! check!!!
per_mtow = 1; % as asked in sizing mission
A = [6, 9, 12]; % copied from slides 67 OUR DECISION I GUESS


V_cruise = 92.5; % taken form slides TO BE ESTIMATED FROM STATS ANALYSIS
eta_p = 0.8; % Assumed absolutely no check has been done TO BE REVISITED

for i = A

    WPto = power_setting./per_mtow.*eta_p.*(rho_isa./rho).^(3/4)*...
        ((cd0_clean.*0.5.*rho_isa.*V_cruise.^(3)./(per_mtow.*WS_x) + ...
        per_mtow.*WS_x./(pi.*i.*e_clean.*0.5.*rho_isa.*V_cruise)).^(-1));
    ylim([0, 0.9])
    plot(WPto, "black")
end


%% Sizing for climb rate performance
% [A] considered 6,9,12 from slides

c = 5; % As from slide 78 TO BE CHANGED RECORDING TO REQUIREMENTS
eta_p = 0.8; % As from slide 78

for i = A

    WP_climb = eta_p./(sqrt(WS_x).*sqrt(2./rho_isa)./(1.345.*(i.*e_clean).^(0.75)./(cd0_clean.^(1./4))));
    plot(WP_climb, "blue")
end


%% Sizing for climb gradient
%%%%% Not completely finished !!! Only the formula is put in. Additional
%%%%% line for changing Cl with changing A should be added!!!1

% A considered : 6, 9 , 12 
grad = 0.083; % taken from slide 93 TO BE CHECKED AGAIN WITH REGULATIONS!

cd = 0.05; % random should be calculated with drag polar!!!!!


WP_grad = eta_p./(sqrt(WS_x).*(grad+cd./Cl_max_clean).*sqrt(2./(rho_isa.*Cl_max_clean)));
plot(WP_grad, "cyan")
xlabel("WS")
ylabel("W/P")


legend({'Stall land. config.','Cl_to = 1.6','Cl_to = 1.9','Cl_to = 2.0'...
    ,'Cl_la = 1.6','Cl_la = 1.8','Cl_la = 2.0','Cl_la = 2.2'...
    ,'V_cruise at A = 6','V_cruise at A = 9','V_cruise at A = 12'...
    ,'c at A = 6','c at A = 9','c at A = 12'...
    ,'c/V at A = 6'},'Location','northeast')
title("W/P to W/S Diagram")

