function [weights, Wl_Wto, Cd0, LD_cruise, W4W5, m_cruise, V_stall] = Fuel_Frac(chosen_MTOW, A, e, V_cruise, V_stall)
%% Implement Oliview code for weight estimation
clc;
% chosen_MTOW = 1835 % taken from the average of reference aircraft
mass_crew = 90.72; % kg
% Result for OEW, trendline param, trendline param, average_MTOW
[OEW_result, Slope, Intercept, average_MTOW] = ref_input_func(chosen_MTOW);
Mres_Mused = 0.25; % Assumption I guess

% A = 8; % <------ INPUT
% e = 0.75; % <------ INPUT

%% Zero lift drag coefficient and LD_cruise
C_fe = 0.0055; % <---- Assumed
Sw_S = 4; % <----- Assumed
Cd0 = C_fe*Sw_S;



LD_cruise = sqrt((pi*A*e)./(4*Cd0))
cp = 0.0000000845; % <----- Assumed
range = 200; % nmi  <----- INPUT
range = 1845 * range; % m
eta_p = 0.8;
g = 9.80065; 

M_4 = 1./exp(range/(LD_cruise.*eta_p./(g.*cp)));


%% Endurance equations
% V_cruise = 180; % kts <---- Input 
V_cruise = V_cruise * 0.51444;
V_loiter = 0.75 * V_cruise;
% V_stall = 61; % kts <---- Input
V_stall = 0.51444* V_stall;

LD_loiter = sqrt(3*Cd0*pi*A*e)./(4*Cd0);

E = 45; % min loiter by requirements
E = 60 * E; % [s]

M_6 = 1./exp(E*g*cp*V_loiter./(LD_loiter));

%% Rest fractions
M_1a = 0.995;
M_1b = 0.997;
M_2 = 0.998;
M_3 = 0.992;
M_5 = 0.993;

M_7 = 0.993;

M_ff = M_1a * M_1b * M_2 * M_3 * M_4 * M_5 * M_6 * M_7;

M_used = 1 - M_ff;

Wl_Wto = M_2 * M_3 * M_5 * M_7 * M_4 * M_6;

OEW_result = OEW_result + mass_crew;    % final OEW with crew

w_fuel_used = M_used*chosen_MTOW
w_fuel_res = Mres_Mused*w_fuel_used
w_fuel_tot = w_fuel_used + w_fuel_res

% Added
m_cruise = chosen_MTOW * M_1a * M_1b * M_2 * M_3;
%%%% if I delete it, delete m_cruise in output and in weights summary

W4W5 = M_4/M_5;

weights = ["OEW: ",OEW_result;  % kg
    "MTOW chosen: ",chosen_MTOW;% kg
    "Fuel used: ",w_fuel_used;  % kg
    "Aspect Ratio: ",A;
    "Oswald's factor: ",e;
    "Stall speed: ", V_stall;
    "Cruise speed: ", V_cruise;
    "Mass cruise: ", m_cruise;
    "Total fuel mass: ", w_fuel_tot]


end
%% WINGLOADING DIAGRAM!!!!

