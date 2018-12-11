function [weights, Cd0, LD_cruise, M_batt, MTOW, OEW] = Electric_range_equation (chosen_MTOW, A, e, V_cruise, V_stall, Payload, OEW, cb)
%% Implement Oliview code for weight estimationeference aircraft
mass_crew = 90.72; % kg
% Result for OEW, trendline param, trendline param, average_MTOW
[OEW_result, Slope, Intercept, average_MTOW] = ref_input_func(chosen_MTOW);
Mres_Mused = 0.25; % Assumption I guess

% A = 8; % <------ INPUT
% e = 0.75; % <------ INPUT
MTOW=chosen_MTOW
%% Zero lift drag coefficient and LD_cruise
C_fe = 0.0055; % <---- Assumed
Sw_S = 4; % <----- Assumed
clc;
% chosen_MTOW = 1835 % taken from the average of r
Cd0 = C_fe*Sw_S;

LD_cruise = sqrt((pi*A*e)./(4*Cd0))

range = 250; % nmi
range = 1845 * range; % m
eta_p = 0.8;
g = 9.80665; 
eta_i=0.95;%inverter efficiency
eta_m=0.95;%motor efficiency
eta_prop=0.8; %propellor efficiency
eta_total=eta_i*eta_m*eta_prop;
M_batt_ratio = range/eta_i/eta_m/eta_prop*g/cb/LD_cruise
M_batt=M_batt_ratio*MTOW

%% Endurance equations

%% Rest fractions

newMTOW=OEW+M_batt+Payload
disp(abs(newMTOW-MTOW))
while abs(newMTOW-MTOW)>0.1
    MTOW=newMTOW
    M_batt=1.2* M_batt_ratio*MTOW;
    OEW=OEW   +M_batt+Payload
    newMTOW=OEW+Payload+M_batt
   
end  
%OEW=M_batt+OEW
disp(MTOW)
%%%% if I delete it, delete m_cruise in output and in weights summary

weights = ["OEW: ",OEW_result;  % kg
    "MTOW chosen: ",MTOW;% kg
    "Battery Mass ", M_batt;  % kg
    "Aspect Ratio: ",A;
    "Oswald's factor: ",e;
    "Stall speed: ", V_stall;
    "Cruise speed: ", V_cruise;
    ]


end
%% WINGLOADING DIAGRAM!!!!

