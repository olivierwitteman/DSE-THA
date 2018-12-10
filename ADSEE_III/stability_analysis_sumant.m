
clear all
% Weights based on engine selection and preliminary design (might change)

vars = load('../ADSEE_I/variables_ADSEE_I.mat');
OEW = vars.OEW;  %operational empty weight [N]
MTOW = vars.MTOW; %maximum take-off weight [N]
W4W5 = vars.W4W5;
W_f = vars.W_fuel_used;
% clear all
% Weights based on engine selection and preliminary design (might change) 
OEW = 1300*9.81;  %operational empty weight [N]
MTOW = 1900*9.81; %maximum take-off weight [N]
Payload  = 350*9.81; %total payload [N]
%Dimensions and parameters (fixed) 
S = vars.S; %reference surface area [m^2]
A = vars.A;   %aspect ratio [-]
b = sqrt(S*A) %span [m]
cg = S/b; %average constant chord [m]
lambda = vars.double(tr) %taper ratio of main wing [-]
Cr = 2*S/((1+lambda).*b); %main wing root chord [m]
Ct = lambda*Cr; %main wing tip chord [m]
v = vars.V_cruise; %cruise velocity% h = 2400m
rho = 0.966632; % [kg/m^3]
T = 272.55; % T[K]
a = sqrt(1.4*287.15*T); % speed of sound % v = 92.6; % [m/s], 180kts
M = v/a; %mach number at cruise [-]
beta=sqrt(1-M^2);
bf = 1.5; %fuselage width [m] assumption
hf = 1.5; %fuselage height [m]assumption
bh = 4.0; %horizontal tail span [m]assumption
Sh = 2.726; %horizontal tail area [m^2]assumption
Sv = 1.32;  %vertical tail area [m^2]assumption

Cr_h = 0.8; %horizontal tail root chord [m] 
Ct_h = 0.563; %horizontal tail tip chord [m]
lfn = 2.059; %distance from nose to leading edge of root chord [m]
ln = 2.425; %distance from engine to quater chord mac [m]
bn = 0.9; %width of nacelles (engines) [m]
lambda_h = Ct_h/Cr_h; %taper ratio of horizontal tail [-]
Ah = bh^2/Sh; %aspect ratio of horizontal tail [-] 
Snet = S-bf*((Cr+1.35)/2); %net area (excluding the eclosed wing area in the fuselage) [m^2]
eta = 0.95; %airfoil efficiency factor [-]
lf = 7.75; % total length of fuselage [m]
l_press = 5.0; % length of pressurized area (assumed) [m]
%-------------------------------------------------------------------------
% Calculate wing sweep angle
sweep_LE = 5.68*pi/180; %sweep at leading edge [rad] %needs to be changed
sweep_14 = atan(tan(sweep_LE) + (Cr/(2*b))*(lambda -1)); %sweep at quater chord [rad]
sweep_12 = atan(tan(sweep_LE) - (4/A)*(0.5*((1-lambda)/(1+lambda)))); %sweep at half chord [rad]
%-------------------------------------------------------------------------
% Calculate horizontal tail sweep angle
sweep_LE_h = 6.27*pi/180; %sweep at leading edge [rad] %needs to be changed
sweep_14_h = atan(tan(sweep_LE_h) + (Cr_h/(2*bh))*(lambda_h - 1)); %sweep at quater chord [rad]
sweep_12_h = atan(tan(sweep_LE_h) - (4/Ah)*(0.5*((1-lambda_h)/(1+lambda_h)))); %sweep at half chord [rad]
%-------------------------------------------------------------------------
% Calculate wing mean aerodynamic chord (mac)
mac = (2/3)*Cr*( (1 + lambda + lambda^2)/(1+lambda)); % mean aerodynamic chord [m]
y_mac = (b/6)*((1 + 2*lambda)/(1 + lambda)); %y location of mac [m]
x_mac = y_mac*tan(sweep_LE);          %x location of mac [m] 
x_datum = 2.059 ;                  % measured from planform for given geometry (from nose to wing) [m]
%--------------------------------------------------------------------------
%Calculate horizontal tail mean aerodynamic chord
mac_h = (2/3)*Cr_h*( (1 + lambda_h + lambda_h^2)/(1+lambda_h)); % mean aerodynamic chord of the horizontal tail [m] 
y_mac_h = (bh/6)*((1 + 2*lambda_h)/(1 + lambda_h)); % y location of mac of horizontal tail [m]
x_mac_h = y_mac_h*tan(sweep_LE_h); %x location of mac of horizontal tail [m] 
x_datum_h =6.9 ;                  % measured from planform for given geometry (from nose to horizontal tail) [m]
%--------------------------------------------------------------------------
lh = x_datum_h + x_mac_h + 0.25*mac_h - (x_datum + x_mac + 0.25*mac); %distance between aerodynamic center of main wing and horizontal tail [m] 
x_lemac = x_datum + x_mac;          % x location of leading edge mean aerodynamic chord [m]
x_OEW = x_lemac + 0.375*mac;        % assumed CG of operational empty weight [m]
x_Cargo = 4.0;                       % assumed CG of cargo in meters [m]
x_Fuel = x_lemac + 0.5*mac;         % CG of fuel [m]
%--------------------------------------------------------------------------

%% scissor plot
%%% Variables Stability (CRUISE CONFIGURATION *- deg Aoa, - deg Incidence, V = m/s)
SM = 0.15;
%stability margin for safety given as percentage/100
CLa_w = (2*pi*A)/(2 + sqrt(4+(A/eta)^2*(1 +(tan(sweep_12)^2/beta^2)) )) ; 
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for main wing
CLa_h = (2*pi*Ah)/(2+ sqrt(4 + (Ah/eta)^2 *(1+(tan(sweep_12_h)^2/beta^2)) ));
%CLa_h = 2.0; %(for verification)
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for horizontal tail
CLa_A_h = CLa_w*(1+(2.15*bf/b))*(Snet/S) + ((pi/2)*(bf^2/S)); 
%CLa_A_h = 8.0; %(for verification)
%dClDalpha for tail-less aicraft 

zh = 0.600 + 0.4881 ;%vertical distance between wing and tail root chord taken from current geometry [m]
m_tv = 2*zh/b; % distance factor between horizontal tail and vortex shed plane of main wing [-]
r = 2*lh/b; % distance factor quarter chord main wing and tail [-]

%ADDITION FOR PROPELLER 
rho = 0.966632; % density at given altitude [kg/m^3]
Pbr= 132; %shaft horse power of one engine 132HP = 99000W (Rotax 915)
Cl= 0.645 ;%lift coefficient at given altitude for AoA = 0 with incidence -2 deg!
phi= asin(m_tv/r)*180/pi; %angle between r and m_tv

delta_s_de_da = 6.5*((rho*Pbr^2*S^3+Cl^3)/(lh^4*MTOW^3))^(0.25)*(sin(6*phi))^(2.5); %downwash propeller factor

ked = ((0.1121+0.1265*sweep_14+0.1766*sweep_14^2) / r^2 ) + 0.1024/r + 2;  %downwash corrective coefficient 
ked0 = 0.1124/r^2 + 0.1024/r + 2; %downwash corrective coefficient 
de_da = delta_s_de_da + (ked/ked0)*( (r/(r^2 + m_tv^2))*(0.4876/sqrt(r^2 + 0.6319 + m_tv^2))+...
    (1+(r^2/(r^2 + 0.7915+5.0734*m_tv^2))^0.3113)*(1-sqrt(m_tv^2/(1+m_tv^2))))...
    *(CLa_w/(pi*A)); %total downwash with added delta_s controbution for propeller

kn = -4.0;% for an engine positioned in front of the lemac/nose propeller
x_ac_w = 0.4; %aerodynamic center of wing (assumed at 0.4mac)
x_ac_c = x_ac_w - ((1.8/CLa_A_h)*(bf*hf*lfn/(S*mac))) +...
    ((0.273/(1+lambda))*((bf*cg*(b-bf))/(mac^2*(b+2.15*bf))))*tan(sweep_14) +...
    2*kn*((bn^2*ln)/(S*mac*CLa_A_h)); %total aircraft aerodynamic center

Vh_V = sqrt(0.85); %flow velocity ratio between H-tail and main wing [-]

%% controllability

CL_A_h = 1.04;% lift coefficient of wing+fuselage (without tail, landing configuration), from Data sheet given by Martin
dCl_max = 1.7; % change from zero to Clmax, from Data sheet given by Martin [-]
CL_h = -0.2;% for adjustable tail, from Data sheet given by Martin [-]
CL_0 = 0.8563; % from Data sheet given by Martin (for plane) [-]
Cm0 = -0.216; % % from Data sheet given by Martin (for main wing) [-]

         %Cmac_w = Cm0*((A*cos(sweep_14)^2)/(A+2*cos(sweep_14))); %pitching moment coefficient at aerodynamic center 
%for wing

         %Cmac_nac = 0.2; %assumed pitching moment coefficient at aerodynamic center (please estimate correctly when engine data available)
%for nacelle

         %Cmac_fus = -1.8*(1-(2.5*bf/lf))*((pi*bf*hf*lf)/(4*S*mac))*(CL_0/CLa_A_h); %pitching moment coefficient at aerodynamic center 
%for fuselage

         %Cm_ac = Cmac_w + Cmac_fus + Cmac_nac; %total moment coefficient at aerodynamic center
         
Cm_ac = -0.081; %From data sheet provided by Martin (ac assumed to be within +-10% of the neutral point)
%------------------------------------------------------------------------

x_cg_c = (-1:0.01:1);

% Stability Curve
Sh_S = ((x_cg_c) + SM - (x_ac_c))/((CLa_h/CLa_A_h)*(1-de_da)*(lh/mac)*Vh_V^2);

% Neutral Stability Curve (including stability margin)
Sh_S_NS = ((x_cg_c) - (x_ac_c))/((CLa_h/CLa_A_h)*(1-de_da)*(lh/mac)*Vh_V^2 );

% Controllablity Curve
Sh_S_C = ((Cm_ac/CL_A_h) + (x_cg_c)-(x_ac_c)) / ((CL_h/CL_A_h)*(lh/mac)*Vh_V^2);

%-------------------------------------------------------------------------

% plots
figure
plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C)
title('Scissors-plot: Stability & Controllability Curve')
xlabel('x_{cg}/c [%]')
ylabel('S_h/S [-]')
axis([-1 1 -0.5 0.6])
legend('Stability','Neutral Stability','Controllability')
