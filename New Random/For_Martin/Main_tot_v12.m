% TODO: THE POTATO PLOT IS SOMEHWAT WRONG, IT DOES NOT INCORPORATE THE NEW
% VALUE FOR XLEMAC, THUS IS A LITTLE BIT OFF. DO A SEPARATE FILE WITH
% RUBENS's CODE FOR THE CORRECT PLOT

clc
clear all
close all
tic
% AR = 10;
% Lambda = - 2.8;
% TR = 0.4;
% tc = 0.18;
% CLmax_Clean = 1.5;
% CLmax_TO = 1.7;
% CLmax_L = 2.1;
% 
% 
% % Oswald factor (no propulsive interaction assumed)
% e_clean=0.82;
% e_TO=0.77;
% e_L=0.73;
%  
% % zero-lift drag coefficient (no propulsive interaction assumed)
% CD0_clean=0.033;
% CD0_TO=0.068;
% CD0_L=0.098;
% save("CD0_matrix_1.mat","CD0_clean")
 

% Input_new(10)

disp("Done")





%% THe new added part
%%% Description
%
% This script combines the wing-loading power-loading diagrams and mission
% analysis to carry out the complete "Class-1.5" sizing of an HEP aircraft
% including aero-propulsive interaction effects. A description of the
% method can be found in the paper by de Vries, Brown and Vos (2018). The
% way in which aero-propulsive interaction effects are accounted for is
% briefly described in "WP_WS_diagram.m". The aero-propulsive models are
% specified in "WingPropDeltas.m". The powertrain model is described in
% "PowerTransmissionComputation_v2". The input parameters of the code are 
% described and specified in "Input.m". 
% 
% After running, the following variables (structures) should appear in the 
% workspace:
%
%   - a: contains aerodynamic properties such as the assumed CLmax or CD0
%       in the different flight conditions evaluated in the WS-WP diagram,
%       as well as wing geometries such as AR or TR. Defined by user in
%       "Input.m".
%   - AC: indicates which constraint is active/limiting for each powertrain
%       component's power loading, as well as wing loading, depending on
%       the design criteria chosen. For example, "AC.minGT.bat" gives the
%       name of the constraint which limits the minimum battery size when
%       selecting the design wing-loading corresponding to minimum gas
%       turbine size. AC is computed in "ComputeDesignPoint.m".
%   - AEROdes: provides the aerodynamic (such as CL, L/D,...) and
%       operational (such as v, M,...) properties of the aircraft in a
%       determined flight condition (= constraint) for a given design
%       point. For example, "AEROdes.minWS.cr.v" returns the velocity
%       required during the cruise constraint at the design wing loading
%       corresponding to minimum wing area (maximum W/S). AEROdes is
%       obtained in "ComputeDesignPoint.m".
%   - c: constants specified by the user in "Input.m", such as gravity or
%       sea-level atmospheric density.
%   - f: contains anonymous functions such as the power lapse of the
%       engine, atmospheric density lapse, or the weight of a powertrain
%       component as a function of installed power. Specified by user in
%       "Input.m".
%   - m: specifies operational/mission parameters, divided into flight
%       conditions (= constraints). For exmaple, "m.TO.h" is the altitude
%       at which take-off is performed. Some variables are specified by the
%       user in "Input.m"; others are computed along the way and added to
%       the structure.
%   - M: mass breakdown of the aircraft, as obtained from the mission
%       analysis. It is obtained from the "ComputeWeights.m" function.
%   - MA: collects results of the mission analysis ("MissionAnalysis.m",
%       unsuprisingly). One field exists per mission segment. For each
%       mission segment, several parameters are specified as an array, with
%       each element of the array corresponding to a timestep along the
%       mission segment. For example, "MA.Dcl.Ef" refers to the fuel energy
%       remaining on the aircraft (as a function of time) during the climb
%       phase of the diversion mission. MA also collects the resulting fuel
%       fractions and degree-of-hybridization of the aircraft.
%   - MA_in: input parameters for the mission analysis specified by the
%       user in "Input.m", such as Range, diversion altitude, and initial
%       guesses for OEM or DOH. Additional parameters are added per mission
%       segment throughout "MissionAnalysis.m".
%   - p: powertrain parameters specified by user in "Input.m", including
%       several which are flight-condition dependent. For example,
%       "p.L.etap1" refers to the propulsive efficiency of the primary
%       propulsors in landing conditions.
%   - s: program settings, such as convergence tolerances, figure numbers,
%       etc. Specified by user in "Input.m". The handles of the figures
%       generated are added to "s.figs".
%   - TW: Thrust-to-weight ratio of the different constraints evaluated in
%       the WP-WS diagram. Created in the different constraint functions
%       called in "WP_WS_diagram.m". The wing loading values at which the
%       TW arrays are specified are given by the structure "WS".
%   - TW_WSdes: Thrust-to-weight ratio evaluated at the different design
%       wing-loadings (specified in "WSdes") for each constraint. For
%       example, "TW_WSdes.minp.cr" refers to the thrust-to-weight ratio
%       obtained during cruise when selecting the wing loading
%       corresponding to minimum propulsive power as design point. TW_WSdes
%       is computed in "ComputeDesignPoint.m".
%   - WP_comp: contains the sizing power-loading values of each component 
%       of the powertrain, for each constraint and as a function of "WS".
%       For example, "WP_comp.cr.GB" gives, as a function of wing loading,
%       the sizing power-loading value (i.e. the maximum value the
%       component has to be able to produce/absorb) of the gearbox in
%       cruise conditions. This structure is generated as output of the
%       constraint functions in "WP_WS_diagram.m".
%   - WP_loss: contains the power-loading losses of each component 
%       of the powertrain, for each constraint and as a function of "WS".
%       For example, "WP_loss.cr.GB" gives, as a function of wing loading,
%       the amount of power lost due to e.g. heat dissipation (expressed
%       as a power-lodaing) in the gearbox in cruise conditions. This 
%       structure is generated as output of the constraint functions in 
%       "WP_WS_diagram.m".
%   - WP_path: contains the power-loading values of the different paths
%       that link the different components of the powertrain, for each 
%       constraint and as a function of "WS". For example, "WP_path.cr.s1" 
%       gives, as a function of wing loading, the power transmitted through
%       the primary shaft, from the gearbox to the primary propulsor 
%       (expressed as a power-loading). This structure is generated as 
%       output of the constraint functions in "WP_WS_diagram.m".
%   - WP_select: a series of WP-arrays obtained from WP_comp, WP_loss, or 
%       WP_path, which are manually selected at the end of 
%       "WP_WS_diagram.m" and passed on to "ComputeDesignPoint.m" in order
%       to evaluate the design points which lead to maximum power-loading
%       values of the fields specified (per constraint) in WP_select.
%   - WPdes: design power-loadings obtained from "ComputeDesignPoint.m" for
%       different design criteria. For example, "WPdes.minWS.GTM" refers to
%       the required power loading of the gas turbine (corrected to
%       TO/SL/max thrust conditions) when selecting the design point
%       corresponding to minimum wing area. The results in this structure
%       determine the installed power of the different powertrain 
%       components during the mission analysis, once the aircraft weight is 
%       known.
%   - WS: contains the wing-loading values at which TW and WP are sampled
%       for each constraint in the WP-WS diagram. These arrays are
%       generated in the constraint functions in "WP_WS_diagram.m".
%   - WSdes: wing-loading values of the design points selected. For
%       example, "WSdes.minbat" gives the wing loading value at which the
%       required battery power is minimized. Obtained from
%       "ComputeDesignPoint.m".
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 04-04-18
%%% Last modified: 16-04-18


%% Initialize
cd0_space = linspace(0.03, 0.05, 10);
save("cd0_matrix.mat", "cd0_space")


% %% Main body
% for index123 = [1:1:length(cd0_space)]
% save("index_number.mat", "index123");
% Load input parameters
Input_new_v2;

% Check input is consistent
disp([s.levelString '> Checking powertrain input settings'])
CheckInput;

% Generate wing-loading power-loading diagrams
disp([s.levelString '> Evaluating W/S-W/P diagram'])
WP_WS_diagram;

% Run mission analysis
disp([s.levelString '> Starting Mission Analysis'])
MissionAnalysis;


%% Generate additional plots if desired

% I. Powertrain plots
% (Note: the power paths used here come from the
% constraints, not the mission analysis! So the powers shown do not
% necessarily coincide with a given point on the MA mission profile)
if s.plotPowertrain == 1
    disp([s.levelString '> Generating powertrain plots'])
    P_path.cr = structfun(@(x) 1./x/m.cr.f*M.TOM*c.g/1e6,...
        WP_path.cr,'UniformOutput',0);
    P_path.TO = structfun(@(x) 1./x/m.TO.f*M.TOM*c.g/1e6,...
        WP_path.TO,'UniformOutput',0);
    P_path.L = structfun(@(x) 1./x/m.L.f*M.TOM*c.g/1e6,...
        WP_path.L,'UniformOutput',0);
    [s.figs(end+1)] = PlotPowertrain(P_path.cr,WS.cr,...
        WSdes.(s.SelDes),s.figStart+size(s.figs,2),'%6.2f',...
        'Cruise constraint [MW]');
    [s.figs(end+1)] = PlotPowertrain(P_path.TO,WS.TO,...
        WSdes.(s.SelDes),s.figStart+size(s.figs,2),'%6.2f',...
        'Take-off constraint [MW]');
    [s.figs(end+1)] = PlotPowertrain(P_path.L,NaN,...
        NaN,s.figStart+size(s.figs,2),'%6.2f',...
        'Landing constraint [MW]');
    clear('P_path')
end

% II. Aerodynamic polar in cruise conditions
if s.Polar.plot == 1
    disp([s.levelString '> Generating cruise polar'])
    [~,~,~,~,~] = CreateLiftDragPolarMap(WSdes.(s.SelDes),'cr',a,p,f,s,1);
end

% III. Power-control envelopes
if s.Env.plot == 1
    disp([s.levelString '> Generating power-control envelope(s)'])
    [s] = CreatePowerControlEnvelope(p,s,f,c,WPdes,AEROdes,MA,M);
end

% IV. Add drag requirements to stall constraint
if s.LandingCheck.plot == 1
    disp([s.levelString '> Adding landing drag requirements'])
    [s] = CreateLandingDragRequirements_v2(a,m,p,f,s,c,WPdes);
end


%% End
disp([s.levelString '> Completed. Run time: ' num2str(toc) ' seconds'])

figures_to_close = [18, 30, 29, 27, 25, 24, 23, 21, 15, 17, 14, 13, 12, 10]

close Figure 18
close Figure 30
% close Figure 29 % landing powertrain
% close Figure 27   % cruise constraints
close Figure 25
close Figure 24
close Figure 23
close Figure 21
close Figure 15
close Figure 17
close Figure 14
close Figure 13
close Figure 12
close Figure 10

save Figure

close all

V_cruise = m.cr.v;
h = m.cr.h;
S = Sw_roelof;
A = a.AR;
frac = m.cr.f;
m_cruise = MTOW * frac

wing_planform_design(V_cruise, A, S, m_cruise, h) % Done untill block 3
[summary_wing] = wing_planform_design(V_cruise, A, S, m_cruise, h); % m_cruise

summary_wing = [summary_wing; ["Wing Area", S]];















%% Drag part
P = MTOW*9.81/WPdes.minGTM.GTM; % ADD THE SECONDARY POINT WHICH WE DO NOT KNOW for some reason
A = a.AR;
MAC = double(summary_wing(11,2));           
OEW = M.OEM + M.EM1 + M.EM2 + M.bat + M.f + M.GT + M.w;
S_ref = S;                        
v = V_cruise;                     
W4W5 = m.cr.f/m.L.f;
W_f = M.f; % total fuel
taper_ratio = double(summary_wing(3,2));    
sweep_c4 = double(summary_wing(4,2));       
sweep_c2 = double(summary_wing(5,2));       
sweep_LE = double(summary_wing(6,2));       
sweep_TE = double(summary_wing(7,2));       
b = double(summary_wing(2,2));              
V_stall = 31.38;
LAMBDA = sweep_c4;    % Wingsweep at 0.25MAC
% Wfiml = 0.97 * MTOW * 9.81;      % Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48 !!!!! ELECTRIC = 1
Wfiml = frac * MTOW * 9.81;
WSsc = MTOW * m.cr.f * 9.81 / S_ref;      % Wing loading at the start of the cruise
WSec = W4W5 * WSsc;     % Wing loading at the end of the cruise

fus_length = 8.352;                       % <------ INPUT
empen_x = 0.95;                          %      <---- INPUT WRONG WRONG WRONG
most_aft_cg = 0.6 * fus_length;         % <----- WRONG WRONG WRONG INITIAL ROUGH CALCULATION

%%%%%%%% ISA CALCULATION
h = m.cr.h;
lambda_isa = 0.0065;
T_0_isa = 288.15;
g_isa = 9.81;
R_isa = 287.1;
P0_isa=101.325*10.^3;
T_isa = T_0_isa - lambda_isa * h;
P_isa = P0_isa * (T_isa / T_0_isa)^(g_isa / (lambda_isa * R_isa));
rho = P_isa / (R_isa * T_isa);
a = sqrt(1.4*R_isa*T_isa); % speed of sound
Mach = v/a;



SF_S = 1.0;
[S_h, S_v] = control_surf_func(MAC, S, b, fus_length, empen_x, most_aft_cg);
S_h = SF_S*S_h;
S_v = SF_S*S_v;

S_ht = S_h;
S_vt = S_v;

A_h = 6.0;    % <----- INPUT   [3, 5] slide 68 lecture 7 ADSEE 1 UsE THIS
A_v = 1.84;  % <----- INPUT   [1, 2] slide 68 lecture 7 ADSEE 1

b_h = sqrt(A_h * S_h);
b_v = sqrt(A_v * S_v);



%inputs (page 477-479 Raymer) (everything is in retard units) (lbs,
%gallons, ft^3, ft^2, inch etc.)
W_dg = MTOW * 0.9 * 2.2; % Design gross weight
N_z = 4.4; % Load factor
N_gear = 3; % Find Raymer!!!!
lambda = taper_ratio; % taper ratio

LAMBDA_ht = sweep_c4; % Sweep at 25% MAC
A_ht = A_h; % Aspect ratio horizontal tailwing                   ??????
H_t_over_H_v = 0; % = 0 for conventional tail, 1 for T-tail
LAMBDA_vt = LAMBDA_ht; % Sweep at 25% of vertical tail MAC                    ??????
A_vt = A_v; % Aspect ratio vertical tail                         ??????


lambda_vt = 1.; % taper raio vertical tail                                    ??????
lambda_h = 1; %Taper ratio horizontal tail                                    ??????
L_t = 4.1176*3.2808; % Tail length, wing quarter MAC to tail quarter MAC in ft   ??????
% N_en = double(vars.N); %nr of engines 
N_en = 2;

W_press = 0 ;%11.9+(V_pr*P_delta)^0.271; %Weight penalty due to pressurization; PROBABLY ZERO FOR OUR DESIGNS BECAUSE WE DON'T PRESSURIZE OUR CABIN
W_l = (MTOW - W_f) * 2.2; %Landing design gross weight

V_t = W_f / (0.840 * 3.79); %Total fuel volume in gallons

L_m = 15.; %Extended length of main landing gear                ??????
L_n = 15.; %Extended nose gear length (inch)                    ??????


%% Engine Weight Inputs
W_en = P/1340 * 2.2; %Total engine weight (each) in pounds              <---- INPUT

%%
V_i = V_t * 1.05; %Integral tanks volume in gallons

N_t = 1; %Number of fuel tanks                                  ??????
W_uav = 0.03 * MTOW * 2.2; %Uninstalled avionics weight in pounds
N_p = 5; %Number of personal onboard

% cl = 0.30647;
cl_cruise  = double(summary_wing(8,2));
clmax = CLmax_L;               % XXXX

cambered = 0; % 1 for True, 0 for False
e = e_clean;
c_avg = sqrt(S_ref/A);  %                           ?????? WHICH CHORD IS THIS????


MAC_ht = b_h/A_ht;
MAC_vt = b_v/A_vt;

tc_avg = double(summary_wing(14, 2)); % (t/c)_avg is the average thickness to chord
xc_max = 0.25; % (x/c)_max is the position of maximum thickness         ????????

% C_f_e = 0.0055; % light AC - single engine
% C_f_e = 0.0045; % light AC - twin engine

% k = 0.152E-5; % polished sheet metal
% k = 0.634E-5; % smooth paint
k = 0.052E-5;   % composites





L1 = 1.7;               %       <--------- INPUT
L2 = 2.534;               %       <--------- INPUT
L2 = 2.8;
L3 = 3.5;                 %       <--------- INPUT

L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate


A_cs = 2.34; % cross sectional area of the fuselage .  <-------- INPUT
D = sqrt(A_cs/pi); % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)



%% ADSEE II - Lecture 1

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, v, WSsc, WSec);
%% ADSEE II - Lecture 2

[CLmax, alpha_stall,CL_alpha_clean] = clalpha(A, clmax, CLdes, sweep_c2, Mach, sweep_LE, tc_avg);



%% ADSEE II - Lecture 3 - Drag coefficient estimation
%% Component contributions propeller AC (fast method)
Component = ['Wing', 'Fuselage multi-engine', 'Fuselage single-engine', 'Nacelles', 'Tail (hor + ver)', 'misc'];
C_D_Cs = [0.007, 0.08, 0.11, 0.06, 0.008, 0.15];
% Change these according to component name (defined above)
A_Cs = [15, 0, 3.0, 0, 0.3];
Fast_Cd0 = ADSEE_II_Drag.fast_sum_C_D_0(C_D_Cs, A_Cs, S_ref);

S_w = ADSEE_II_Drag.S_wet_c(S_ref, S_ht, S_vt, D, L1, L2, L3);

%% Component method
% [Fuselage, Wing, horizontal tail, vertical tail]
C_f_c_fuselage = ADSEE_II_Drag.fp_skin_friction(0.25, k, rho, v, L2, mu, a); % increase
C_f_c_wingtail = ADSEE_II_Drag.fp_skin_friction(0.5, k, rho, v, L3, mu, a);
C_f_cs = [C_f_c_fuselage, C_f_c_wingtail, C_f_c_wingtail, C_f_c_wingtail];
option = [2, 1, 1, 1]; % length? it was 2 1 1 1
FF_cs = [ADSEE_II_Drag.form_factor(option(1), L2, D, tc_avg, xc_max, LAMBDA, v, a),...
    ADSEE_II_Drag.form_factor(option(2), L2, D, tc_avg, xc_max, LAMBDA, v, a),...
    ADSEE_II_Drag.form_factor(option(3), L2, D, tc_avg, xc_max, LAMBDA, v, a),...
    ADSEE_II_Drag.form_factor(option(4), L2, D, tc_avg, xc_max, LAMBDA, v, a)];
IF_cs = [1.0, 1.25, 1.05, 1.05];
S_cs = ADSEE_II_Drag.S_wet_c(S_ref, S_ht, S_vt, D, L1, L2, L3);
% ---------------------------------------------------------->
%                                                           |
cd0_c = ADSEE_II_Drag.tot_comp_drag0(C_f_cs, FF_cs, IF_cs, S_cs, S_ref, 0);
% misc = ADSEE_II_Drag.cD_misc0(0.034, A_cs, L2*D*0.1, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c, c);
misc = ADSEE_II_Drag.cD_misc0(0.034, A_cs, L2*D*0, v, a, 0.6, 1., 0, 0.1*S_ref, S_ref, 0.1*c_avg, c_avg);


total_cD0 = cd0_c + misc;
save("CD0_matrix_2.mat","total_cD0")

cD = total_cD0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2;
% cD = Fast_Cd0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2

cD
L_D = CLdes/cD
total_cD0
CD0_clean


%% Payload range
payload_new = 363;
percent_emptiness_payload = 0.5;
PR_func(MTOW, OEW, M.f, payload_new, L_D, W4W5, percent_emptiness_payload)


%% Aileron calculation
c_r = double(summary_wing(9,2));
c_t = double(summary_wing(10,2));
c_l_alpha = 0.32;
[aileron_length] = AileronNEW(c_r, c_t, sweep_LE, sweep_TE, c_l_alpha,...
    S_ref, total_cD0, V_stall, b)








%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART
%%%% BROKEN UP PART%%%% BROKEN UP PART%%%% BROKEN UP PART
%%%% BROKEN UP PART

















%% Scissor and Potato plots
thick = 2;
Fuel = M.f;
Batteries = M.bat*9.81;
%Dimensions and parameters (fixed) 
Snet=S*2; %wetted area
cg = S/b; %average constant chord [m]
Mach_1 = 0.27; %mach number at cruise [-]
beta=sqrt(1-Mach_1^2);
bf = 1.416; %fuselage width [m]                              % <----- INPUT           
hf = 1.65; %fuselage height [m]                            % <----- INPUT                      
lf = L1 + L2 + L3; % total length of fuselage [m]       
ln = 0.25*0; %distance from engine to quater chord mac [m] . % <----- INPUT   
lfn = 2.85; %nose to leading edge [m] (GUESS) .            % <----- INPUT   
bn = 0.5; %width of nacelles (engines) [m] .               % <----- INPUT   
Ct = c_t;%tip chord [m]
Cr = c_r;%root chord [m]
lambda = Ct/Cr; %taper ratio of main wing [-]
Ah = 6 - 0; %aspect ratio of horizontal tail [-] .           % <----- INPUT   
bh= 2.82; %horizontal tail span [m] .                      % <----- INPUT   
Sh = bh^2/Ah; %Area of the horizontal tail wing [m] .      % <----- INPUT   NOT NEEDED I THINK ? MAYBE NO SWEEP AND LEAVE IT OUT
cg_h= Sh/bh; %average chord
lambda_h = 0.39 %taper ratio of tail wing [-] (
Cr_h = 2*Sh/((1+lambda_h)*bh);  %     NOT NEEDED
Ct_h = lambda_h*Cr_h;           %  INPUT   NOT NEEDED
eta = 0.95; %airfoil efficiency factor [-]


% Calculate wing sweep angle
sweep_LE = 2.9*pi/180;  %sweep at leading edge [rad]  
sweep_4 = atan(tan(sweep_LE) + (Cr/(2*b))*(lambda -1)); %sweep at quater chord [rad]
sweep_2 = -2.9*pi/180;%atan(tan(sweep_LE) - (4/A)*(0.5*((1-lambda)/(1+lambda)))); %sweep at half chord [rad]

% Calculate tail wing sweep angle (could be ifferent)
sweep_LE_h = 0*pi/180;
sweep_4_h = atan(tan(sweep_LE_h) + (Cr_h/(2*bh))*(lambda_h -1));
sweep_2_h = -2.9*pi/180;

% measured from planform for given geometry (from nose to horizontal tail) [m]
x_datum_h = 0.9 * fus_length ;   %(ASSUMED)                  % <----- INPUT               

%distance between aerodynamic center of main wing and horizontal tail [m] 
% lh = L_t/3.2808; %(ASSUMED) .                          % <----- INPUT   
lh = L3 +L2/2;
lh = 4.0 + 0.5;
 

% % % % % x locations
% % % % x_lemac = 3.55;          % x location of leading edge mean aerodynamic chord [m] (GUESS)
% % % % x_OEW = x_lemac + 0.375*MAC;        % assumed CG of operational empty weight [m]
% % % % x_Cargo = 4.9;                       % assumed CG of cargo in meters [m]
% % % % x_Fuel = x_lemac + 0.5*MAC;         % CG of fuel [m]

%% Variables Stability 

SM = 0.05; 
%stability margin for safety given as percentage/100
CLaw = (2*pi*A)/(2 + sqrt(4 + (A/eta)^2 *(1 + (tan(sweep_2)^2/beta^2)) )) ; 
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for main wing
CLah = (2*pi*Ah)/(2+ sqrt(4 + (Ah/eta)^2 *(1 + (tan(sweep_2_h)^2/beta^2)) ));
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for horizontal tail
CLaAh = CLaw*(1+(2.15*bf/b))*(Snet/S) + ((pi/2)*(bf^2/S)); 
%dClDalpha for tail-less aicraft 

%% Downwash

zh = 0.825;%vertical distance between wing and tail root chord taken from current geometry [m] (ASSUMED) % <----- INPUT   
m_tv = 2*zh/b; % distance factor between horizontal tail and vortex shed plane of main wing [-]
r = 2*lh/b; % distance factor quarter chord main wing and tail [-]

%For a propeller on the wing 
% rho= 0.9; % density at given altitude [kg/m^3]
Pbr= 40; %shaft horse power of one engine 132HP = 99000W (Rotax 915) (ASSUMED) % <----- INPUT INPUT
Cl= double(summary_wing(8,2));%lift coefficient at given altitude  (cl_cruize)
phi= asin(m_tv/r)*180/pi; %angle between r and m_tv
%assume 0, as the engine is located on the tip of the fuselage for the fuel
%configuration
delta_s_de_da = 6.5*((rho*Pbr^2*S^3+Cl^3)/(lh^4*MTOW^3))^(0.25)*(sin(6*phi))^(2.5); %downwash propeller factor

ked = ((0.1121+0.1265*sweep_4+0.1766*sweep_4^2) / r^2 ) + 0.1024/r + 2;  %downwash corrective coefficient 
ked0 = 0.1124/r^2 + 0.1024/r + 2; %downwash corrective coefficient 
de_da = delta_s_de_da + (ked/ked0)*( (r/(r^2 + m_tv^2))*(0.4876/sqrt(r^2 + 0.6319 + m_tv^2))+...
    (1+(r^2/(r^2 + 0.7915+5.0734*m_tv^2))^0.3113)*(1-sqrt(m_tv^2/(1+m_tv^2))))...
    *(CLaw/(pi*A)); %total downwash with added delta_s controbution for propeller !!!!!! delta_s_de_da 000

%% Aerodynamic center
kn = -4;% for an engine positioned in front of the lemac/nose propeller
x_ac_w = 0.25; %aerodynamic center of wing (assumed at 0.4mac)

x_ac_c = x_ac_w - ((1.8/CLaAh)*(bf*hf*lfn/(S*MAC))) +...
    ((0.273/(1+lambda))*((bf*cg*(b-bf))/(MAC^2*(b+2.15*bf))))*tan(sweep_4) +...
    2*kn*((bn^2*ln)/(S*MAC*CLaAh)); %total aircraft aerodynamic center


%% Speed on the tail an wing
Vh_V = sqrt(0.85); %flow velocity ratio between H-tail and main wing [-]

%% Controllability
CLAh = 1.25;% lift coefficient of wing+fuselage (without tail, landing configuration) . % <----- INPUT   1.2
% lecture 4 slide 37
CLh = -0.35*(Ah)^(1/3);% for fixed .                   % <----- INPUT   
CL0 = 0.8563 -0.3; % Zero incident while flaps out . % <----- WRONG WRONG
Cm0 = -0.216; %  (for main wing) [-]
CL=CLmax_clean;

mu_1=0.17; % accounts for the effect of the camber increase associated with the deflection angle of the flaps
mu_2=0.5;%represent the 2D effect lift contribution generated by the high lift devices
mu_3=0.55;% accounts for the effect of sweep angle
cf_c= 1.0; %extended chord due to flaps and normal chord length
S_wf_S=1.0; %ratio between flapped wing area and reference wing area 
CLmax= CLmax_L; %aircraft lift at landing
dClmax = CLmax - CL;
%% Cmac zero-lift pitching moment
Cmac_w = Cm0*((A*cos(sweep_4)^2)/(A+2*cos(sweep_4))); %pitching moment coefficient at aerodynamic center for wing

Cmac_nac = 0.2; %assumed pitching moment coefficient at aerodynamic center (please estimate correctly when engine data available)
%for nacelle

Cmac_fus = -1.8*(1-(2.5*bf/lf))*((pi*bf*hf*lf)/(4*S*MAC))*(CL0/CLaAh); %pitching moment coefficient at aerodynamic center 
% for fuselage

Delta_f_Cmac = mu_2*(-mu_1*dClmax*cf_c-[CL+dClmax*(1-S_wf_S)]*cf_c/8*(cf_c-1))+0.7*A/(1+2/A)*mu_3*dClmax*tan(2.8)% flaps

Cm_ac = Cmac_w + Cmac_fus + Cmac_nac + Delta_f_Cmac %total moment coefficient at aerodynamic center
Cm_ac = -0.62;
         
%Cm_ac = -0.65; %(ac assumed to be within +-10% of the neutral point)

%% Equations
x_cg=(-1:0.01:1);

% Stability curve
Sh_S= (x_cg-x_ac_c-SM-0.01)/((CLah/CLaAh)*(1-de_da)*((Vh_V)^(2))*(lh/MAC)); %Stability gradient

% Controllablity Curve
Sh_S_C = ((Cm_ac/CLAh) -(x_ac_c)) / ((CLh/CLAh)*(lh/MAC)*Vh_V^2)+ (x_cg)/ ((CLh/CLAh)*(lh/MAC)*Vh_V^2);

%% Rotation Take off setting
% xG = 4.6;%new obtaine position of the undercarriage [m]
% 
% VS1=31.6; % vstall speed (requirement);
% Vh = 31; %random assumption
% VR =33.18; % vstall * 1.05 (ASSUMPTION)
% nh = (x_datum_h-xG)/lh*(Vh/VR)^2; 
% theta = 0.5;
% nq= 1+ CLah*theta*(x_datum_h-xG)/(CLh*VR);
% 
% Sh_S_R=(CLmax/(nh*nq*CLh))*((Cm_ac/CLmax)-(VS1/VR)^2*((xG-x_cg)/MAC))+(CLAh/CLh*(xG/MAC-0.25));


fscissor = figure
yyaxis right
plot(x_cg,Sh_S,x_cg,Sh_S_C, "LineWidth", thick)
xlabel('x_{cg}/MAC [%]')
ylabel('S_h/S [-]')
axis([0.1 0.5 0 0.4])


%% POTATO PLOT
% ADD ROUND UP TO THIRD DECIMAL POINT AND EXTRACT THE VALUE OF THE CHOSEN
% LEMAC 
% clear all
x_lemac = [2: 0.01: 4.2];
cg_mat = zeros(length(x_lemac),2);
counter = 1
for i  = x_lemac
    lbs_to_kg = 0.45359237;
    mass_pax=175;                    %lbs               
    mass_pax = mass_pax*lbs_to_kg;

    mass_bags = 25;                                
    mass_bags = mass_bags * lbs_to_kg;

    mass_fuel = M.f;

    W_OEW = OEW;                                     
    cg_OEW = 2.865-0.2;     % 2.665                    % <--------- INPUT INPUT INPUT INPUT INPUT INPUT INPUT INPUT
    cg_OEW = 3.25;

    seat_pilot=1.595 ;                    %c.g. Position Pilot   <----- INPUT 
    seat_row1=2.721  ;                    %c.g. Position Row 1   <----- INPUT 
    seat_row2=4.063  ;                     %c.g. Position Row 2   <----- INPUT 
    location_cargo = 1.923;                %c.g. Position Baggage <----- INPUT  1.923
    location_fuel = 3.762;              %c.g. Position Baggage <----- INPUT
    location_fuel = i + 0.5*MAC
    %location_batteries%c.g. Position Fuel    <----- INPUT
%     location_fuel = cg_OEW + 0.10*l_fus;
%     location_feul = x_Fuel;

    %cargo
    W_OEW_cargo = W_OEW+4*mass_bags;
    cg_OEW_cargo=((cg_OEW*W_OEW)+(location_cargo*4*mass_bags))/(W_OEW_cargo);

   
    W_OEW_1pax=1*mass_pax+W_OEW_cargo;
    W_OEW_2pax=2*mass_pax+W_OEW_cargo;
    W_OEW_3pax=3*mass_pax+W_OEW_cargo;
    W_OEW_4pax=4*mass_pax+W_OEW_cargo;
    %back to front
    cg_btf_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row2*mass_pax))/(W_OEW_1pax);
    cg_btf_2=((cg_btf_1*W_OEW_1pax)+(seat_row2*mass_pax))/(W_OEW_2pax);
    cg_btf_3=((cg_btf_2*W_OEW_2pax)+(seat_row1*mass_pax))/(W_OEW_3pax);
    cg_btf_4=((cg_btf_3*W_OEW_3pax)+(seat_row1*mass_pax))/(W_OEW_4pax);


    %front to back
    cg_ftb_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row1*mass_pax))/(W_OEW_1pax);
    cg_ftb_2=((cg_ftb_1*W_OEW_1pax)+(seat_row1*mass_pax))/(W_OEW_2pax);
    cg_ftb_3=((cg_ftb_2*W_OEW_2pax)+(seat_row2*mass_pax))/(W_OEW_3pax);
    cg_ftb_4=((cg_ftb_3*W_OEW_3pax)+(seat_row2*mass_pax))/(W_OEW_4pax);


    %include fuel
    cg_nofuel=cg_btf_4;
    cg_fuel=((cg_nofuel*W_OEW_4pax)+(location_fuel*mass_fuel))/(mass_fuel+W_OEW_4pax);


    cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);
    cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);

    cg_max = (cg_max - i)/MAC;
    cg_min = (cg_min - i)/MAC;
    cg_mat(counter, 1) = round(cg_min,3);
    cg_mat(counter, 2) = round(cg_max,3);
    
    counter = counter + 1;
    
end

%figure
yyaxis left
plot([cg_mat(:,1), cg_mat(:,2)], x_lemac/lf, "LineWidth", thick)
ylim([0.25 0.34])
xlim([0.1 0.5])
hold on
% plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C)
xlabel("x_{cg}/MAC", "FontSize", 30)
ylabel("x_{LEMAC}/L_{FUS}", "FontSize", 30)

legend("FORWARD CG", "AFT CG", "Stability", "Controlability")
set(gca,'FontSize',25);

% close all





prompt_xlemac = 'X_lemac position: ';
x_lemac_scissor = double(input(prompt_xlemac));


prompt_aftcg = 'Most AFT cg: ';
most_aft_cg = double(input(prompt_aftcg));

prompt_forwardcg = 'Most FORWARD cg: ';
most_forward_cg = double(input(prompt_forwardcg));

prompt_ShS = 'Chosen ShS: ';
ShS_ratio = double(input(prompt_ShS));

Sh_final = S*ShS_ratio;
A_h = 6.0;    % <----- INPUT   [3, 5] slide 68 lecture 7 ADSEE 1
A_v = 1.5;  % <----- INPUT   [1, 2] slide 68 lecture 7 ADSEE 1

b_h = sqrt(A_h * S_h);
b_v = sqrt(A_v * S_v);


% wing_planform_design(V_cruise, A_ht, Sh_final, m_cruise, h) % Done untill block 3
[summary_horiz_tail] = wing_planform_design(V_cruise, A_ht, Sh_final, m_cruise, h); % m_cruise
Cr_h_final = double(summary_horiz_tail(9,2));
Ct_h_final = double(summary_horiz_tail(10,2));
MAC_h_final = double(summary_horiz_tail(11,2));
b_h_final = b_h;



ind = find(cg_mat(:,2) == most_aft_cg);
% get the index and the run the loop untill this index with the ne x_lemac



%% POTATO PLOT VERSION two
% ADD ROUND UP TO THIRD DECIMAL POINT AND EXTRACT THE VALUE OF THE CHOSEN
% LEMAC 
% clear all
% x_lemac = [1: 0.01: 4.2];
cg_mat = zeros(length(x_lemac),2);
counter = 1;
x_lemac = x_lemac_scissor*lf

for i  = x_lemac
%     i = x_lemac_scissor*lf; % PROBABLY NOT CORRECT
    
    lbs_to_kg = 0.45359237;
    mass_pax=175;                    
    mass_pax = mass_pax*lbs_to_kg;

    mass_bags = 25;                                
    mass_bags = mass_bags * lbs_to_kg;

    mass_fuel = M.f;

    cg_OEW = 3.25;

    seat_pilot=1.595 ;                    %c.g. Position Pilot   <----- INPUT 
    seat_row1=2.721  ;                    %c.g. Position Row 1   <----- INPUT 
    seat_row2=4.063  ;                     %c.g. Position Row 2   <----- INPUT 
    location_cargo = 1.923;                %c.g. Position Baggage <----- INPUT  1.923
    location_fuel = 3.762;              %c.g. Position Baggage <----- INPUT
    location_fuel = i + 0.5*MAC
    %location_batteries%c.g. Position Fuel    <----- INPUT
%     location_fuel = cg_OEW + 0.10*l_fus;
%     location_feul = x_Fuel;

    %cargo
    W_OEW_cargo = W_OEW+4*mass_bags;
    cg_OEW_cargo=((cg_OEW*W_OEW)+(location_cargo*4*mass_bags))/(W_OEW_cargo);

   
    W_OEW_1pax=1*mass_pax+W_OEW_cargo;
    W_OEW_2pax=2*mass_pax+W_OEW_cargo;
    W_OEW_3pax=3*mass_pax+W_OEW_cargo;
    W_OEW_4pax=4*mass_pax+W_OEW_cargo;
    %back to front
    cg_btf_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row2*mass_pax))/(W_OEW_1pax);
    cg_btf_2=((cg_btf_1*W_OEW_1pax)+(seat_row2*mass_pax))/(W_OEW_2pax);
    cg_btf_3=((cg_btf_2*W_OEW_2pax)+(seat_row1*mass_pax))/(W_OEW_3pax);
    cg_btf_4=((cg_btf_3*W_OEW_3pax)+(seat_row1*mass_pax))/(W_OEW_4pax);


    %front to back
    cg_ftb_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row1*mass_pax))/(W_OEW_1pax);
    cg_ftb_2=((cg_ftb_1*W_OEW_1pax)+(seat_row1*mass_pax))/(W_OEW_2pax);
    cg_ftb_3=((cg_ftb_2*W_OEW_2pax)+(seat_row2*mass_pax))/(W_OEW_3pax);
    cg_ftb_4=((cg_ftb_3*W_OEW_3pax)+(seat_row2*mass_pax))/(W_OEW_4pax);


    %include fuel
    cg_nofuel=cg_btf_4;
    cg_fuel=((cg_nofuel*W_OEW_4pax)+(location_fuel*mass_fuel))/(mass_fuel+W_OEW_4pax);


    cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);
    cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);

    cg_max = (cg_max - i)/MAC;
    cg_min = (cg_min - i)/MAC;
    cg_mat(counter, 1) = round(cg_min,3);
    cg_mat(counter, 2) = round(cg_max,3);
    
    counter = counter + 1;
%     if counter-1 == ind
    if i == x_lemac
        disp(i)
        disp(i)
        disp(i)
        disp(i)
        disp(i)
        fpotato = figure
%         x_lemac = x_lemac * lf; %%% AFTER x_lemac was chosen
        line([([cg_OEW,cg_OEW_cargo]- i)/MAC],[[W_OEW,W_OEW_cargo]],'Color','green');
        line([([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4]-i)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');
        line([([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4]-i)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');
        line([([cg_nofuel,cg_fuel]-i)/MAC],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');

        xlabel("x_{cg}/MAC")
        ylabel("Mass [kg]")
        legend("Cargo", "Front to back", "Back to front", "Fuel")
        set(gca, "FontSize", 20)
        saveas(fpotato,'PotatoDiagram.fig');
        close all
        
        break
    end
end

%% LANDING GEAR POSITION
perc_mac = 0.475; % 40% of the mac is the landing gear


[length_strut, l_gear_n, z, Ymlg ] = lg_position(perc_mac, lf, L3, ...
     x_lemac_scissor*lf + most_aft_cg*MAC, x_lemac_scissor*lf, MAC, bf, MTOW) 

%% DCLMAX 
[dCLMAX, TVOL] = planformlayout(b, c_r, c_t, bf, aileron_length, payload_new, M.f, S);
dCLMAX



% 