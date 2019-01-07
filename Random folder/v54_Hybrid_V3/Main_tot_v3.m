clc
clear all
close all
tic
AR = 10;
Lambda = - 2.8;
TR = 0.4;
tc = 0.18;
CLmax_Clean = 1.5;
CLmax_TO = 1.7;
CLmax_L = 2.1;


% Oswald factor (no propulsive interaction assumed)
e_clean=0.82;
e_TO=0.77;
e_L=0.73;
 
% zero-lift drag coefficient (no propulsive interaction assumed)
CD0_clean=0.033;
CD0_TO=0.068;
CD0_L=0.098;
save("CD0_matrix_1.mat","CD0_clean")
 

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



%% Main body

% Load input parameters
Input_new;

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
close Figure 29
close Figure 27
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
% vars = load('../ADSEE_I/variables_ADSEE_I.mat');





% P = double(vars.P); % FIND IT FROM ROELOF
P = WPdes.minGTM.GTM; % ADD THE SECONDARY POINT WHICH WE DO NOT KNOW for some reason
A = a.AR
MAC = double(summary_wing(11,2));           % <---- CHANGE FOR ELECTRIC/HYBRID
OEW = M.OEM + M.EM1 + M.EM2 + M.bat + M.f + M.GT + M.w;
S_ref = S;                        % <---- FROM ROELOF
v = V_cruise;        % <---- CHANGE FOR ELECTRIC/HYBRID
W4W5 = m.cr.f/m.L.f
% W_f = double(vars.W_fuel_total);  % <---- CHANGE FOR ELECTRIC/HYBRID
W_f = M.f
taper_ratio = double(summary_wing(3,2));    % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_c4 = double(summary_wing(4,2)); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_c2 = double(summary_wing(5,2)); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_LE = double(summary_wing(6,2)); % <---- CHANGE FOR ELECTRIC/HYBRID
sweep_TE = double(summary_wing(7,2)); % <---- CHANGE FOR ELECTRIC/HYBRID
b = double(summary_wing(2,2));               % <---- CHANGE FOR ELECTRIC/HYBRID
V_stall = 31.38;
LAMBDA = sweep_c4;    % Wingsweep at 0.25MAC
% Wfiml = 0.97 * MTOW * 9.81;      % Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48 !!!!! ELECTRIC = 1
Wfiml = frac * MTOW * 9.81;
WSsc = MTOW * m.cr.f * 9.81 / S_ref;      % Wing loading at the start of the cruise
WSec = W4W5 * WSsc;     % Wing loading at the end of the cruise

fus_length = 7.5; % <------ INPUT
empen_x = 0.9;  %   <---- INPUT WRONG WRONG WRONG
most_aft_cg = 0.6 * fus_length;     % <----- WRONG WRONG WRONG

%%%%%%%% added martin isa
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
M = v/a;



SF_S = 1.0;
[S_h, S_v] = control_surf_func(MAC, S, b, fus_length, empen_x, most_aft_cg);
S_h = SF_S*S_h;
S_v = SF_S*S_v;

S_ht = S_h;
S_vt = S_v;

A_h = 5.6;    % <----- INPUT   [3, 5] slide 68 lecture 7 ADSEE 1
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
L_t = 3.9*3.2808; % Tail length, wing quarter MAC to tail quarter MAC in ft   ??????
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





L1 = 1.7;
L2 = 2.8;
L3 = 3;

L = (L1+L2+L3)*3.281 ; %Fuselage structural length in ft for lecture 6 raymer pls dont hate


A_cs = 2.9;
D = sqrt(A_cs/pi); % derived from frontal area (even though fuselage may not be cilindrical)

mu = 1.7331332E-5; % viscosity of standard air at h=2400m (T=272K)












%% ADSEE II - Lecture 1

[Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, v, WSsc, WSec);
%% ADSEE II - Lecture 2

[CLmax, alpha_stall,CL_alpha_clean] = clalpha(A, clmax, CLdes, sweep_c2, M, sweep_LE, tc_avg);



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

cD = total_cD0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2
% cD = Fast_Cd0 + ADSEE_II_Drag.k_f(A, LAMBDA, CLdes) * (CLdes)^2


L_D = CLdes/cD

check_matrix = zeros(30, 1);    % MATRIX FOR SAVING CD0_OLIVERS
if total_cD0 - CD0_clean < 0.5
    disp("Difference is small")
end



