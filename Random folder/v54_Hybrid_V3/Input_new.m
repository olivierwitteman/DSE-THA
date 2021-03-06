% function InputFunc = input_new(AR, Lambda, TR, tc,...
%    CLmax_clean, CLmax_TO, CLmax_L,...
%    e_clean, e_TO, e_L,...
%    CD0_clean, CD0_TO, CD0_L)
%% Aerodynamic / wing properties
ind  = load("index_number.mat");
cd0s = load("cd0_matrix.mat");



% General
a.AR = 10;                                      % Wing aspect ratio [-]
a.Lambda = -2.8;                                   % Half-chord sweep angle of wing [deg]
a.TR = 0.4;                                     % Taper ratio of wing [-]
a.tc = 0.18;                                    % Thickness-to-chord ratio of root section [-]
a.nUlt = 4.4;                                   % Ultimate load factor [-]
 
% maximum lift coefficient (no propulsive interaction assumed)
CLmax_clean=1.55;
CLmax_TO=1.7;
CLmax_L=2.1;
 
% Oswald factor (no propulsive interaction assumed)
e_clean=0.82;
e_TO=0.77;
e_L=0.73;
 
% zero-lift drag coefficient (no propulsive interaction assumed)
CD0_clean = 0.0357;
CD0_TO= CD0_clean + 0.015;
CD0_L= CD0_clean + 0.035;
 
% Cruise
a.cr.CD0 = CD0_clean;                         % Cruise zero-lift drag coefficient [-]
a.cr.e = e_clean;                             % Cruise oswald factor [-]
 
% Landing
a.L.CD0 = CD0_L;                              % Landing zero-lift drag coefficient [-]
a.L.e = e_L;                                  % Landing oswald factor [-]
a.L.CLmax = CLmax_L;                          % Landing maximum lift coefficient of isolated wing [-]
 
% Take off
a.TO.CD0 = CD0_TO;                            % TO zero-lift drag coefficient [-]
a.TO.e = e_TO;                                % TO oswald factor [-]
a.TO.CLmax = CLmax_TO;                        % TO maximum lift coefficient of isolated wing [-]
 
% OEI Balked landing
a.bL.CD0 = CD0_L;                           % OEI balked landing zero-lift drag coefficient (LG retracted) [-]
                                              % CLmax, e, assumed to be the same as in landing configuration.
% OEI Ceiling
a.cI.CD0 = CD0_clean;                         % Ceiling (OEI) zero-lift drag coefficient [-]
a.cI.e = e_clean ;                            % Ceiling (OEI) oswald factor [-]
a.cI.CLmax = CLmax_clean;                     % Clean configuration max lift coefficient of isolated wing [-]
 
% Start-of-climb
a.cl.CD0 = CD0_clean;                         % Start-of-climb zero-lift drag coefficient [-]
a.cl.e = e_clean;                             % Start-of-climb oswald factor [-]
a.cl.CLmax = CLmax_clean;                     % Start-of-climb maximum lift coefficient of isolated wing [-]
 
% Top-of-climb
a.ct.CD0 = CD0_clean;                         % Top-of-climb zero-lift drag coefficient [-]
a.ct.e = e_clean;                             % Top-of-climb oswald factor [-]
a.ct.CLmax = CLmax_clean;                         % Top-of-climb maximum lift coefficient of isolated wing [-]
 
 
%% Propulsion System
 
% Propulsion system layout
p.config = 'SPPH';                          % Powertrain architecture ('conventional', 'turboelectric', 'serial',
                                                %   'parallel', 'PTE', 'SPPH', 'e-1', 'e-2', or 'dual-e')
p.b_dp = 0.5;                                   % Fraction of wing span occupied by DP system [-]
p.dy = 0.01;                                       % Spacing between adjacent DP propulsors, as fraction of propulsor diameter [-]
p.N1 = 2;                                       % Number of chains in primary powertrain [-]
p.N2 = 6;                                       % Number of chains in secondary powertrain[-]
p.DP = 2;                                       % Which PS has an effect on wing performance? (1 = primary, 2 = secondary, 0 = none)
p.xp = -0.25;                                     % Axial position of propellers as a fraction of chord
                                                % xp < 0: tractor
                                                % 0 < xp < 1: OTW
                                                % xp > 1: pusher
                                                % xp = Inf: No effect of prop on wing
% Component properties (excl. propulsive)
p.eta_EM1 = 0.95;                               % Conversion efficiency of (electro-) generators
p.eta_EM2 = 0.95;                               % Conversion efficiency of electromotors
p.eta_PM = 0.95;                                % Conversion efficiency of PMAD
p.eta_GB = 0.95;                                % Transmission efficiency of gearboxes
p.eta_GT = 0.35;                                % Conversion (thermal) efficiency of gas turbine
p.SE.bat = 1.8e6;                               % Battery specific energy [J/kg] 
p.SE.f = 42.8e6;                                % Fuel specific energy [J/kg]
p.SP.EM = 7.7e3;                                % Electrical machine specific power [W/kg]
p.SP.bat = 2000;                                % Battery pack specific power [W/kg]
p.minSOC_miss = 0.2;                            % Minimum SOC (maximum discharge) of batteries after nominal mission [-]
p.minSOC_tot = 0;                               % Minimum SOC (maximum discharge) of batteries after diversion mission [-]
 
% Cruise
p.cr.etap1 = 0.8;                               % Primary propulsors' propulsive efficiency in cruise (of ISOLATED propulsors) [-]
p.cr.etap2 = 0.8;                              % Secondary propulsors' propulsive efficiency in cruise (of ISOLATED propulsors) [-]
p.cr.Gamma = 0;                                 % Thrust vectoring in cruise [deg]
 
% Landing
p.L.etap1 = 0.65;                                % Primary propulsors' propulsive efficiency in landing conditions (of ISOLATED propulsors) [-]
p.L.etap2 = 0.75;                                % Secondary propulsors' propulsive efficiency in landing conditions (of ISOLATED propulsors) [-]
p.L.Gamma = 0;                                  % Thrust vectoring in landing configuration [deg]
 
% Take off
p.TO.etap1 = 0.80;                              % Primary propulsors' propulsive efficiency in TO conditions (of ISOLATED propulsors) [-]
p.TO.etap2 = 0.70;                              % Secondary propulsors' propulsive efficiency in TO conditions (of ISOLATED propulsors) [-]
p.TO.Gamma = 0;                                 % Thrust vectoring in TO configuration [deg]
 
% OEI Balked landing
p.bL.etap1 = 0.8;                               % Primary propulsors' propulsive efficiency in balked landing conditions (of ISOLATED propulsors) [-]
p.bL.etap2 = 0.7;                               % Secondary propulsors' propulsive efficiency in balked landing conditions (of ISOLATED propulsors) [-]
p.bL.Gamma = 0;                                 % Thrust vectoring in balked landing configuration [deg]
 
% OEI ceiling
p.cI.etap1 = 0.8;                               % Primary propulsors' propulsive efficiency in ceiling (OEI) conditions (of ISOLATED propulsors) [-]
p.cI.etap2 = 0.8;                               % Secondary propulsors' propulsive efficiency in ceiling (OEI) conditions (of ISOLATED propulsors) [-]
p.cI.Gamma = 0;                                 % Thrust vectoring in ceiling (OEI) configuration [deg]
 
% Start-of-climb
p.cl.etap1 = 0.75;                              % Primary propulsors' propulsive efficiency in climb conditions (of ISOLATED propulsors) [-]
p.cl.etap2 = 0.75;                              % Secondary propulsors' propulsive efficiency in start-of-climb conditions (of ISOLATED propulsors) [-]
p.cl.Gamma = 0;                                 % Thrust vectoring in start-of-climb configuration [deg]
 
% Top-of-climb
p.ct.etap1 = 0.8;                               % Primary propulsors' propulsive efficiency in top-of-climb conditions (of ISOLATED propulsors) [-]
p.ct.etap2 = 0.8;                               % Secondary propulsors' propulsive efficiency in top-of-climb conditions (of ISOLATED propulsors) [-]
p.ct.Gamma = 0;                                 % Thrust vectoring in top-of-climb configuration [deg]
 
 
 
%% Mission/operational requirements for WP diagram
% Note: the throttle, phi and Phi values used to evaluate the constraints
% should be consistent with the power-control profiles specified in the MA.
% In future revisions, t, phi and Phi should automatically be selected from
% the MA_in structure (by e.g. evaluating the maximum and minimum per
% segment).
 
% Cruise
m.cr.h = 2400;                                 % Cruise altitude [m]
m.cr.M = 0.2797;                                  % Cruise Mach number [-]
m.cr.f = 0.999;                                  % Cruise weight fraction W/MTOW [-]
m.cr.t = 0.8;                                   % Cruise throttle setting P/P_max [-] (see note at end)
m.cr.phi = 0;                                 % Cruise supplied power ratio [-]
m.cr.Phi = 0;                                 % Cruise shaft power ratio [-]
 
% Landing
m.L.h = 0;                                      % Landing altitude [m]
m.L.f = 0.97;                                   % Landing weight fraction W/MTOW [-]
m.L.vs = 31.4;                                  % Stall speed requirement in landing conditions [m/s]
m.L.vApp = 1.23;                                % Stall margin during approach/landing, vApp/vs [-] (see Patterson 2017)
m.L.vAppIso = 1.05;                             % Stall margin of isolated wing during approach/landing, vApp/vsIso [-]
m.L.t = 0.5;                                      % Landing throttle setting P/P_max [-] (see note at end)
m.L.phi = 0.05;                                 % Landing supplied power ratio [-]
m.L.Phi = 0.05;                                 % Landing shaft power ratio [-]
 
% Take off
 
m.TO.h = 0;                                     % TO altitude [m]
m.TO.f = 1;                                     % TO weight fraction W/MTOW [-]
m.TO.s = 762;                                   % TO runway length [m]
m.TO.t = 1;                                     % TO throttle setting P/P_max [-] (see note at end)
m.TO.phi = 0.1;                                % TO supplied power ratio [-]
m.TO.Phi = 0.3;                                % TO shaft power ratio [-]
 
% OEI Balked landing
m.bL.G = 0.021;                                 % OEI balked landing climb gradient [-] (CS25.121d)
m.bL.f = 0.97;                                   % Max landing weight (MLW) as a fraction of MTOW [-]
m.bL.vMargin = 1.4;                             % Stall margin in balked-landing conditions
m.bL.t = 1;                                     % Balked landing throttle setting P/P_max [-] (see note at end)
m.bL.phi = 0.06;                                % Balked landing supplied power ratio [-]
m.bL.Phi = 0.05;                                   % Balked landing shaft power ratio [-]
 
% OEI ceiling
m.cI.h = 1200;                                  % OEI ceiling [m]
m.cI.f = 0.999;                                  % OEI-ceiling weight fraction W/MTOW [-]
m.cI.c = 0.5;                                   % Ceiling climb rate [m/s] (also used for cruise ceiling and cruise speed!)
m.cI.vMargin = 1.25;                            % Stall margin in OEI-ceiling conditions (also used for cruise ceiling) CHECK!
m.cI.t = 1;                                   % OEI-ceiling throttle setting P/P_max [-] (see note at end)
m.cI.phi = 0.02;                                % OEI-ceiling landing supplied power ratio [-]
m.cI.Phi = 0.1;                                   % OEI-ceiling landing shaft power ratio [-]
 
% Start-of-climb
m.cl.h = 0;                                     % Altitude for start-of-climb constraint [m]
m.cl.f = 1;                                     % Start-of-climb weight fraction W/MTOW [-]
m.cl.v = 37.4;                                    % Velocity at start-of-climb (shoud be equal to V2 obtained from TO constraint) [m/s]
m.cl.G = 0.02;                                   % Start-of-climb climb gradient [-] (based on MA observances)
m.cl.dVdt = 0.5;                                % Start-of-climb acceleration [m/s2]
m.cl.t = 1;                                   % Start-of-climb throttle setting P/P_max [-] (see note at end)
m.cl.phi = 0.07;                                % Start-of-climb supplied power ratio [-]
m.cl.Phi = 0.3;                                   % Start-of-climb shaft power ratio [-]
 
% Top-of-climb
m.ct.h = 2400;                                 % Altitude for top-of-climb constraint [m]
m.ct.f = 0.999;                                % Start-of-climb weight fraction W/MTOW [-]
m.ct.M = 0.2797;                               % Mach at top-of-climb (shoud be equal to cruise Mach) [-]
m.ct.G = 0.01;                                 % Top-of-climb climb gradient [-] (based on MA observances)
m.ct.dVdt = 0.05;                              % Top-of-climb acceleration [m/s2]
m.ct.t = 1;                                    % Top-of-climb throttle setting P/P_max [-] (see note at end)
m.ct.phi = 0.05;                               % Top-of-climb supplied power ratio [-]
m.ct.Phi = 0.3;                               % Top-of-climb shaft power ratio [-]
 
 
%% Mission Analysis input
 
% Mission characteristics
MA_in.PL = 363;                               % Payload [kg]
MA_in.R = 463000;                               % Range [m]
MA_in.R_div = 187515;                            % Diversion range [m]
MA_in.h_div = 914;                             % Diversion altitude [m]
MA_in.M_div = 0.207;                              % Diversion cruise Mach number [-]
                                                % Nominal mission M and h are specified in the
                                                % "m" structure (cruise constraint)
% Initial guesses for convergence loop
MA_in.OEM = 1347;                              % Operative empty mass incl. powertrain and wing, excl. bat [kg]
MA_in.FF_tot0 = 0.1;                            % Fuel fraction (excl. batteries, incl. diversion) of aircraft [-]
MA_in.FF_miss0 = 0.1;                           % Fuel fraction (excl. batteries, excl. diversion) of nominal mission [-]
MA_in.DOH_tot0 = 0.1;                           % Degree-of-hybridization, Ebat/(Efuel + Ebat) of aircraft [-]
MA_in.DOH_miss0 = 0.1;                          % Degree-of-hybridization, Ebat/(Efuel + Ebat) of nominal mission [-]
 
% Mission analysis power control settings. Linear interpolation used 
% between [start of segment, end of segment]. For climb and descent, 
% interpolation is carried out versus altitude, and for cruise, versus 
% range flown.
% Climb (all parameters must be specified; SEP is computed)
MA_in.cl.xi = [0.6 0.8];
MA_in.cl.phi = [0.1 0.05];
MA_in.cl.Phi = [0.1 0.1];
 
% Cruise (level flight is specified, so one DOF must be kept free)
MA_in.cr.xi = [NaN NaN];
MA_in.cr.phi = [0 0]; 
MA_in.cr.Phi = [0 0];
 
% Descent (all parameters must be specified; SEP is computed)
MA_in.de.xi = [0.01 0.01];
MA_in.de.phi = [0.01 0.01];
MA_in.de.Phi = [0.01 0.01];
 
% Diversion climb (all parameters must be specified; SEP is computed)
MA_in.Dcl.xi = [0.5 0.5];
MA_in.Dcl.phi = [0.05 0.05];
MA_in.Dcl.Phi = [0.05 0.05];
 
% Diversion cruise (level flight is specified, so one DOF must be kept free)
MA_in.Dcr.xi = [NaN NaN];
MA_in.Dcr.phi = [0 0]; 
MA_in.Dcr.Phi = [0 0]; 
 
% Diversion descent (all parameters must be specified; SEP is computed)
MA_in.Dde.xi = [0.01 0.01]; 
MA_in.Dde.phi = [0.01 0.1]; 
MA_in.Dde.Phi = [0.05 0.05]; 
 
 
%% Constants
 
c.g = 9.81;                                     % Gravity acceleration [m/s2]
c.rho_SL = 1.225;                               % Sea-level density [kg/m3]
c.T_SL = 288.15;                                % Sea-level temperature [K]
 
 
%% Program settings
 
% Design considerations
s.SelDes = 'minGT';                             % Selected design condition ('minWS','minGT',...)
s.Tcmax = 2;                                    % Maximum thrust coefficient (defined as Tc = T/rho/v^2/D^2) that 
                                                %   each individual propulsor should not surpass                                       
% Convergence settings                                              
s.n = 100;                                      % Number of points sampled per constraint [-]
s.TWmax = 1.0;                                  % Maximum thrust loading evaluated [-]
s.WSmax = 10000;                                 % Maxumum wing loading evaluated [N/W]
s.WSmin = 0;                                    % Minimum wing loading evaluated for landing constraint [N/W]
s.WPmax = 0.3;                                  % Maximum power loading shown in diagram [N/W]
s.itermax = 500;                                % Maximum number of iterations [-]
s.errmax = 1e-4;                                % Convergence criterion
s.NWS = 300;                                    % Number of wing loading points to sample when computing design point
s.rf = 0.6;                                     % Relaxation factor for convergence, recommended values [0.1 - 1.0]. 
                                                %   Lower RF = slower, but generally more chance of convergence
s.dt.cl = 15;                                   % Time step in MA, climb phase [s]
s.dt.cr = 40;                                   % Time step in MA, cruise phase [s]
s.dt.de = 20;                                   % Time step in MA, descent phase [s]
s.dt.Dcl = 5;                                   % Time step in MA, diversion climb phase [s]
s.dt.Dcr = 30;                                  % Time step in MA, diversion climb phase [s]
s.dt.Dde = 15;                                  % Time step in MA, diversion climb phase [s]
 
% Presentation of results                                                
s.levelString = [];                             % String inserted at the start of each displayed message
s.options = 0;                                  % Plot figures etc. in subroutines (careful with loops!)
s.figStart = 10;                                % Number of first figure generated
s.plotWPWS = 1;                                 % Plot WS-WP diagrams? 
s.plotMA = 1;                                   % Plot mission analysis graphs? 
s.plotPowertrain = 1;                           % Plot powertrain diagrams? 
s.plotTc = 1;                                   % Plot thrust coefficient constraints?
 
% Polar characteristics
s.Polar.plot = 1;                               % Plot aerodynamic polar?
s.Polar.TcInterval = [0 2];                     % Thrust coeff. interval sampled when creating aero polar
s.Polar.MInterval = [0 0.9];                    % Mach number interval sampled when creating aero polar
s.Polar.CLisoInterval = [0.3 1.8];              % Airframe lift coeff. interval sampled when creating aero polar. Avoid extremely low/high values
s.Polar.N = 20;                                 % Number of N and Tc points sampled in polar
s.Polar.N_CLiso = 100;                          % Number of CL_iso points sampled in polar
 
% Landing constraint check characteristics
s.LandingCheck.plot = 1;                        % Plot detailed landing constraint? 
s.LandingCheck.CL_map = 0.1:0.2:3.0;            % Isolated wing CL values plotted  during landing constraint check [-]
s.LandingCheck.CD0_map = [0.04:0.02:0.1 ...     % CD0 values plotted during landing constraint check [-]
                       0.15:0.05:0.4 0.5:0.1:1];
 
% Power-control envelope characteristics
s.Env.plot = 0;                                 % Plot power-control envelope? 
s.Env.Nphi = 50;                                % Number of phi/Phi values sampled
s.Env.Nxi = 40;                                 % Number of xi values sampled
s.Env.Nh = 30;                                  % Number of altitudes sampled
s.Env.con = 'cr';                               % Condition used for propulsive efficiency
s.Env.SPPH_phis = [0.05 0.2 NaN NaN];           % Constant phi values, for each of one an envelope is created
s.Env.SPPH_Phis = [NaN NaN 0.4 0.6];            % Constant Phi values, for each of one an envelope is created
                                                % The length of SPPH_phis same. For each element i, only ONE of the two can
                                                % be specified. The other must be NaN.
                                  
                                                
%% Functions & Dependencies
 
% Density lapse [kg/m3]
f.rho = @(h) c.rho_SL*(c.T_SL./(c.T_SL-0.0065*h)).^(1-9.81/287/0.0065);
 
% Speed of sound [m/s]
f.a = @(h) ((c.T_SL-0.0065*h)*1.4*287)^0.5;
 
% Take-off parameter correlation [N2/m2/W], s in [m] 
% First option: from AE1222-II course
% A = 0.0812;
% B = 8.52;
% Second option: from Raymer
f.TOP = @(s) 0.084958084*s+6.217443298;
        
% Normalized rotor sizing: D^2/W [m2/N]
f.D2W = @(WS,b_dp,N,dy,AR) (b_dp/N/(1+dy))^2*AR/WS;
 
% Thrust lapse: T_max/T_max_SL [-]
f.Alpha = @(rho) (rho/c.rho_SL)^0.75;
 
% Drag model: currently assuming a symmetric parabolic drag polar [-]
f.CD = @(CD0,CL_iso,AR,e) CD0 + CL_iso^2/(pi*AR*e);
 
% Weight correlations (per component instance!)
% Turboshaft weight in [kg] as a function of shaft
% power in [W], based on Roskam Part 5, Figure 6.2.
f.W.GT = @(P) 0.45359*10.^((log10(P/745.7)-0.011405)/1.1073);
 
% OEM in [kg] as a function of MTOM [kg], based on Roskam Part 1,
% Table 2.15 (in lb: 10.^((log10(MTOM)-AA)/BB))
% AA = 0.3774; BB = 0.9647;   % Regional turboprop aircraft
% AA = 0.0833; BB = 1.0383;   % Transport jets
% AA = 0.0966; BB = 1.0298;   % twin turbo prop 
f.W.OEM = @(MTOM) 0.45359*10.^((log10(MTOM/0.45359)-0.0966)/1.0298);
 
% Electrical machine weight in [kg] as a function of installed power in
% [W]. Using a constant power density for now.
f.W.EM = @(P) P/p.SP.EM;
 
% Ratio between SEP used to climb and total SEP (i.e. used to climb and
% accelerate). X will be in the interval [0,1].
f.SEPsplit = @(X) 0.83+0.1*cos(X*pi);
 
 
%% Notes
% Note regarding current nomenclature: the symbols used to designate 
% some of the variables has changed over the course of developing this
% code. This should be fixed throughout the code, but for now:
%   - Throttle is indicated with "t", or "xi" (in some subroutines). 
%   - "phi" is the supplied power ratio.
%   - "Phi" is the shaft power ratio, although this variable is represented 
%       with "Psi" in literature (see the paper of de Vries, Brown & Vos, 
%       2018)
%   - The thrust share provided by the propulsors (which is related to Phi)
%       is indicated using "T", which is not (only) the thrust of the 
%       aircraft, but the thrust produced by the DP propulsors divided by
%       the total thrust of the aircraft. In the paper this is represented 
%      with a "Chi".
%
% Note regarding the "throttle":
%   - For powertrain architectures containing a gas turbine (conventional,
%       turboelectric, serial, parallel, PTE or SPPH), this refers to the
%       throttle setting of the gas turbine, P_gt/P_gt_available, where
%       P_gt_available is the maximum power available in the given flight
%       condition.
%   - For fully electric architectures:
%       - For an e-1 or dual-e architecture, throttle refers to the
%           throttle setting of the primary electrical machine, P_gb/P_EM1,
%           where P_EM1 is the installed power of the primary electrical
%           machine. Since it is assumed to present no power lapse with
%           altitude or Mach, it is equal to the available power.
%       - For an e-2 architecture, throttle refers to the throttle setting
%           of the secondary electrical machine, P_s2/P_EM2, since this
%           layout contains no 

% end