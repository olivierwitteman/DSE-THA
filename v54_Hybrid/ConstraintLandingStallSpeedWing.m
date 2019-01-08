function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintLandingStallSpeedWing(a,m,p,f,s,c)
%% Landing operating conditions

% Altitude
m.Liso.h = m.L.h;

% Density
m.Liso.rho = f.rho(m.Liso.h);

% Approach speed
v_approach = m.L.vs*m.L.vApp;       

% Speed at which isolated wing (no DP effect) must stall
m.Liso.v =3* v_approach/m.L.vAppIso; 

% Dynamic pressure
q = 0.5*m.Liso.rho*m.Liso.v^2; 

% Maximum Lift coefficient
a.Liso.CL = a.L.CLmax;
a.Liso.CLmax = a.L.CLmax;


%% Compute required wing loading

% Account for weight fraction 
WS = a.L.CLmax*q/m.L.f;

% Other output variables
TW = NaN;
WP_path.p = NaN;
WP_path.p1 = NaN;
WP_path.p2 = NaN;
WP_path.s2 = NaN;
WP_path.s1 = NaN;
WP_path.e2 = NaN;
WP_path.f = NaN;
WP_path.bat = NaN;
WP_path.gt = NaN;
WP_path.gtm = NaN;
WP_path.e1 = NaN;
WP_path.gb = NaN;
WP_comp.GT = NaN;
WP_comp.GTM = NaN;
WP_comp.GB = NaN;
WP_comp.P1 = NaN;
WP_comp.P2 = NaN;
WP_comp.PM = NaN;
WP_comp.EM1 = NaN;
WP_comp.EM2 = NaN;
WP_comp.EM1M = NaN;
WP_comp.EM2M = NaN;
WP_loss = WP_comp;


