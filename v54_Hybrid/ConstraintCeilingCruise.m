function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintCeilingCruise(a,m,p,f,s,c)
%% Cruise ceiling operating conditions

% Settings are the same as cruise speed constraint, but without vel/M
% specified, and with CLmax specified instead.
m.cc = m.cr;
m.cc.v = [];
m.cc.M = [];
p.cc = p.cr;
a.cc = a.cr;

% Climb rate at ceiling should be standard 100ft/min, same as OEI ceiling
% Max lift coefficient is the clean wing CLmax, same as OEI ceiling. The
% velocity margin is also the same as OEI ceiling.
a.cc.CLmax = a.cI.CLmax;
m.cc.c = m.cI.c;
m.cc.vMargin = m.cI.vMargin;

% Freestream density, velocity and dynamic pressure
m.cc.rho = f.rho(m.cc.h);

% No horizontal acceleration
dvdt = 0;  

% Given climb rate (same as OEI climb rate)
climb = [m.cc.c NaN];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];


%% Compute power loading

% Initialize variables
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP = NaN(size(WS));
a.cc.CL = NaN(size(WS));
a.cc.dCL = NaN(size(WS));
a.cc.dCDi = NaN(size(WS));
p.cc.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.cc.f;
      
    % Initial guess for flight speed
    vr = m.cI.vMargin*(WS_in*2/m.cc.rho/a.cc.CLmax)^0.5;
    
    % Initial guess for dynamic pressure
    q = 0.5*m.cc.rho*vr^2;
    
    % Initial guess to speed up convergence
    TW_in = q*a.cc.CD0./WS_in + WS_in/pi/a.AR/a.cc.e/q;

    % Compute power loading and flight speed
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc,v] = ...
        ComputeThrustLoading_vMargin('cc',TW_in,WS_in,climb,...
                                     maneuver,dvdt,a,m,p,f,s,c);
    
    % Correct to MTOW
    TW(i) = TW_out*m.cc.f;
    WP(i) = WP_out/m.cc.f;
    
    % Save aerodynamic variables
    a.cc.CL(i) = CL;
    a.cc.dCL(i) = dCL;
    a.cc.dCDi(i) = dCDi;
    p.cc.Tc(i) = Tc;
    m.cc.v(i) = v;
    p.cc.detap(i) = detap;
    
end


%% Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input
indices = isnan(WP);
WP(indices) = Inf;

[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'cc',OEI,m,p,f,s);

% Change -Inf to +Inf to avoid warnings; Inf is a possible solution for
% zero power flow
pathnames = fieldnames(WP_path);
for i = 1:size(pathnames,1)
    WP_path.(pathnames{i})(WP_path.(pathnames{i})==-Inf) = Inf;
    WP_path.(pathnames{i})(indices) = NaN;
end
compnames = fieldnames(WP_comp);
for i = 1:size(compnames,1)
    WP_loss.(compnames{i})(WP_loss.(compnames{i})==-Inf) = Inf;
    WP_comp.(compnames{i})(WP_comp.(compnames{i})==-Inf) = Inf;
    WP_loss.(compnames{i})(indices) = NaN;
    WP_comp.(compnames{i})(indices) = NaN;
end



