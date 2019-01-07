function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintCruiseSpeed(a,m,p,f,s,c)
%% Cruise operating conditions

% No horizontal acceleration
dvdt = 0;  

% Climb rate specified: we want the aircraft to be able to reach cruise
% altitude (= top of climb), so it must also be able to climb at this
% speed. Using standard 100ft/min.
climb = [0 NaN];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];

% Freestream density, velocity and dynamic pressure
m.cr.rho = f.rho(m.cr.h);
m.cr.v = m.cr.M*f.a(m.cr.h);
q = 0.5*m.cr.rho*m.cr.v^2;


%% Thrust and (required) power loading 

% Intialize variables for loop
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP = NaN(size(WS));
a.cr.CL = NaN(size(WS));
a.cr.dCL = NaN(size(WS));
a.cr.dCDi = NaN(size(WS));
p.cr.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.cr.f;
    
    % Initial guess to speed up convergence
    TW_in = q*a.cr.CD0./WS_in + WS_in/pi/a.AR/a.cr.e/q;
    
    % Compute thrust and power loading
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc] = ComputeThrustLoading_vKnown(...
                        'cr',TW_in,WS_in,climb,maneuver,dvdt,a,m,p,f,s,c);
       
    % Correct to MTOW
    TW(i) = TW_out*m.cr.f;
    WP(i) = WP_out/m.cr.f;
     
    % Save aerodynamic variables
    a.cr.CL(i) = CL;
    a.cr.dCL(i) = dCL;
    a.cr.dCDi(i) = dCDi;
    p.cr.Tc(i) = Tc;
    p.cr.detap(i) = detap;
    
end


%% Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input.
% Keep indices to revert changes later
indices = isnan(WP);
WP(indices) = Inf;

% Call sizing routine
[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'cr',OEI,m,p,f,s);

% Change -Inf to +Inf to avoid warnings; Inf is a possible solution for
% zero power flow. Convert back to NaNs if necessary
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




