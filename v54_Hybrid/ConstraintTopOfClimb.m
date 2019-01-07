function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintTopOfClimb(a,m,p,f,s,c)
%% Operating conditions at start-of-climb

% Horizontal acceleration
dvdt = m.ct.dVdt;  

% Climb gradient specified
climb = [NaN m.ct.G];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];

% Freestream density, velocity and dynamic pressure
m.ct.rho = f.rho(m.ct.h);
m.ct.v = m.ct.M*f.a(m.ct.h);
q = 0.5*m.ct.rho*m.ct.v^2;


%% Thrust and (required) power loading 

% Intialize variables for loop
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP = NaN(size(WS));
a.ct.CL = NaN(size(WS));
a.ct.dCL = NaN(size(WS));
a.ct.dCDi = NaN(size(WS));
p.ct.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.ct.f;
    
    % Initial guess to speed up convergence
    TW_in = q*a.ct.CD0./WS_in + WS_in/pi/a.AR/a.ct.e/q;
    
    % Compute thrust and power loading
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc] = ComputeThrustLoading_vKnown(...
                        'ct',TW_in,WS_in,climb,maneuver,dvdt,a,m,p,f,s,c);
       
    % Correct to MTOW
    TW(i) = TW_out*m.ct.f;
    WP(i) = WP_out/m.ct.f;
     
    % Save aerodynamic variables
    a.ct.CL(i) = CL;
    a.ct.dCL(i) = dCL;
    a.ct.dCDi(i) = dCDi;
    p.ct.Tc(i) = Tc;
    p.ct.detap(i) = detap;
    
end


%% Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input.
% Keep indices to revert changes later
indices = isnan(WP);
WP(indices) = Inf;

% Call sizing routine
[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'ct',OEI,m,p,f,s);

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




