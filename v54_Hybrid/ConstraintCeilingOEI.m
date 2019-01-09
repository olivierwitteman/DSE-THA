function [WS1,WS2,WP_path1,WP_comp1,WP_loss1,WP_path2,WP_comp2,WP_loss2,TW,a,m,p] = ConstraintCeilingOEI(a,m,p,f,s,c)
%% OEI-2 ceiling operating conditions

% Freestream density, velocity and dynamic pressure
m.cI.rho = f.rho(m.cI.h);

% No horizontal acceleration
dvdt = 0;  

% Given climb rate
climb = [m.cI.c NaN];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];


%% Compute power loading

% Initialize variables
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP2 = NaN(size(WS));
a.cI.CL = NaN(size(WS));
a.cI.dCL = NaN(size(WS));
a.cI.dCDi = NaN(size(WS));
p.cI.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.cI.f;
      
    % Initial guess for flight speed
    vr = m.cI.vMargin*(WS_in*2/m.cI.rho/a.cI.CLmax)^0.5;
    
    % Initial guess for dynamic pressure
    q = 0.5*m.cI.rho*vr^2;
    
    % Initial guess to speed up convergence
    TW_in = q*a.cI.CD0./WS_in + WS_in/pi/a.AR/a.cI.e/q;

    % Compute power loading and flight speed
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc,v] = ...
        ComputeThrustLoading_vMargin('cI',TW_in,WS_in,climb,...
                                     maneuver,dvdt,a,m,p,f,s,c);
    
    % Correct to MTOW
    TW(i) = TW_out*m.cI.f;
    WP2(i) = WP_out/m.cI.f;
    
    % Save aerodynamic variables
    a.cI.CL(i) = CL;
    a.cI.dCL(i) = dCL;
    a.cI.dCDi(i) = dCDi;
    p.cI.Tc(i) = Tc;
    m.cI.v(i) = v;
    p.cI.detap(i) = detap;
end


%% Compute HEP component power loading: OEI in secondary powertrain

% Assign same required power for both cases
WP1 = WP2;

% Secondary powertrain failure
OEI = 2;

% Replace NaNs with Inf's so that the code understands that WP is an input
indices = isnan(WP2);
WP2(indices) = Inf;

% Evaluate powertrain model
[WP_path2,WP_comp2,WP_loss2] = SizePowertrain(WP2,'cI',OEI,m,p,f,s);

% Switch back to NaNs
pathnames = fieldnames(WP_path2);
for i = 1:size(pathnames,1)
    WP_path2.(pathnames{i})(WP_path2.(pathnames{i})==-Inf) = Inf;
    WP_path2.(pathnames{i})(indices) = NaN;
end
compnames = fieldnames(WP_comp2);
for i = 1:size(compnames,1)
    WP_loss2.(compnames{i})(WP_loss2.(compnames{i})==-Inf) = Inf;
    WP_comp2.(compnames{i})(WP_comp2.(compnames{i})==-Inf) = Inf;
    WP_loss2.(compnames{i})(indices) = NaN;
    WP_comp2.(compnames{i})(indices) = NaN;
end

% Primary powertrain failure
OEI = 1;

% Replace NaNs with Inf's so that the code understands that WP is an input
indices = isnan(WP1);
WP1(indices) = Inf;

% Evaluate powertrain model
[WP_path1,WP_comp1,WP_loss1] = SizePowertrain(WP1,'cI',OEI,m,p,f,s);

% Switch back to NaNs
pathnames = fieldnames(WP_path1);
for i = 1:size(pathnames,1)
    WP_path1.(pathnames{i})(WP_path1.(pathnames{i})==-Inf) = Inf;
    WP_path1.(pathnames{i})(indices) = NaN;
end
compnames = fieldnames(WP_comp2);
for i = 1:size(compnames,1)
    WP_loss1.(compnames{i})(WP_loss1.(compnames{i})==-Inf) = Inf;
    WP_comp1.(compnames{i})(WP_comp1.(compnames{i})==-Inf) = Inf;
    WP_loss1.(compnames{i})(indices) = NaN;
    WP_comp1.(compnames{i})(indices) = NaN;
end

% Duplicate wing loading output
WS1 = WS;
WS2 = WS;


