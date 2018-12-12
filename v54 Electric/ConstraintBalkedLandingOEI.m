function [WS1,WS2,WP_path1,WP_comp1,WP_loss1,WP_path2,WP_comp2,WP_loss2,TW,a,m,p] = ConstraintBalkedLandingOEI(a,m,p,f,s,c)
%% Balked landing operating conditions

% No horizontal acceleration
dvdt = 0;  

% Given climb gradient
climb = [NaN m.bL.G];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];

% Assign variables shared with landing condition
a.bL.e = a.L.e;
a.bL.CLmax = a.L.CLmax;
m.bL.rho = m.L.rho;
m.bL.h = m.L.h;


%% Compute power loading

% Initialize variables
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP2 = NaN(size(WS));
a.bL.CL = NaN(size(WS));
a.bL.dCL = NaN(size(WS));
a.bL.dCDi = NaN(size(WS));
p.bL.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.bL.f;
    
    % Initial guess for flight speed
    % CS25.121d: Flight speed has to be 1.1*1.4 times V_SR, the reference stall
    % speed. The reference stall speed is taken to be equal to the stall speed
    % in landing configuration for  the balked landing, since in essence the AC
    % is in landing configuration.
    m.bL.vr = m.bL.vMargin*(WS_in*2/m.L.rho/a.bL.CLmax)^0.5;
    
    % Initial guess for dynamic pressure
    q = 0.5*m.L.rho*m.bL.vr^2;
         
    % Initial guess to speed up convergence
    TW_in = q*a.bL.CD0./WS_in + WS_in/pi/a.AR/a.cr.e/q;

    % Compute power loading and flight speed
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc,v] = ...
        ComputeThrustLoading_vMargin('bL',TW_in,WS_in,climb,...
                                     maneuver,dvdt,a,m,p,f,s,c);
    
    % Correct to MTOW
    TW(i) = TW_out*m.bL.f;
    WP2(i) = WP_out/m.bL.f;
    
    % Save aerodynamic variables
    a.bL.CL(i) = CL;
    a.bL.dCL(i) = dCL;
    a.bL.dCDi(i) = dCDi;
    p.bL.Tc(i) = Tc;
    m.bL.v(i) = v;
    p.bL.detap(i) = detap;
    
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
[WP_path2,WP_comp2,WP_loss2] = SizePowertrain(WP2,'bL',OEI,m,p,f,s);

% Switch -Infs
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
[WP_path1,WP_comp1,WP_loss1] = SizePowertrain(WP1,'bL',OEI,m,p,f,s);

% Switch -Infs
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




