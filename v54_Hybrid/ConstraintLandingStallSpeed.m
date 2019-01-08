function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintLandingStallSpeed(a,m,p,f,s,c)
%% Landing operating conditions

% Density
m.L.rho = f.rho(m.L.h);

% Velocity/dynamic pressure
m.L.v = m.L.vs;
m.L.M = m.L.v/f.a(m.L.h);
q = 0.5*m.L.rho*m.L.v^2; 
           
% No horizontal acceleration
dvdt = 0;               

% No climb rate, no climb gradient
climb = [0 NaN];        

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN]; 


%% Start loop to compute power required for cruise speed with DP enabled
WS = linspace(s.WSmin,s.WSmax,s.n);
TW = NaN(size(WS));
WP = NaN(size(WS));
a.L.CL = NaN(size(WS));
a.L.dCL = NaN(size(WS));
a.L.dCDi = NaN(size(WS));
p.L.Tc = NaN(size(WS));
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.L.f;
    
    % Initial guess to speed up convergence
    TW_in = q*a.L.CD0./WS_in + WS_in/pi/a.AR/a.L.e/q;
    
    % Compute thrust and power loading
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc] = ComputeThrustLoading_vKnown(...
                        'L',TW_in,WS_in,climb,maneuver,dvdt,a,m,p,f,s,c);

    % Correct to MTOW
    TW(i) = TW_out*m.L.f;
    WP(i) = WP_out/m.L.f;
       
    % If the wing loading is unrealistically high and it is starting to
    % diverge, stop
    if TW(i) > 10e2 
        TW(i) = NaN;
        WP(i) = NaN;
        break
    end
    
    % Save aerodynamic variables
    a.L.CL(i) = CL;
    a.L.dCL(i) = dCL;
    a.L.dCDi(i) = dCDi;
    p.L.Tc(i) = Tc;
    p.L.detap(i) = detap;
    
end

% For unrealistically high WS values the solution will not converge,
% remove indexes where any of the results has NaN value
NaNarray = sum([TW; WP; a.L.CL; a.L.dCL; a.L.dCDi; p.L.Tc],1);
TW = TW(~isnan(NaNarray));
WP = WP(~isnan(NaNarray));
WS = WS(~isnan(NaNarray));
a.L.CL = a.L.CL(~isnan(NaNarray));
a.L.dCL = a.L.dCL(~isnan(NaNarray));
a.L.dCDi = a.L.dCDi(~isnan(NaNarray));
p.L.Tc = p.L.Tc(~isnan(NaNarray));

WS_store = WS;
TW_store = TW;
WP_store = WP;
CL_store = a.L.CL;

% Interpolate at isolated wing CLmax value to obtain max wing loading.
WS = interp1(a.L.CL,WS,a.L.CLmax,'linear');
TW = interp1(a.L.CL,TW,a.L.CLmax,'linear');
WP = interp1(a.L.CL,WP,a.L.CLmax,'linear');
a.L.dCL = interp1(a.L.CL,a.L.dCL,a.L.CLmax,'linear');
a.L.dCDi = interp1(a.L.CL,a.L.dCDi,a.L.CLmax,'linear');
p.L.Tc = interp1(a.L.CL,p.L.Tc,a.L.CLmax,'linear');
a.L.CL = a.L.CLmax;


%%  Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input,
% but save indices to revert later
indices = isnan(WP);
WP(indices) = Inf;

% Call sizing function
[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'L',OEI,m,p,f,s);

% Change -Inf to +Inf to avoid warnings; Inf is a possible solution for
% zero power flow. Also revert back to NaN if necessary
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







