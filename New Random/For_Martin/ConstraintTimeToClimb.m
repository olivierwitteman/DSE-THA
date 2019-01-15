function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintTimeToClimb(a,m,p,f,s,c)
%% Cruise ceiling operating conditions

% Compute rate of climb required at bottom of climb. Assume constant climb
% rate
m.cl.c = (m.cr.h-m.cl.h)/m.cl.ttc;

% Freestream density, velocity and dynamic pressure
m.cl.rho = f.rho(m.cl.h);

% No horizontal acceleration
dvdt = 0;  

% Given climb rate (same as OEI climb rate)
climb = [m.cl.c NaN];       

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN];


%% Compute power loading

% Initialize variables
WS = linspace(0,s.WSmax,s.n);
TW = NaN(size(WS));
WP = NaN(size(WS));
a.cl.CL = NaN(size(WS));
a.cl.dCL = NaN(size(WS));
a.cl.dCDi = NaN(size(WS));
p.cl.Tc = NaN(size(WS));

% Loop over WS values
for i = 1:length(WS)
    
    % Wing loading in flight condition
    WS_in = WS(i)*m.cl.f;
        
    % Initial guess for flight speed
    vr = m.cl.vMargin*(WS_in*2/m.cl.rho/a.cl.CLmax)^0.5;
    
    % Initial guess for dynamic pressure
    q = 0.5*m.cl.rho*vr^2;

    % Initial guess to speed up convergence
    TW_in = q*a.cl.CD0./WS_in + WS_in/pi/a.AR/a.cl.e/q;

    % Compute power loading and flight speed
    [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc,v] = ...
        ComputeThrustLoading_vMargin('cl',TW_in,WS_in,climb,...
                                     maneuver,dvdt,a,m,p,f,s,c);
    % Correct to MTOW
    TW(i) = TW_out*m.cl.f;
    WP(i) = WP_out/m.cl.f;
    
    % Save aerodynamic variables
    a.cl.CL(i) = CL;
    a.cl.dCL(i) = dCL;
    a.cl.dCDi(i) = dCDi;
    p.cl.Tc(i) = Tc;
    m.cl.v(i) = v;
    p.cl.detap(i) = detap;
    
end

% Remove imaginary numbers, which may appear at extreme wing loading values
if ~isreal(WP)
    array = [];
    for i = 1:length(WP)
        if ~isreal(WP(i))
            array = [array i];
        end
    end
    disp([s.levelString ...
            ' > Detected imaginary WP values at WS indices ['...
            num2str(array) ']. Selecting real part only.'])
    WP = real(WP);
    TW = real(TW);
end


%% Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input.
% Keep indices to revert changes later
indices = isnan(WP);
WP(indices) = Inf;

% Call sizing routine
[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'cl',OEI,m,p,f,s);

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



