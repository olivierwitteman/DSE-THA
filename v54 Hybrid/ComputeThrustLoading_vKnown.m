function [WP,TW,CL,dCL,dCDi,dCD0,detap,Tc] = ComputeThrustLoading_vKnown(con,TW_in,WS,climb,maneuver,dvdt,a,m,p,f,s,c)
% This function computes the thrust (and power) loading for a given wing
% loading and flight conditions. It calculates useful (NOT shaft) power. 
% Powers and thrusts are NOT corrected to MTOW, SL or max throttle
% settings. 
%
% Input:
%   - con: Flight condition being evaluated. String which is used as name 
%       for the field within the structure variables (e.g. "TO", "cr")
%   - TW_in: initial guess for thrust loading (to speed up convergence)
%   - WS: wing loading [N/m2]
%   - climb: array of two elements, (1) climb rate [m/s] and (2) climb
%       gradient [-]. Only one can be specified, the other has to be NaN.
%       If both are NaN, zero climb is assumed (c = ROC = 0). If both have
%       values, an error is returned.
%   - maneuver: array of three elements, (1) bank angle [deg], (2) load
%       factor [-] and (3) turn radius [m]. Only one can be specified, the 
%       other two have to be NaN. If all are NaN, zero bank angle is 
%       assumed (mu = 0, n = 1, r = Inf). If two or more are numbers, an
%       error is returned.
%   - dvdt: horizontal acceleration of aircraft [m/s2]
%   - a,m,p,s,f,c: structures containing aicraft parameters, mission
%       parameters, powertrain parameters, programme settings, anonymous
%       functions, and constants, respectively. See input of
%       WP_WS_diagram_DP for more info.
%
% Output:
%   - WP: required flight power loading [N/W]
%   - TW: required flight thrust loading [-]
%   - CL: isolated wing lift coefficient [-]
%   - dCL: wing lift coefficient increase due to DP [-]
%   - dCDi: wing (thrust-induced) drag coefficient increase due to DP [-]
%   - dCD0: wing parasite drag coefficient increase due to DP [-]
%   - detap: change in propulsive efficiency due to DP [-]
%   - Tc: propulsors' thrust coefficient [-]
%
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 23-08-17
%%% Last modified: 05-02-18


%% Input 

% Retrieve input from structures
CD0 = a.(con).CD0; AR = a.AR; e = a.(con).e; Lambda = a.Lambda;
rho = m.(con).rho; v = m.(con).v; M = m.(con).M;
Gamma = p.(con).Gamma; T = p.(con).T; xp = p.xp; b_dp = p.b_dp; ...
    dy = p.dy; N = p.N; etap = p.(con).etap;
itermax = s.itermax; errmax = s.errmax; options = s.options; rf = s.rf;
g = c.g;

% Disk loading [m2/N]
D2W = f.D2W(WS,b_dp,N,dy,AR);

% Check input values for climb rate/gradient
if length(climb)~=2
    error(['The inpute array for "climb" has to have two elements,'...
           'one of which has to be NaN. Specify either climb rate '...
           'or climb gradient.'])
else
    % If both elements are NaNs
    if sum(isnan(climb))==2
        c = 0;
        
    % If only one element is specified    
    elseif sum(isnan(climb))==1
        
        % If ROC is specified
        if ~isnan(climb(1))
            c = climb(1);
            
        % If climb gradient is specified
        elseif ~isnan(climb(2))
            c = climb(2)*v;
        end
    
    % If both elements have values
    else
        error(['The inpute array for "climb" has to have two elements,'...
               'one of which has to be NaN. Specify either climb rate '...
               'or climb gradient.'])
    end
end

% Check input values for maneuvers
if length(maneuver)~=3
    error(['The inpute array for "maneuver" has to have three elements,'...
           'two of which have to be NaN. Specify either bank angle, '...
           'load factor or turn radius.'])
else
    % If all elements are NaNs
    if sum(isnan(maneuver))==3
        mu = 0;
    
    % If only one element is specified    
    elseif sum(isnan(maneuver))==2
        
        % If bank angle is specified
        if ~isnan(maneuver(1))
            computeMu = 1;
            mu = maneuver(1);
            
        % If load factor is specified
        elseif ~isnan(maneuver(2))
            computeMu = 2;
            n = maneuver(2);
            
        % If turn radius is specified
        elseif ~isnan(maneuver(3))
            computeMu = 3;
            r = maneuver(3);
        end
    
    % If two or more elements have values
    else
        error(['The inpute array for "maneuver" has to have three '...
               'elements, one of which has to be NaN. Specify either '...
               'bank angle, load factor or turn radius.'])
    end
end


%% Compute equilibrium flight point

% Dynamic pressure
q = 0.5*rho*v^2;

% Initial values
TW0 = TW_in;                                  
CL0 = WS/q;
err = 1;                                           
iter = 0;                                          

% Start convergence loop on TW, CL
while err > errmax
    iter = iter+1;
    
    % Compute thrust coefficient, defined as Tc = T/(rho*v^2*D^2)
    Tc = T*TW0/N/rho/v^2/D2W;
    
    % Compute prop radius/wing chord ratio
    Rc = 0.5*(D2W*WS*AR)^0.5;

    % Compute delta terms
    [dCl,dCd0,dCdi,detap] = WingPropDeltas(Tc,Rc,xp,AR,CL0,M,Lambda,...
                                            Gamma,etap,options);
    dCL = dCl*b_dp;
    dCD0 = dCd0*b_dp;
    dCDi = dCdi*b_dp;
       
    % Compute bank angle if necessary
    if computeMu == 2
        mu = rad2deg(acos((1-(c/v)^2)^0.5 / (n + T*TW0*sind(Gamma))));
    elseif computeMu == 3
        mu = rad2deg(asin(v^2/r/g/...
            (q/WS*(CL0 + dCl) + T*TW0*sind(Gamma))));
    end
    
    % Compute updated isolated wing lift coefficient (2nd law of Newton
    % along Z' direction)
    CL1 = WS/q*(1/cosd(mu)*(1-(c/v)^2)^0.5 - T*TW0*sind(Gamma)) - dCL;
    
    % Compute updated thrust loading (2nd law of Newton along X' direction)
    TW1 = 1/(1-T+T*cosd(Gamma))*(c/v + dvdt/g + ...
            q/WS*(f.CD(CD0,CL1,AR,e) + dCD0 + dCDi));
    
    % Unrealistic input values the model may return negative,
    %   and afterwards complex, values, which will cause
    %   erroneous solutions. In this case provide NaN as answer
    CL_tot = CL1 + dCL;
    if CL_tot < 0                                   
        CL1 = NaN;                               
    end                                             
    
    % Calculate error
    err1 = abs(TW1-TW0)/abs(TW0);
    err2 = abs(CL1-CL0)/abs(CL0);
    err = err1 + err2;
        
    % Update values, add small offset to avoid reiterating NaN values. Use
    % relaxation factor to avoid divergence
    if isnan(TW1) || isnan(CL1)                                  
        TW0 = 1.1*TW_in;
        CL0 = 1.1*WS/q;
    else
        TW0 = TW0+rf*(TW1-TW0);
        CL0 = CL0+rf*(CL1-CL0);
    end
    
    % Limit number of iterations
    if iter >= itermax
        TW1 = NaN;
        CL1 = NaN;
        break
    end
    
end

% Compute required (flight, not shaft) power loading
WP = 1/TW1/v;

% Assign output
TW = TW1;
CL = CL1;



