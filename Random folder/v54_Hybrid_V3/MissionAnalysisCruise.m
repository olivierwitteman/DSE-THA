function [A, Sw_give] = MissionAnalysisCruise(a,p,f,s,c,con,aircraft,A,R_target,dt) % added the Sw_give
%% Initialize variables

% Loop counter
k = 1;

% Initial guess for delta eta_p
detap = 0;

% Aircraft characteristics
D2 = aircraft.D2;
Sw = aircraft.Sw;
Sw_give = Sw % added
Rc = aircraft.Rc;

% Powertrain component efficiencies
etas.GT = p.eta_GT;     % Gas turbine
etas.GB = p.eta_GB;     % Gearbox
etas.EM1 = p.eta_EM1;   % Primary electrical machine
etas.PM = p.eta_PM;     % PMAD (power management and distribution)
etas.EM2 = p.eta_EM2;   % Secondary electrical machine

% Gravity
g = c.g;

% Operating conditions
h = A.h(1);
rho = f.rho(h);
v = A.v(1);
M = v/f.a(h);
q = 0.5*rho*v^2;
G = 0;

% Initial guess for required propulsive power
Pp = 0.5*rho*v^3*Sw*(a.(con).CD0 + 0.5^2/pi/a.AR/a.(con).e);

% Compute available power at altitude
if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
    Pa = aircraft.Pdes.EM1;
elseif strcmp(p.config,'e-2')
    Pa = aircraft.Pdes.EM2;
else
    Pa = aircraft.Pdes.GTM*f.Alpha(rho);
end

% First value of array
A.dhdt(1) = 0;
A.M(1) = M;


%% Loop over time steps until required range is reached
while A.R(k)-A.R(1) < R_target
    
    % Update weight & range
    W = A.W(k);
    R = A.R(k);
    
    % Select power-control parameters
    if R-A.R(1) >= 0 && R-A.R(1) <= R_target
    xi = interp1([0 R_target],A.xi,R-A.R(1),'linear');
    phi = interp1([0 R_target],A.phi,R-A.R(1),'linear');
    Phi = interp1([0 R_target],A.Phi,R-A.R(1),'linear');
    else
        error('Range out of bounds')
    end
    chi = p.(con).T;
  
    % Calculate lift required and break down into delta CL and CLiso
    if k == 1
        CL_iso0 = 0.5;
    else
        CL_iso0 = A.aero.CLiso(k-1);
    end
    err = 1;
    while err > s.errmax
        
        % Compute thrust
        T = Pp/v;
        Tc = chi*T/p.N/rho/v^2/D2;
        
        % Calculate required total lift coefficient
        CL = W/Sw/q*((1-G^2)-chi*sind(p.(con).Gamma)*T/W);
        
        % Obtain deltas
        [dCl,dCd0,dCdi,detap] = WingPropDeltas(Tc,Rc,p.xp,a.AR,CL_iso0,...
            M,a.Lambda,p.(con).Gamma,p.(con).etap,0);
        dCL = dCl*p.b_dp;
        dCD0 = dCd0*p.b_dp;
        dCDi = dCdi*p.b_dp;
        
        % Check convergence of lift coefficient and update
        CL_iso = CL-dCL;
        err = abs(CL_iso-CL_iso0)/CL_iso0;
        CL_iso0 = CL_iso;
    end
    
    % Compute drag
    D = (a.(con).CD0 + dCD0 + CL_iso^2/pi/a.AR/a.(con).e + dCDi)*q*Sw;
    
    % Update required flight power
    Pp = D*v;
        
    % Update propulsive efficiency
    if p.DP == 2
        etas.P1 = p.(con).etap1;
        etas.P2 = p.(con).etap2+detap;
    elseif p.DP == 1
        etas.P1 = p.(con).etap1+detap;
        etas.P2 = p.(con).etap2;
    else
        etas.P1 = p.(con).etap1;
        etas.P2 = p.(con).etap2;
    end
    
    % Calculate how much power is required from each component
    [P_out,~,phi_out,Phi_out,~] = PowerTransmissionComputation_v2...
        (p.config,etas,phi,Phi,xi,Pp,[],[],Pa);
    
    % Compute xi_out seperately 
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
        xi_out = -P_out.e1/Pa;
    elseif strcmp(p.config,'e-2')
        xi_out = P_out.e2/Pa;
    else
        xi_out = P_out.gt/Pa;
    end

    % Deltas corresponding to this timestep
    deltah = 0;
    deltav = 0;
    deltaR = v*dt;
    deltaEbat = -P_out.bat*dt;% Energy consumed: negative if discharging
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e') ||...
            strcmp(p.config,'e-2')
        deltaEf = 0;
    else
        deltaEf = -P_out.f*dt;
    end
    deltaW = g*deltaEf/p.SE.f;
    
    % Update data for next iteration
    A.t(k+1) = A.t(k) + dt;
    A.h(k+1) = A.h(k) + deltah;
    A.W(k+1) = A.W(k) + deltaW;
    A.v(k+1) = A.v(k) + deltav;
    A.G(k+1) = 0;
    A.dVdt(k+1) = 0;
    A.dhdt(k+1) = 0;
    A.R(k+1) = A.R(k) + deltaR;
    A.M(k+1) = A.v(k+1)/f.a(A.h(k+1));
    A.Ebat(k+1) = A.Ebat(k) + deltaEbat;
    A.Ef(k+1) = A.Ef(k) + deltaEf;
    
    % Save aerodynamic variables
    A.aero.CLiso(k) = CL_iso;
    A.aero.CL(k) = CL_iso + dCL;
    A.aero.CDiso(k) = a.(con).CD0 + CL_iso^2/pi/a.AR/a.(con).e;
    A.aero.CD(k) = a.(con).CD0 + dCD0 + CL_iso^2/pi/a.AR/a.(con).e + dCDi;
    A.aero.LD(k) = A.aero.CL(k)/A.aero.CD(k);
    A.aero.etap_DP(k) = p.(con).etap + detap;
    
    % Save powers
    A.P.x(k) = NaN;
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e') ||...
            strcmp(p.config,'e-2')
        A.P.availableGT(k) = NaN;
    else
        A.P.availableGT(k) = Pa;
    end
    A.P.drag(k) = D*v;
    A.P.acceleration(k) = 0;
    A.P.climb(k) = 0;
    if k == 1; names = fieldnames(P_out); end
    for j = 1:size(names,1)
        A.P.(names{j})(k) = P_out.(names{j});
    end
    
    % Save power-control variables
    A.P.xi(k) = xi_out;
    if strcmp(p.config,'serial') || strcmp(p.config,'parallel') ||...
            strcmp(p.config,'SPPH')
        A.P.phi(k) = phi_out;
    else
        A.P.phi(k) = NaN;
    end
    if strcmp(p.config,'PTE') || strcmp(p.config,'SPPH') || ...
            strcmp(p.config,'dual-e')
        A.P.Phi(k) = Phi_out;
    else
        A.P.Phi(k) = NaN;
    end
    
    % Update counter, avoid excessive iterations
    k = k+1;
    if k > 5*s.itermax
        if strcmp(con,'Dcr')
            error('Could not divert at given altitude/Mach number')
        else
            error('Could not cruise at given altitude/Mach number')
        end
    end
end

% Add last element to aerodynamic/power variables to make plotting arrays 
% same length
aeronames = fieldnames(A.aero);
for i = 1:size(aeronames,1)
    A.aero.(aeronames{i})(k) = NaN;
end
Pnames = fieldnames(A.P);
for i = 1:size(Pnames,1)
    A.P.(Pnames{i})(k) = NaN;
end
