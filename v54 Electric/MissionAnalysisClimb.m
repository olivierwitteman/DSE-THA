function [A] = MissionAnalysisClimb(a,p,f,s,c,con,aircraft,A,h_target,M_target,dt)
%% Initialize variables

% Loop counter
k = 1;

% Initial guess for delta eta_p
detap = 0;

% SEP split profile: should be optimized!
% x_array = 0.8+0.1*cos(linspace(0,pi,100));
h_array = linspace(0,h_target,100);
x_array = f.SEPsplit(h_array/max(h_array));

% Aircraft characteristics
D2 = aircraft.D2;
Sw = aircraft.Sw;
Rc = aircraft.Rc;

% Powertrain component efficiencies
etas.GT = p.eta_GT;     % Gas turbine
etas.GB = p.eta_GB;     % Gearbox
etas.EM1 = p.eta_EM1;   % Primary electrical machine
etas.PM = p.eta_PM;     % PMAD (power management and distribution)
etas.EM2 = p.eta_EM2;   % Secondary electrical machine

% Gravity
g = c.g;

% First value of array
A.dhdt(1) = A.G(1)*A.v(1);
A.M(1) = A.v(1)/f.a(A.h(1));


%% Loop over time steps until cruise altitude and speed are reached
while (A.h(k) < h_target || A.v(k) < M_target*f.a(h_target)) 
    
    % Operating conditions
    h = A.h(k);
    v = A.v(k);
    rho = f.rho(h);
    M = v/f.a(h);
    q = 0.5*rho*v^2;
    G = A.G(k);
    W = A.W(k);
    
    % Select power-control parameters
    if h >= 0 && h <= h_target
        xi = interp1([0 h_target],A.xi,h,'linear');
        phi = interp1([0 h_target],A.phi,h,'linear');
        Phi = interp1([0 h_target],A.Phi,h,'linear');
    else
        error('Altitude out of bounds')
    end
    chi = p.(con).T;
    
    % Compute available power at altitude
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
        Pa = aircraft.Pdes.EM1;
    elseif strcmp(p.config,'e-2')
        Pa = aircraft.Pdes.EM2;
    else
        Pa = aircraft.Pdes.GTM*f.Alpha(rho);
    end
    
    % Calculate lift required and break down into delta CL and CLiso
    if k == 1
        CL_iso0 = 0.5;
    else
        CL_iso0 = A.aero.CLiso(k-1);
    end
    err = 1;
    while err > s.errmax
        
        % Update propulsive efficiency (use delta from previous iteration)
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
        
        % Compute thrust
        [P_out,xi_out,phi_out,Phi_out,~] = ...
                PowerTransmissionComputation_v2(p.config,...
                                        etas,phi,Phi,xi,[],[],[],Pa);
        T = P_out.p/v;
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

    % Compute excess power
    SEP = ((1-chi*(1-cosd(p.(con).Gamma)))*T - D)*v/W;
        
    % Distinguish cases: cruise altitude reached
    if h >= h_target
        x = 0;
      
    % Cruise velocity reached
    elseif v >= M_target*f.a(h_target) && v/f.a(h) >= M_target;
        x = 1;

    % Any other point: follow x-profile    
    else
        x = interp1(h_array,x_array,h,'linear');
        if isnan(x); x = 0.5; end
    end
    
    % Calculate climb rate and horizontal acceleration
    dhdt = SEP*x;
    dVdt = SEP*(1-x)*g/v;
  
    % Deltas corresponding to this timestep
    deltah = dhdt*dt;
    deltav = dVdt*dt;
    deltaR = v*dt*cos(G);
    deltaEbat = -P_out.bat*dt;% Energy consumed: negative if discharging
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e') ||...
            strcmp(p.config,'e-2')
        deltaEf = 0;
    else
        deltaEf = -P_out.f*dt;
    end
    deltaW = g*deltaEf/p.SE.f;
    if deltaW > 0; error('Fuel mass has increased'); end
    
    % Update data for next iteration
    A.t(k+1) = A.t(k) + dt;
    A.h(k+1) = min([h_target A.h(k) + deltah]);
    A.W(k+1) = A.W(k) + deltaW;
    if (A.v(k) + deltav)/f.a(A.h(k+1)) > M_target
        A.v(k+1) = M_target*f.a(A.h(k+1));
    else
        A.v(k+1) = A.v(k) + deltav;
    end
    A.G(k+1) = dhdt/v;
    A.dVdt(k+1) = dVdt;
    A.dhdt(k+1) = dhdt;
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
    
    % Save power-related variables
    A.P.x(k) = x;
    if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e') ||...
            strcmp(p.config,'e-2')
        A.P.availableGT(k) = NaN;
    else
        A.P.availableGT(k) = Pa;
    end
    A.P.drag(k) = D*v;
    A.P.acceleration(k) = W*dVdt*v/g;
    A.P.climb(k) = W*dhdt;
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
        if strcmp(con,'Dcl')
            error('Could not reach diversion altitude/Mach number')
        else
            error('Could not reach cruise altitude/Mach number')
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
