function [WS,WP_path,WP_comp,WP_loss,TW,a,m,p] = ConstraintTakeOffDistance(a,m,p,f,s,c)
%% Take-off operating conditions
% Note: this function assumes both propulsive systems are generating thrust
% in TO conditions (i.e. Phi in [0,1]).
% Note 2: velocity margin changed in Initiatior!

% Density 
m.TO.rho = f.rho(m.TO.h);

% The TOP relation is based on shaft power ratio. 
% Compute weighted average of propulsive efficiency using shaft power ratio
% if applicable, if not select corresponding powertain (primary/secondary).
% If Phi is applicable but unknown, estimate eta_p by taking the average.
% This value will be corrected by adding detap in proportion to Phi or
% 1-Phi, depending on p.DP.
% Note that "etap" here refers to the overall propulsive efficiency, while
% p.TO.etap refers to the propulsive efficiency of the DP propulsors.
if ~isnan(m.TO.Phi)
    etap = (1-m.TO.Phi) * p.TO.etap1 + (m.TO.Phi) * p.TO.etap2;
else
    if strcmp(p.config,'PTE') || strcmp(p.config,'SPPH') || ...
       strcmp(p.config,'dual-e')
           etap = 0.5*(p.TO.etap1+p.TO.etap1);
           disp([s.levelString '  > Phi is unknown at this stage, so '...
                 'propulsive efficiency in TO conditions is estimated '...
                 'as the average of the primary and secondary '...
                 'propulsors'' efficiencies'])
    elseif strcmp(p.config,'conventional') || ...
           strcmp(p.config,'e-1') || strcmp(p.config,'parallel')
        etap = p.TO.etap1;
    elseif strcmp(p.config,'turboelectric') || ...
           strcmp(p.config,'e-2') || strcmp(p.config,'serial') 
        etap = p.TO.etap2;
    end
end


%% Compute power loading (not guaranteeing equilibrium flight)

% Initialize loop variables
WS = linspace(0,s.WSmax,s.n); 
WP = NaN(1,s.n);                           
TW = NaN(1,s.n);    
m.TO.v = NaN(1,s.n);         

% Loop over all wing loading values
for i = 1:s.n
    
    % Reset convergence loop values
    err = 1;               
    iter = 0;           
    
    % Initial guess for stall power loading in TO config
    WP0s = 0.2;                                         
    
    % Initial guess for actual power loading in TO config
    WS0 = WS(i);
    
    % Inital guess for stall speed in TO config (at screen height)
    vs = (WS0/a.TO.CLmax/0.5/m.TO.rho)^0.5;         
    
    % Initial guess for delta-eta_p
    detap = 0;
    
    % Rotor sizing
    D2W_TO = f.D2W(WS0,p.b_dp,p.N,p.dy,a.AR);           

    % Compute stall speed for each wing loading
    while err>s.errmax
        
        % Loop counter
        iter = iter+1;
        
        % Update effective propulsive efficiency (see comment in first
        % section regarding what happens if Phi is unknown: generic share 
        % of 0.5 is used)
        if ~isnan(m.TO.Phi); Phi = m.TO.Phi; else Phi = 0.5; end
        if p.DP == 1;       eta_eff = etap+detap*(1-Phi);
        elseif p.DP == 2;   eta_eff = etap+detap*Phi;
        else                eta_eff = etap; 
        end
        
        % Thrust loading (in TO conditions at stall speed, not at TO-speed)
        TW0s = eta_eff/WP0s/vs;                       
        
        % Compute thrust coefficient, defined as Tc = T/(rho*v^2*D^2)
        Tcs = p.TO.T*TW0s/p.N/m.TO.rho/vs^2/D2W_TO;
        
        % Compute prop radius/wing chord ratio
        Rc = 0.5*(D2W_TO*WS0*a.AR)^0.5;
        
        % Estimate Mach number
        M = vs/f.a(m.TO.h);
        
        % Compute delta terms in stall conditions
        [dCl,~,~,detap] = WingPropDeltas(Tcs,Rc,p.xp,a.AR,a.TO.CLmax,...
                                     M,a.Lambda,p.TO.Gamma,p.TO.etap,0);
                   
        % Effective wing cl increase due to DP TV in stall conditions
        dCL_TV = p.TO.T*TW0s*WS0*sind(p.TO.Gamma)/0.5/m.TO.rho/vs^2;
        
        % Total maximum lift coefficient
        CL_tot = a.TO.CLmax + p.b_dp*dCl + dCL_TV;      
        
        % Updated stall velocity
        vs = (WS0/CL_tot/0.5/m.TO.rho)^0.5;         
        
        % Shaft power loading at stall speed
        WP1s = f.TOP(m.TO.s)*CL_tot*(m.TO.rho/c.rho_SL)/WS0...
                    *(1 - p.TO.T + p.TO.T*cosd(p.TO.Gamma)); 
        
        % Calculate error & update values
        if isnan(WP1s)
            err = 1;
            
            % Add small offset to avoid reiterating NaN values
            WP0s = 0.21;                                
        else
            err = abs(WP1s-WP0s)/abs(WP0s); 
            WP0s = WP1s;
        end    
        
        % Stop at max iterations
        if iter>=s.itermax
            CL_tot = NaN;
            break
        end
    end
    
    % Shaft power loading at take-off speed
    % Note: Not sure of the rho/rho_SL factor should be included when
    %   calculating shaft power (gas turbine altitude lapse should not
    %   be included). For TO at SL the results is the same, but caution
    %   should be exerted if TO/L is carried out at higher altitudes.
    WP1 = f.TOP(m.TO.s)*CL_tot/1.1^2*(m.TO.rho/c.rho_SL)/WS0...
        *(1 - p.TO.T + p.TO.T*cosd(p.TO.Gamma));
    
    % Convert to required power
    WP(i) = WP1/eta_eff;

    % Store speed, lift coefficient, and propulsive efficiency change
    m.TO.v(i) = vs*1.1;
    m.TO.M(i) = vs*1.1/f.a(m.TO.h);
    a.TO.CL(i) = CL_tot/1.1^2;
    p.TO.detap(i) = detap;
    
    % Save thrust loading. Note: may not be very representative!
    TW(i) = 1/WP(i)/m.TO.v(i);    
    
end

      
%% Compute HEP component power loading

% All engines operative
OEI = 0;

% Replace NaNs with Inf's so that the code understands that WP is an input.
% Keep indices to revert changes later
indices = isnan(WP);
WP(indices) = Inf;

% Call sizing routine
[WP_path,WP_comp,WP_loss] = SizePowertrain(WP,'TO',OEI,m,p,f,s);

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












