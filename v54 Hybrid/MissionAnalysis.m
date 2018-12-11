%%% Description
%
% This routine calculates the weight breakdown of the aircraft by analyzing
% a sizing mission using the wing- and power-loading values obtained from
% the WS-WP diagram. The mission is assumed to be comprised of climb +
% cruise + descent + climb + diversion + descent.
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 29-03-17
%%% Last modified: 16-04-18


%% Initialize

% Copy climb characteristics to descent characteristics. Descent phase not
% used in WP-Ws diagram, only in MA. Currently "hard coded" here to avoid
% excessive input variables, but can be adapted here if the characteristics
% of the aircraft/mission have to be modified. Also copy parameters of
% nominal mission to diversion mission.
a.de = a.cl;
p.de = p.cl;
a.Dcr = a.cr; p.Dcr = p.cr; 
a.Dcl = a.cl; p.Dcl = p.cl; 
a.Dde = a.de; p.Dde = p.de; 
            

%% Start loop on MTOM

% Initial guess of take-off mass [kg]. 
if electric == 0
    TOM0 = (MA_in.PL+MA_in.OEM)/(1-FF_tot0*(1+p.SE.f/p.SE.bat*DOH_tot0/(1-DOH_tot0)));
else
    TOM0 = (MA_in.PL+MA_in.OEM)/(1-BF_tot0);
end

% Initialize loop variables
count = 0;
err = 1;
conv.DOH = DOH_tot0;
if electric == 0; conv.FF = FF_tot0; else conv.FF = BF_tot0; end;
conv.TOM = TOM0;
conv.err = err;
    
% Loop till TOM convergence
while err > s.errmax
    
    % Progress
    count = count + 1;
    disp([s.levelString '  > Iteration ' num2str(count)])
    
    % Installed powers [W] 
    aircraft.Pdes = structfun(@(x) TOM0*c.g./x,WPdes.(s.SelDes),'UniformOutput',0); 
    
    % DP Disk area (D^2) [m2]
    aircraft.D2W = f.D2W(WSdes.(s.SelDes),p.b_dp,p.N,p.dy,a.AR);
    aircraft.D2 = c.g*TOM0*aircraft.D2W;
    
    % Wing area [m2]
    aircraft.Sw = c.g*TOM0/WSdes.(s.SelDes);
    
    % Prop radius/wing chord ratio
    aircraft.Rc = 0.5*(aircraft.D2W*WSdes.(s.SelDes)*a.AR)^0.5;
    
    % Save aircraft weight
    aircraft.TOM = TOM0;
    
    % Compute battery energy and fuel energy status at end of nominal
    % mission
    % During first iteration, use initial guesses
    if count == 1
        
        % Fuel mass used during nominal mission and diversion
        Mf_miss0 = TOM0*FF_miss0;
        Mf_tot0 = TOM0*FF_tot0;
        Mf_div = Mf_tot0-Mf_miss0;
        
        % Battery enrgy used during nominal mission and diversion
        if electric == 0
            Ebat_miss0 = Mf_miss0*DOH_miss0/(1-DOH_miss0)*p.SE.f;
            Ebat_tot0 = Mf_tot0*DOH_tot0/(1-DOH_tot0)*p.SE.f;
        else
            Ebat_miss0 = TOM0*BF_miss0*p.SE.bat;
            Ebat_tot0 = TOM0*BF_tot0*p.SE.bat;
        end
        Ebat_div = Ebat_tot0-Ebat_miss0;
        
        % Initial guess for energy remaining at end-of-mission: no battery
        % energy left
        Ebat_end = 0;
    
    % During remaining iterations, use info of previous iteration
    else
    
        % Fuel mass used during nominal mission and diversion
        Mf_miss0 = (MA.cl.Ef(1)-MA.de.Ef(end))/p.SE.f;
        Mf_tot0 = (MA.cl.Ef(1)-MA.Dde.Ef(end))/p.SE.f;
        Mf_div = Mf_tot0-Mf_miss0;
        
        % Battery energy used during nominal mission and diversion
        Ebat_miss0 = MA.cl.Ebat(1) - MA.de.Ebat(end);
        Ebat_tot0 = MA.cl.Ebat(1) - MA.Dde.Ebat(end);
        Ebat_div = Ebat_tot0-Ebat_miss0;
        Ebat_end = MA.Dde.Ebat(end);
    end  
    clear('MA')
    
    
    %% 1. Climb segment
    % During climb the output power is specified, while the flight path is not
    
    % Initial conditions
    MA_in.cl.t(1) = 0;                          % Start of time vector
    MA_in.cl.R(1) = 0;                          % No range flown yet
    MA_in.cl.h(1) = m.TO.h;                     % Altitude of TO runway
    MA_in.cl.W(1) = TOM0*c.g;                   % Start at TOM
    MA_in.cl.v(1) = AEROdes.(s.SelDes).TO.v;    % Take V2 speed as initial speed
    MA_in.cl.G(1) = 0.01;                       % Give generic initial climb gradient
    MA_in.cl.dVdt(1) = 0.1;                     % Give generic initial acceleration
    MA_in.cl.Ebat(1) = Ebat_tot0;               % All battery energy left
    MA_in.cl.Ef(1) = Mf_tot0*p.SE.f;            % All fuel left

    % Call function
    [MA.cl] = MissionAnalysisClimb(a,p,f,s,c,'cl',aircraft,MA_in.cl,...
        m.cr.h,m.cr.M,s.dt.cl);
    
    
    %% 2. Descent segment
    % For the descent segment, the program iterates backwards, starting at
    % landing altitude, with zero fuel and zero battery energy left
    
    % Initial conditions: first index corresponds to last timestep. Output
    % arrays are flipped such that the first index corresponds to first
    % timestep.
    MA_in.de.t(1) = 0;                          % Start of time vector
    MA_in.de.R(1) = MA_in.R;                    % End of nominal mission
    MA_in.de.h(1) = m.L.h;                      % Altitude of landing runway
    MA_in.de.W(1) = (TOM0-Mf_miss0)*c.g;        % Only diversion fuel weight left at landing
    MA_in.de.v(1) = m.L.vs*1.3;                 % Approach speed
    MA_in.de.G(1) = 0.0;                        % Generic final climb gradient
    MA_in.de.dVdt(1) = 0.0;                     % Generic final deceleration
    MA_in.de.Ebat(1) = Ebat_div;                % Only diversion battery energy left
    MA_in.de.Ef(1) = Mf_div*p.SE.f;             % Only diversion fuel energy left
    
    % Call function
    [MA.de] = MissionAnalysisDescent(a,p,f,s,c,'de',aircraft,MA_in.de,...
        m.cr.h,m.cr.M,s.dt.de);

                         
    %% 3. Cruise segment
    % During cruise the flight path is specified, while one of the power control
    % parameters is not
    
    % Range flown during cruise segment (depends on ranges covered in climb and
    % descent phases) [m]
    R_climb = max(MA.cl.R);
    R_descent = max(MA.de.R)-min(MA.de.R);
    R_cruise = MA_in.R - R_climb - R_descent;
    if R_cruise < 0
        error(['Climb and descent require more than the '...
                'specified nominal mission range'])
    end
    
    % Initial conditions
    MA_in.cr.t(1) = MA.cl.t(end);
    MA_in.cr.R(1) = MA.cl.R(end);
    MA_in.cr.h(1) = m.cr.h;
    MA_in.cr.W(1) = MA.cl.W(end);
    MA_in.cr.v(1) = m.cr.M*f.a(m.cr.h);
    MA_in.cr.G(1) = 0;
    MA_in.cr.dVdt(1) = 0;
    MA_in.cr.Ebat(1) = MA.cl.Ebat(end);
    MA_in.cr.Ef(1) = MA.cl.Ef(end);
    
    % Call function
    [MA.cr] = MissionAnalysisCruise(a,p,f,s,c,'cr',aircraft,...
        MA_in.cr,R_cruise,s.dt.cr);

    % Shift time-axis of descent phase
    MA.de.t = MA.de.t + abs(min(MA.de.t)) + max(MA.cr.t);
    
      
    %% 4. Diversion climb segment
    % During climb the output power is specified, while the flight path is not
    
    % Initial conditions
    MA_in.Dcl.t(1) = max(MA.de.t);              % Start of time vector
    MA_in.Dcl.R(1) = MA_in.R;                   % Nominal mission range flown
    MA_in.Dcl.h(1) = m.L.h;                     % Altitude of L runway
    MA_in.Dcl.W(1) = (TOM0-Mf_miss0)*c.g;       % Only diversion fuel weight left
    MA_in.Dcl.v(1) = AEROdes.(s.SelDes).bL.v;   % Balked landing speed
    MA_in.Dcl.G(1) = 0.01;                      % Give generic initial climb gradient
    MA_in.Dcl.dVdt(1) = 0.1;                    % Give generic initial acceleration
    MA_in.Dcl.Ebat(1) = Ebat_div;               % Only diversion battery energy left
    MA_in.Dcl.Ef(1) = Mf_div*p.SE.f;            % Diversion fuel left

    % Call function
    [MA.Dcl] = MissionAnalysisClimb(a,p,f,s,c,'Dcl',aircraft,MA_in.Dcl,...
        MA_in.h_div,MA_in.M_div,s.dt.Dcl);
    
    
    %% 5. Diversion descent segment
    % For the descent segment, the program iterates backwards, starting at
    % landing altitude, with zero fuel and zero battery energy left
    
    % Initial conditions: first index corresponds to last timestep. Output
    % arrays are flipped such that the first index corresponds to first
    % timestep.
    MA_in.Dde.t(1) = 0;                         % Start of time vector
    MA_in.Dde.R(1) = MA_in.R+MA_in.R_div;       % End-of-mission
    MA_in.Dde.h(1) = m.L.h;                     % Altitude of landing runway
    MA_in.Dde.W(1) = (TOM0-Mf_tot0)*c.g;        % No fuel weight left at landing
    MA_in.Dde.v(1) = m.L.vs*1.3;                % Approach speed
    MA_in.Dde.G(1) = 0.0;                       % Generic final climb gradient
    MA_in.Dde.dVdt(1) = 0.0;                    % Generic final deceleration
    MA_in.Dde.Ebat(1) = Ebat_end;               % Some battery energy left
    MA_in.Dde.Ef(1) = 0;                        % No fuel left
    
    % Call function
    [MA.Dde] = MissionAnalysisDescent(a,p,f,s,c,'Dde',aircraft,MA_in.Dde,...
        MA_in.h_div,MA_in.M_div,s.dt.Dde);
    
    
    %% 6. Diversion cruise segment
    % During cruise the flight path is specified, while one of the power control
    % parameters is not
    
    % Range flown during cruise segment (depends on ranges covered in climb and
    % descent phases) [m]
    R_Dclimb = max(MA.Dcl.R)-min(MA.Dcl.R);
    R_Ddescent = max(MA.Dde.R)-min(MA.Dde.R);
    R_Dcruise = MA_in.R_div - R_Dclimb - R_Ddescent;
    if R_Dcruise < 0
        error(['Diversion climb and descent require more than the '...
                'specified diversion range'])
    end
    
    % Initial conditions
    MA_in.Dcr.t(1) = MA.Dcl.t(end);
    MA_in.Dcr.R(1) = MA.Dcl.R(end);
    MA_in.Dcr.h(1) = MA_in.h_div;
    MA_in.Dcr.W(1) = MA.Dcl.W(end);
    MA_in.Dcr.v(1) = MA_in.M_div*f.a(MA_in.h_div);
    MA_in.Dcr.G(1) = 0;
    MA_in.Dcr.dVdt(1) = 0;
    MA_in.Dcr.Ebat(1) = MA.Dcl.Ebat(end);
    MA_in.Dcr.Ef(1) = MA.Dcl.Ef(end);
    
    % Call function
    [MA.Dcr] = MissionAnalysisCruise(a,p,f,s,c,'cr',aircraft,...
        MA_in.Dcr,R_Dcruise,s.dt.Dcr);

    % Shift time-axis of descent phase
    MA.Dde.t = MA.Dde.t + abs(min(MA.Dde.t)) + max(MA.Dcr.t);
    
    
    %% Evaluate convergence 

    % Recompute TOM
    [M,MA,DOHmiss,DOHtot] = ComputeWeights(a,m,p,f,s,c,aircraft,MA_in.PL,MA);
    TOM1 = M.TOM;
    Mf_tot1 = M.f;
    Mf_miss1 = M.f_miss;
    Ef_tot1 = Mf_tot1*p.SE.f;
    Ebat_tot1 = M.bat*p.SE.bat;
    
    % Update DOH and FF (swap NaN with 0 for configurations which do not
    % use batteries, to get a coherent error value). Also update nominal
    % mission DOH and FF, not to use in convergence, but to store. This
    % will allow the user to get a more accurate initial guess for the next
    % time the code is run
    if electric == 0
        if strcmp(p.config,'conventional') || ...
                strcmp(p.config,'turboelectric') || strcmp(p.config,'PTE')
            DOH_tot1 = 0;
        else
            DOH_tot1 = DOHtot;
        end
        FF_tot1 = Mf_tot1/TOM1;
    else
        DOH_tot1 = 1;
        BF_tot1 = Mbat_tot1/TOM1;
    end
    
    % Calculate error: if denominator is zero, use 1 (since FF and DOH are
    % normally of order 1 anyway) to avoid a 0/0 indeterminate
    if DOH_tot0 == 0; DOH_totref = 0.1; else DOH_totref = DOH_tot0; end;
    if electric == 0
        err = 1/3*(abs((TOM1-TOM0)/TOM0) + abs((FF_tot1-FF_tot0)/...
            FF_tot0) + abs((DOH_tot1-DOH_tot0)/DOH_totref));
    else
        err = 1/3*(abs((TOM1-TOM0)/TOM0) + abs((BF_tot1-BF_tot0)/...
            BF_tot0) + abs((DOH_tot1-DOH_tot0)/DOH_totref));
    end
    
    % Break if not converging
    if count > s.itermax;
        error('Mission analysis did not converge');
    end    
    
    % Save convergence of variables
    conv.DOH = [conv.DOH DOH_tot1];
    if electric == 0
        conv.FF = [conv.FF FF_tot1];
    else
        conv.FF = [conv.FF BF_tot1];
    end
    conv.TOM = [conv.TOM TOM1];
    conv.err = [conv.err err];
        
    % Update
    TOM0 = TOM1;
    if electric == 0
        FF_tot0 = FF_tot1;
        DOH_tot0 = DOH_tot1;
    else
        BF_tot0 = BF_tot1;
    end
end


%% Check power envelope

% Powertrain components/paths to be verified (fields of WPdes.(SelDes)
% structure)
Pnames = {'GTM','EM1','EM2','bat','s1'};

% String for warning
Pstrings = {'gas turbine','primary electrical machine',...
            'secondary electrical machine','battery','primary shaft'};

% Loop over mission segments
names = fieldnames(MA);
Ncon = size(names,1);
for i = 1:Ncon
    
    % Save an additional field in the MA.() structure which indicates
    % whether any power limit has been exceeded (for plotting)
    MA.(names{i}).limits = zeros(length(MA.(names{i}).t),1);
    
    % Evaluate limit of each component
    for j = 1:size(Pnames,2)
       
        % Generate a field in the MA.().P.Limit structure, containing an
        % array with value 1 for timesteps where the maximum allowable
        % power of the given component has be surpassed, and zeros
        % elsewhere.
        MA.(names{i}).P.limits.(Pnames{j}) = ...
                zeros(length(MA.(names{i}).t),1);
        
        % Loop over time steps of current segment
        for k = 1:length(MA.(names{i}).t)
            
            % Define limting power depending on the powertrain
            % component/path considered
            switch Pnames{j}
                
                % For the gas turbine, compute current available power
                case 'GTM'
                    rho = f.rho(MA.(names{i}).h(k));
                    Pavail = aircraft.Pdes.GTM*f.Alpha(rho);
                    Preq = MA.(names{i}).P.gt(k);

                % For primary electric machine, select sizing power
                case 'EM1'
                    Pavail = aircraft.Pdes.EM1;
                    Preq = max([MA.(names{i}).P.gb(k)...
                        MA.(names{i}).P.e1(k)]);
                    
                % For primary secondary machine, select sizing power    
                case 'EM2'
                    Pavail = aircraft.Pdes.EM2;
                    Preq = max([MA.(names{i}).P.s2(k)...
                        MA.(names{i}).P.e2(k)]);
                    
                % Battery power
                case 'bat'
                    Pavail = aircraft.Pdes.bat;
                    Preq =MA.(names{i}).P.bat(k);
                    
                % Primary shaft power    
                case 's1'
                    Pavail = aircraft.Pdes.s1;
                    Preq =MA.(names{i}).P.s1(k);
            end
            
            % If limits have been surpassed, save.
            if Preq > Pavail
                MA.(names{i}).P.limits.(Pnames{j})(k) = 1;
            end
        end
        
        % Select mission segment string
        switch names{i}
            case 'cl'; Sstring = 'climb';
            case 'cr'; Sstring = 'cruise';
            case 'de'; Sstring = 'descent';
            case 'Dcl'; Sstring = 'diversion climb';
            case 'Dcr'; Sstring = 'diversion cruise';
            case 'Dde'; Sstring = 'diversion descent';
        end 
        
        % Issue warning if any power has been surpassed
        if any(MA.(names{i}).P.limits.(Pnames{j}))==1
            disp([s.levelString '  > Warning: maximum installed ' ...
                  Pstrings{j} ' power has been exceeded during the '...
                  Sstring ' segment! Check power-control envelope.'])
        end
        
        % If the current component power is exceeded, change corresponding
        % index in overall mission-segment limit array
        MA.(names{i}).limits(MA.(names{i}).P.limits.(Pnames{j})==1) = 1;
    end
end


%% Organize output

% Save variables to structures. Note that DOH_tot will now include any 
% possible battery energy not used during the mission, but available due to
% power requirements
if electric == 1
    MA.BF_tot = BF_tot1;
    MA.BF_miss = M.bat_miss/M.TOM;
else
MA.FF_tot = FF_tot1;
MA.FF_miss = M.f_miss/M.TOM;
end
MA.DOH_tot = M.bat*p.SE.bat/(M.f*p.SE.f + M.bat*p.SE.bat);
MA.DOH_miss = DOHmiss;
MA.conv = conv;

% Clear loop variables (can be commented for debugging)
clear('FF_tot0','FF_miss0','DOH_tot0','DOH_miss0','DOHmiss','DOH_miss',...
      'FF_tot1','FF_miss1','DOH_tot1','DOH_miss1','DOHtot',...
      'electric','count','err','TOM0','TOM1','i','conv',...
      'Ebat_miss1','Ebat_tot1','EbatMiss','EbatTot',...
      'Ebat_div','Ebat_end','Ebat_miss0','Ebat_tot0',...
      'Ef_miss1','Ef_tot1','Mbat_div','Mbat_miss0','Mbat_miss1',...
      'Mbat_tot0','Mbat_tot1','Mf_div','Mf_miss0','Mf_miss1',...
      'Mf_tot0','Mf_tot1','names','Ncon','R_climb','R_cruise',...
      'R_Dclimb','R_Dcruise','R_descent','R_Ddescent','j','k',...
      'DOH_missref','DOH_totref','Pavail','Pnames','Preq','Pstrings',...
      'rho','Sstring');

% Call plotting script
if s.plotMA == 1;
    [s] = PlotMissionAnalysis(p,m,s,c,aircraft,M,MA);
end
clear('aircraft')

