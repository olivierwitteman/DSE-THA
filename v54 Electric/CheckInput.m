%%% Description
% This script checks the parameters specified in Input.m. If it
% automatically corrects a parameter (e.g. a value has been specified but
% it should be a NaN) it shows a message in the command window. If it
% cannot correct a parameter (e.g. a NaN has been given, but a
% specific value is required), it returns an error. Only some of the
% powertrain-related parameters are checked. Additionally, it selects which
% of the two powertrain branches' (primary or secondary) properties should
% be applied to the DP system.
%
% 
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 21-06-17
%%% Last modified: 11-04-18


%% Check powertrain input settings compatibility

% Identify constraints that are being used
phases = fieldnames(m);
np = length(phases);

% Check control laws and DP selection are correctly specified
switch p.config
    
    % 1. Conventional powertrain
    case 'conventional'
        
        % Only throttle can be defined by user
        
        % Only primary propulsors exist, overwrite shaft power ratio 
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % No batteries exist, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on primary powertrain
        if p.DP == 2
            p.DP = 1;
            disp([s.levelString '  > Changed DP location to primary '...
                  'powertrain (' p.config ' configuration).'])
        end

    % 2. Turboelectric powertrain
    case 'turboelectric'
        
        % Only throttle can be defined by user
        
        % Only secondary propulsors exist, overwrite shaft power ratio 
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % No batteries exist, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on secondary powertrain
        if p.DP == 1
            p.DP = 2;
            disp([s.levelString '  > Changed DP location to secondary '...
                  'powertrain (' p.config ' configuration).'])
        end

    % 3. Serial powertrain
    case 'serial'
        
        % Throttle and supplied power ratio can be defined by user
        
        % Only secondary propulsors exist, overwrite shaft power ratio
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on secondary powertrain
        if p.DP == 1
            p.DP = 2;
            disp([s.levelString '  > Changed DP location to secondary '...
                  'powertrain (' p.config ' configuration).'])
        end
        
    % 4. Parallel powertrain
    case 'parallel'
        
        % Throttle and supplied power ratio can be defined by user
        
        % Only primary propulsors exist, overwrite shaft power ratio
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on primary powertrain
        if p.DP == 2
            p.DP = 1;
            disp([s.levelString '  > Changed DP location to primary '...
                  'powertrain (' p.config ' configuration).'])
        end
           
    % 5. Partial turboelectric powertrain
    case 'PTE'
        
        % Throttle and shaft power ratio can be defined by user
                
        % No batteries exist, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
    % 6. S-P partial hybrid powertrain
    case 'SPPH'
        
        % Throttle, shaft power ratio and supplied powe ratio can be
        % specified by user. DP can be on both powertrains.
       
    % 7. Full-electric (e-1) powertrain
    case 'e-1'
        
        % Only throttle can be defined by user, which refers to
        % electromotor throttle!
        
        % Only primary propulsors exist, overwrite shaft power ratio
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % No fuel exists, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on primary powertrain
        if p.DP == 2
            p.DP = 1;
            disp([s.levelString '  > Changed DP location to primary '...
                  'powertrain (' p.config ' configuration).'])
        end
        
    % 8. Full-electric (e-2) powertrain
    case 'e-2'
        
        % Only throttle can be defined by user, which refers to
        % electromotor throttle!
        
        % Only secondary propulsors exist, overwrite shaft power ratio
        if any(structfun(@(p) ~isnan(p.Phi), m))
            for i = 1:np; m.(phases{i}).Phi = NaN; end
            disp([s.levelString '  > Phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % No fuel exists, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end
        
        % If DP is used, it must be on secondary powertrain
        if p.DP == 1
            p.DP = 2;
            disp([s.levelString '  > Changed DP location to secondary '...
                  'powertrain (' p.config ' configuration).'])
        end
        
    % 9. Full-electric (dual-e) powertrain
    case 'dual-e'
        
        % Electromotor throttle and shaft power ratio can be defined by
        % user
                
        % No fuel exists, overwrite supplied power ratio 
        if any(structfun(@(p) ~isnan(p.phi), m))
            for i = 1:np; m.(phases{i}).phi = NaN; end
            disp([s.levelString '  > phi is not applicable to a '...
                  p.config ' powertrain. Changed value to NaN.'])
        end   
        
    % Give error for incorrectly specified powertrian configurations    
    otherwise
        error('Incorrect powertrain configuration specified.')
end


%% Assign parameters used for DP-evaluation to corresponding powertrain.   

% Powertrains featuring both primary and secondary propulsors   
if strcmp(p.config,'PTE') || strcmp(p.config,'SPPH') || ...
   strcmp(p.config,'dual-e')

    % If primary powertrain has DP
    if p.DP == 1                
        p.N = p.N1;             % Number of propulsors in DP system
                
        % Fraction of total required thrust produced by DP system, and
        % propulsive efficiency of DP propulsors 
        for i=1:np
            p.(phases{i}).T = 1/( m.(phases{i}).Phi / (1-m.(phases{i}).Phi) ...
                * p.(phases{i}).etap2/p.(phases{i}).etap1 + 1);
            p.(phases{i}).etap = p.(phases{i}).etap1;
        end
        
    % If secondary powertrain has DP    
    elseif p.DP == 2            
        p.N = p.N2;
        
        for i=1:np
            p.(phases{i}).T = 1/( (1-m.(phases{i}).Phi) / m.(phases{i}).Phi ...
                * p.(phases{i}).etap1/p.(phases{i}).etap2 + 1);
            p.(phases{i}).etap = p.(phases{i}).etap2;
        end
    
    % If none of the powertrains have DP
    elseif p.DP == 0;           
        p.N = NaN;
        for i=1:np
            p.(phases{i}).T = 0;
            p.(phases{i}).etap = (1-m.(phases{i}).Phi)* ...
                p.(phases{i}).etap1 + (m.(phases{i}).Phi)*...
                p.(phases{i}).etap2;
        end
    
    % Wrong configuration specified
    else
        error('Incorrect DP-powertrain (primary/secondary/none) assigned')
    end
    
% Powertrains featuring only primary propulsors    
elseif strcmp(p.config,'conventional') || strcmp(p.config,'parallel') ||...
        strcmp(p.config,'e-1')
    p.N = p.N1;
    for i=1:np;
        p.(phases{i}).T = 1;
        p.(phases{i}).etap = p.(phases{i}).etap1;
    end
        
% Powertrains featuring only secondary propulsors           
elseif strcmp(p.config,'turboelectric') || strcmp(p.config,'serial') || ...
        strcmp(p.config,'e-2')
    p.N = p.N2;
    for i=1:np;
        p.(phases{i}).T = 1;
        p.(phases{i}).etap = p.(phases{i}).etap2;
    end
end

% If DP==0, set xp = Inf, such that WingPropDeltas does not affect
% anything.
if p.DP == 0 && p.xp ~= Inf
    p.xp = Inf;
    disp([s.levelString '  > Since no DP powertrain has been specified '...
                        '(p.DP = 0), p.xp has been set to "Inf" to '...
                        ' neglect the aero-propulsive interaction '...
                        'effects.'])
end


%% Check initial conditions for MA

% Make variable names shorter (will be cleared at end of program)
FF_tot0 = MA_in.FF_tot0;
FF_miss0 = MA_in.FF_miss0;
DOH_tot0 = MA_in.DOH_tot0;
DOH_miss0 = MA_in.DOH_miss0;

% Make sure initial point is valid
if FF_tot0*(1+p.SE.f/p.SE.bat*DOH_tot0/(1-DOH_tot0)) > 1
    error(['Invalid initial guess of FF and DOH: this combination '...
           'leads to a total battery and fuel mass greater than MTOM'])
end

% Mission segments 
names = {'cl','cr','de','Dcl','Dcr','Dde'};
N = size(names,2);

% Check if initial guess for DOH is coherent with powertrain selected
if strcmp(p.config,'e-1') || strcmp(p.config,'e-2') ...
                          || strcmp(p.config,'dual-e')
    electric= 1;
    disp([s.levelString '  > Full-electric configuration. Using fuel mass'...
                ' fraction input as initial guess for battery mass'...
                ' fraction'])
            
    % For full-electric configurations (DOH_tot = 1), ignore 
    % "fuel fraction" and converge on a battery fraction instead.
    if DOH_tot0 ~= 1
        disp([s.levelString '  > Setting DOH to 1 for '...
                'fully electric configurations'])
            DOH_tot0 = 1;
            DOH_miss = 1;
    end
    
    % Define a "battery fraction" instead
    BF_tot0 = FF_tot0;
    BF_miss0 = FF_miss0;
    FF_tot0 = 0;
    FF_miss0 = 0;
    
%     % Make sure phi profile is always 1
%     for i = 1:size(names,1)
%         if prod(MA_in.(names{i}).phi==1)==0
%             disp([s.levelString '  > Setting phi-profile to NaN for '...
%                 'fully electric configurations (' names{i} ')'])
%             MA_in.(names{i}).phi = [NaN NaN];
%         end
%     end

% For gas-turbine only configurations    
elseif strcmp(p.config,'conventional')
    electric = 0;
    if DOH_tot0 ~= 1
        disp([s.levelString '  > Setting DOH to 0 for '...
                'conventional configurations'])
            DOH_tot0 = 0;
            DOH_miss = 0;
    end
    
%     % Make sure phi profile is always 0 (NaN)
%     for i = 1:size(names,1)
%         if prod(MA_in.(names{i}).phi==0)==0
%             disp([s.levelString '  > Setting phi-profile to NaN for '...
%                 'conventional configurations (' names{i} ')'])
%             MA_in.(names{i}).phi = [NaN NaN];
%         end
%     end
    
% For hybrid configurations    
else
    electric = 0;
end

% Make sure that some variables are within bounds
vars = [FF_tot0 FF_miss0 DOH_tot0 DOH_miss0 p.minSOC_miss p.minSOC_tot ...
        MA_in.cl.xi MA_in.cr.xi MA_in.de.xi ...
        MA_in.Dcl.xi MA_in.Dcr.xi MA_in.Dde.xi];
if any(vars>1) || any(vars<0)
    error(['Fuel fraction, DOH, SOC, and throttle must be between 0 '...
           'and 1. Check input settings'])
end


%% MA: Check that power-control profiles are consistent with configuration

% Power control names
Pnames = {'xi','phi','Phi'};

% Mission segments where all parameters have to be specified (climb and
% descent)
segments1 = {'cl','de','Dcl','Dde'};
N1 = size(segments1,2);

% Mission segments where one DOF must remain free (cruise)
segments2 =  {'cr','Dcr'};
N2 = size(segments2,2);

% First make sure that the control arrays have either two NaN values or no
% NaN values
for i = 1:N
    for j = 1:size(Pnames,2)
        if (isnan(MA_in.(names{i}).(Pnames{j})(1)) && ...
            ~isnan(MA_in.(names{i}).(Pnames{j})(2))) || ...
            (isnan(MA_in.(names{i}).(Pnames{j})(2)) && ...
            ~isnan(MA_in.(names{i}).(Pnames{j})(1)))
            error(['Either both or none of the fields in MA_in.' ...
                    names{i} '.' Pnames{j} ' must have a NaN value. '...
                    'Check mission analysis input profiles.'])
        end
    end
end

% Check input per powertrain configuration
switch p.config
    
    % Conventional: Throttle must be specified during climb/descent
    case 'conventional'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi))
                error(['Throttle must be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Throttle cannot be specified during cruise
        xiWarning = 0;
        for i = 1:N2
            if ~isnan(sum(MA_in.(segments2{i}).xi))
                xiWarning = 1;
                MA_in.(segments2{i}).xi = [NaN NaN];
            end
        end
        if xiWarning == 1
            disp([s.levelString '  > Throttle cannot be specified '...
                      'during cruise phases for a ' p.config ...
                      ' powertrain. Changed value to NaN.']) 
        end
        
        % Ignore phi or Phi if they have been specified
        phiWarning = 0;
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        
    % Turboelectric: same as conventional configuration
    case 'turboelectric'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi))
                error(['Throttle must be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Throttle cannot be specified during cruise
        xiWarning = 0;
        for i = 1:N2
            if ~isnan(sum(MA_in.(segments2{i}).xi))
                xiWarning = 1;
                MA_in.(segments2{i}).xi = [NaN NaN];
            end
        end
        if xiWarning == 1
            disp([s.levelString '  > Throttle cannot be specified '...
                      'during cruise phases for a ' p.config ...
                      ' powertrain. Changed value to NaN.']) 
        end
        
        % Ignore phi or Phi if they have been specified
        phiWarning = 0;
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        
    % Serial: throttle and phi must be specified during climb/descent 
    case 'serial'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi)) || ...
                    isnan(sum(MA_in.(segments1{i}).phi))
                error(['Throttle and phi must both be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Only one of the two can be specified during cruise
        for i = 1:N2
            if (isnan(sum(MA_in.(segments2{i}).xi)) && ...
                isnan(sum(MA_in.(segments2{i}).phi))) || ...
               (~isnan(sum(MA_in.(segments2{i}).xi)) && ...
                ~isnan(sum(MA_in.(segments2{i}).phi)))
               error(['Either throttle OR phi must be specified during '...
                     'cruise phases for a ' p.config ' powertrain'])
            end
        end
        
        % Ignore Phi if it has been specified
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        
    % Parallel: same as serial configuration
    case 'parallel'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi)) || ...
                    isnan(sum(MA_in.(segments1{i}).phi))
                error(['Throttle and phi must both be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Only one of the two can be specified during cruise
        for i = 1:N2
            if (isnan(sum(MA_in.(segments2{i}).xi)) && ...
                isnan(sum(MA_in.(segments2{i}).phi))) || ...
               (~isnan(sum(MA_in.(segments2{i}).xi)) && ...
                ~isnan(sum(MA_in.(segments2{i}).phi)))
               error(['Either throttle OR phi must be specified during '...
                     'cruise phases for a ' p.config ' powertrain'])
            end
        end
        
        % Ignore Phi if it has been specified
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                'the ' p.config ' powertrain''s MA. Changed value '...
                'to NaN.'])
        end
        
    % PTE: throttle and Phi must be specified during climb/descent
    case 'PTE'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi)) || ...
                    isnan(sum(MA_in.(segments1{i}).Phi))
                error(['Throttle and Phi must both be specified during '...
                    'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Only one of the two can be specified during cruise
        for i = 1:N2
            if (isnan(sum(MA_in.(segments2{i}).xi)) && ...
                    isnan(sum(MA_in.(segments2{i}).Phi))) || ...
                    (~isnan(sum(MA_in.(segments2{i}).xi)) && ...
                    ~isnan(sum(MA_in.(segments2{i}).Phi)))
                error(['Either throttle OR Phi must be specified during '...
                    'cruise phases for a ' p.config ' powertrain'])
            end
        end
        
        % Ignore phi if it has been specified
        phiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                'the ' p.config ' powertrain''s MA. Changed value '...
                'to NaN.'])
        end
        
    % SPPH: throttle, phi and Phi must be specified during climb/descent
    case 'SPPH'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi)) || ...
                    isnan(sum(MA_in.(segments1{i}).phi)) || ...
                    isnan(sum(MA_in.(segments1{i}).Phi))
                error(['Throttle, phi and Phi must all be specified '...
                    'during climb/descent phases for a ' p.config ...
                    ' powertrain'])
            end
        end
        
        % Only two out of three can be specified during cruise
        for i = 1:N2
            xiStatus = isnan(sum(MA_in.(segments2{i}).xi));
            phiStatus = isnan(sum(MA_in.(segments2{i}).phi));
            PhiStatus = isnan(sum(MA_in.(segments2{i}).Phi));
            if sum([xiStatus phiStatus PhiStatus]) ~= 1
                error(['Exactly two of the three power-control '...
                       'parameters (throttle, phi and Phi) must be '...
                       'specified during cruise phases for a ' p.config ...
                       ' powertrain'])
            end
        end
        
    % e-1: same as conventional, although in this case throttle is electric
    case 'e-1'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi))
                error(['Throttle must be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Throttle cannot be specified during cruise
        xiWarning = 0;
        for i = 1:N2
            if ~isnan(sum(MA_in.(segments2{i}).xi))
                xiWarning = 1;
                MA_in.(segments2{i}).xi = [NaN NaN];
            end
        end
        if xiWarning == 1
            disp([s.levelString '  > Throttle cannot be specified '...
                      'during cruise phases for a ' p.config ...
                      ' powertrain. Changed value to NaN.']) 
        end
        
        % Ignore phi or Phi if they have been specified
        phiWarning = 0;
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        
    % e-2: same as e-1
    case 'e-2'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi))
                error(['Throttle must be specified during '...
                     'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Throttle cannot be specified during cruise
        xiWarning = 0;
        for i = 1:N2
            if ~isnan(sum(MA_in.(segments2{i}).xi))
                xiWarning = 1;
                MA_in.(segments2{i}).xi = [NaN NaN];
            end
        end
        if xiWarning == 1
            disp([s.levelString '  > Throttle cannot be specified '...
                      'during cruise phases for a ' p.config ...
                      ' powertrain. Changed value to NaN.']) 
        end
        
        % Ignore phi or Phi if they have been specified
        phiWarning = 0;
        PhiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
            if ~isnan(sum(MA_in.(names{i}).Phi))
                PhiWarning = 1;
                MA_in.(names{i}).Phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        if PhiWarning == 1
            disp([s.levelString '  > Phi-profile is not applicable to '...
                  'the ' p.config ' powertrain''s MA. Changed value '...
                  'to NaN.']) 
        end
        
    % dual-e: same as PTE, but throttle refers to electrical machines
    case 'dual-e'
        for i = 1:N1
            if isnan(sum(MA_in.(segments1{i}).xi)) || ...
                    isnan(sum(MA_in.(segments1{i}).Phi))
                error(['Throttle and Phi must both be specified during '...
                    'climb/descent phases for a ' p.config ' powertrain'])
            end
        end
        
        % Only one of the two can be specified during cruise
        for i = 1:N2
            if (isnan(sum(MA_in.(segments2{i}).xi)) && ...
                    isnan(sum(MA_in.(segments2{i}).Phi))) || ...
                    (~isnan(sum(MA_in.(segments2{i}).xi)) && ...
                    ~isnan(sum(MA_in.(segments2{i}).Phi)))
                error(['Either throttle OR Phi must be specified during '...
                    'cruise phases for a ' p.config ' powertrain'])
            end
        end
        
        % Ignore phi if it has been specified
        phiWarning = 0;
        for i = 1:N
            if ~isnan(sum(MA_in.(names{i}).phi))
                phiWarning = 1;
                MA_in.(names{i}).phi = [NaN NaN];
            end
        end
        if phiWarning == 1
            disp([s.levelString '  > phi-profile is not applicable to '...
                'the ' p.config ' powertrain''s MA. Changed value '...
                'to NaN.'])
        end
end

% Clean workspace
clear('np','i','phases','vars')
clear('N','N1','N2','phiWarning','xiWarning','PhiWarning','segments1',...
      'segments2')





