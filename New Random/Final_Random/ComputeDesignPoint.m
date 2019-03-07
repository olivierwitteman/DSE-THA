function [WSdes,WPdes,TW_WSdes,AEROdes,AC,D2Wdes,a,m,p,figwp] = ComputeDesignPoint(WS,WP,TW,a,m,p,f,s)
%%% Description
% 
% This function computes the design point (in terms of WS, and WP per 
% component) for different design criteria: min wing size, min gas turbine
% size, min battery size, etc. (one per component + min wing size). It also
% computes the aerodynamic parameters corresponding to that design point,
% and plots the resulting WS-WP diagrams per component of the powertrain.
% It is assumed that there are two vertical constraints, which are the
% stall speed constraints "L" and "Liso". If more vertical lines are added,
% search for <'L'> and <'Liso'> and add the new vertical constraint there.
% Moreover, the color/marker/line settings for the graphs are defined for a
% given number of plot lines. If many constraints are added, these arrays
% have to be extended. The legend strings of the figure also have to be 
% updated. For the rest, the function can handle a generic
% number of constraints, which it detects automatically in the input
% structures WS and WP.
%
% Input:
%   - WS: structure WS.[constraintName], where each constraint contains an
%       array of values (horizontal constraint) or a single value (vertical
%       constraint). Wing loading [N/m2].
%   - WP: structure WP.[constraintName].[componentName]. Each sub-field has
%       the same amount of elements as the corresponding WS field. Power
%       loading [N/W].
%   - TW: structure TW.[constraintName2], where "constraintName2" contains
%       all constraints without distinguishing between a OEI1 and a OEI2
%       case (e.g. "bL" is a field instead of "bL1" and "bL2". Thrust
%       loading [-].
%   - a,m,p,f,s: structures containing aircraft, mission, powertrain,
%       function and settings data respectively (see WP_WS_Diagram_DP.m)
%
% Output:
%   - WSdes: structure WS.[designCriterion: provides design wing loading
%       per design criteria [N/m2]
%   - WPdes: structure WP.[designCriterion].[component]: provides design
%       power loading per component per design criteria [N/W]
%   - TWdes: structure TW.[designCriterion].[constraint]: contains thrust
%       loading values per constraint evaluated at the design wing loading
%       of a given design criteria.
%   - AEROdes: structure AEROdes.[designCriterion].[constraint].[variable]:
%       contains aerodynamic variables at the flight condition
%       corresponding to each constraint, for a given design criteria. The
%       variables included are:
%           - v: flight speed [m/s]
%           - CL_iso: isolated wing lift coefficient [-]
%           - CD_iso: isolated wing drag coefficient (parasite and 
%               lift-induced)[-]
%           - Tc: individual propulsor thrust coefficient, defined as
%               Tc = T_dp/N/(rho*v^2*D^2) [-]
%           - dCL: lift coeff. increase due to DP [-]
%           - dCD0: zero-thrust drag coeff. increase due to DP [-]
%           - dCDi: thrust-induced drag coeff. increase due to DP [-]
%           - dCL_TV: effective lift coefficient increase due to thrust
%               vectoring, i.e. T_dp*sin(gamma)/q/S [-]
%           - CL: total aircraft lift coefficient, excl. TV [-]
%           - CD: total aircraft drag coefficient, excl. TV [-]
%           - LD_iso: isolated wing lift-to-drag ratio [-]
%           - LD: powered wing lift-to-drag ratio [-]
%   - dCDi: wing (thrust-induced) drag coefficient increase due to DP [-]
%   - dCD0: wing parasite drag coefficient increase due to DP [-]
%   - AC: struct AC.[designCriterion].[component]: label of active
%       constraint at design point, per component (incl. constraint
%       limiting max wing loading value). Only provides ONE of the limiting
%       constraints at the design point, and not both in case of an
%       intersection.
%   - D2Wdes: structure D2Wdes.[designCriterion]: D2W parameter (normalized
%       propeller diameter) for the design wing loading computed for each
%       design criterion [m2/N]
%   - a,m,p: updated input a,m,p structures.
%   - figwp: figure handles
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 31-10-17
%%% Last modified: 11-04-18


%% Input settings

% Assign values for constraints that did not require certain parameters as
% input
m.Liso.f = m.L.f;
a.Liso.e = a.L.e;
a.bL.e = a.L.e;
a.Liso.CD0 = a.L.CD0;
p.Liso.T = p.L.T;


%% Compute absolute value of powertrain component WP values

% Identify number of fields
names = fieldnames(WP);
n = size(names,1);

% Number of power elements/fields (assumed to be the same for all
% constraints)
namesp = fieldnames(WP.(names{1}));
np = size(namesp,1);

% Loop over all fields and check if any negative values are found. For the
% parallel, e-1 and dual-e configuration, P_eg and P_gb should nominally be
% negative, so don't issue a warning.
for i = 1:n
    for j = 1:np
        if (strcmp(p.config,'parallel') || strcmp(p.config,'e-1') ||...
                strcmp(p.config,'dual-e')) && (strcmp(namesp{j},'e1') ||...
                strcmp(namesp{j},'gb'))
            if ~isempty(WP.(names{i}).(namesp{j})(WP.(names{i}).(namesp{j})>0))
                disp([s.levelString '    > Warning: positive WP values found '...
                      'for component ''' namesp{j} ''' in constraint ''' ...
                      names{i}  '''. Changing to positive values for sizing '...
                      'criteria.'])
            else
                WP.(names{i}).(namesp{j}) = abs(WP.(names{i}).(namesp{j}));
            end
        else
            if ~isempty(WP.(names{i}).(namesp{j})(WP.(names{i}).(namesp{j})<0))
                disp([s.levelString '    > Warning: negative WP values found '...
                      'for component ''' namesp{j} ''' in constraint ''' ...
                      names{i}  '''. Changing to positive values for sizing '...
                      'criteria.'])
                WP.(names{i}).(namesp{j}) = abs(WP.(names{i}).(namesp{j}));
            end
        end
    end
end


%% Separate into vertical and horizontal constraints

% Separate constraints into vertical and curved constraints
nh = 0;
nv = 0;
for i = 1:size(names,1)
    
    % Vertical WS constraints
    if strcmp(names{i},'L') || strcmp(names{i},'Liso')
        nv = nv + 1;
        namesv{nv} = names{i};
        WSv.(names{i}) = WS.(names{i});
        
    % Curved constraints
    else
        nh = nh + 1;
        namesh{nh} = names{i};
        WPh.(names{i}) = WP.(names{i});
        WSh.(names{i}) = WS.(names{i});
    end
end


%% Design point for max wing loading

% Wing loading limited by iso wing and powered stall constraints
for i = 1:nv
    WSCandidates(i) = WSv.(namesv{i});
end
[WSdes.minWS,idxL] = min(WSCandidates);

% Name of active constraint for wing size
AC.minWS.WS = namesv{idxL};

% For each power element, evaluate WP at the design WS for different
% constraints and select minimum
for i = 1:np
    
    % Evaluate each curved constraint
    WP_array = NaN(1,nh);
    for j = 1:nh
        
        % If the element has two or more non-NaN values for this
        % constraint, it means it has been evaluated, and thus it is
        % possible to interpolate
        if sum(~isnan(WPh.(namesh{j}).(namesp{i}))) > 1
            WP_array(j) = interp1(WSh.(namesh{j}),...
                WPh.(namesh{j}).(namesp{i}),WSdes.minWS);
        else
            WP_array(j) = NaN;
        end              
    end
    
    % Select minimum value as limiting constraint
    [WPdes.minWS.(namesp{i}),idx] = min(WP_array);
    AC.minWS.(namesp{i}) = namesh{idx};
       
end


%% Interpolate all curved constraints to common fine WS grid for comparison

% Wing loading values
WSintp = linspace(0,min([WSdes.minWS s.WSmax]),s.NWS);

% Loop over all elements
for i = 1:np
    
    % Loop over all constraints
    for j = 1:nh
        if sum(~isnan(WPh.(namesh{j}).(namesp{i}))) > 1
            WPintp.(namesh{j}).(namesp{i}) = interp1(WSh.(namesh{j}),...
                    WPh.(namesh{j}).(namesp{i}),WSintp);
        else
            WPintp.(namesh{j}).(namesp{i}) = NaN(1,s.NWS);
        end              
    end
end    


%% Design point for min component power
% The structure containing design points will be labelled as such:
% WPdes.min<<Powertrain element for which min power is achieved>>.
%          <<Name of powertrain element for which the design WP is given>>

% Loop over design criteria (one per component)
for i = 1:np
    
    % Generate label for structure
    namesdes{i} = ['min' namesp{i}];
    
    % Compute power loading for each component
    for j = 1:np
        
        % Collect all WP values in rows, one row per curved constraint
        WP_array = NaN(nh,s.NWS);
        for k = 1:nh
            WP_array(k,:) = WPintp.(namesh{k}).(namesp{j});
        end
        
        % Generate WP(WS) curve containing the minimum (i.e. the limiting)
        % values of the WP_array. Save index of the constraint which is 
        % limiting in each case.
        WP_limit.(namesp{j}) = NaN(1,s.NWS);
        idx_limit.(namesp{j}) = NaN(1,s.NWS);
        for l = 1:s.NWS
            [WP_limit.(namesp{j})(l), idx_limit.(namesp{j})(l)] = ...
                        min(WP_array(:,l));
        end
    end
    
    % Select index corresponding to the maximum of the limiting curve of 
    % the component for which the design criteria is applied, ignoring the
    % Infs
    WP_limit.(namesp{i})(WP_limit.(namesp{i})==Inf) = NaN;
    [~,idxWS] = max(WP_limit.(namesp{i}));

    % Assign wing loading value
    WSdes.(namesdes{i}) = WSintp(idxWS);

    % If component is not defined for this powertrain, assign NaN value
    if WSdes.(namesdes{i}) == 0; WSdes.(namesdes{i}) = NaN; end
    
    % Assign limiting power loading values of all components, corresponding
    % to the WS point where component i has its maximum possible WP'
    % Also save corresponding limiting constraint tags
    for j = 1:np
        WPdes.(namesdes{i}).(namesp{j}) = WP_limit.(namesp{j})(idxWS);
        AC.(namesdes{i}).(namesp{j}) = namesh{idx_limit.(namesp{j})};
    end
end


%% Compute aircraft characteristics per design criteria and per constraint

% Number of constraints ignoring the OEI1/OEI2 cases
namesaux = fieldnames(WP);
counter = 0;
for i = 1:size(namesaux,1)
    toggle1 = find(namesaux{i}=='1',1);
    toggle2 = find(namesaux{i}=='2',1);
    if ~isempty(toggle1)
        namesaux{i} = namesaux{i}(1:toggle1-1);
    end
    if isempty(toggle2)
        counter = counter+1;
        namesa{counter,1} = namesaux{i};
    end
end
na = size(namesa,1);

% Number of design criteria (one per component plus wing loading)
namesdes = fieldnames(WPdes);
nd = size(namesdes,1);

% Loop over design criteria (one per component)
for i = 1:nd
    
    % Loop over number of constraints (excl. different OEI cases)
    for j = 1:na
        
        % Generate labels to select corresponding constraint: if the label
        % doesn't exist, add a "1" to select the first of the two OEI
        % constraints
        if isfield(WP,namesa{j})
            namesb{j} = namesa{j};
        else
            namesb{j} = [namesa{j} '1'];
        end
        
        % Compute velocity
        if length(m.(namesa{j}).v) == 1
            AEROdes.(namesdes{i}).(namesa{j}).v = m.(namesa{j}).v;
        else
            WSx = WS.(namesb{j});
            v = m.(namesa{j}).v;
            WSquery = WSdes.(namesdes{i});
            AEROdes.(namesdes{i}).(namesa{j}).v = ...
                interp1(WSx(~isnan(WSx)),v,WSquery,'linear');
        end
        
        % Compute Mach number
        AEROdes.(namesdes{i}).(namesa{j}).M = ...
            AEROdes.(namesdes{i}).(namesa{j}).v / f.a(m.(namesa{j}).h);
        
        % Compute isolated wing lift coefficient
        if length(a.(namesa{j}).CL) == 1
            AEROdes.(namesdes{i}).(namesa{j}).CL_iso = a.(namesa{j}).CL;
        else
            WSx = WS.(namesb{j});
            CL_iso = a.(namesa{j}).CL;
            WSquery = WSdes.(namesdes{i});
            AEROdes.(namesdes{i}).(namesa{j}).CL_iso = ...
                interp1(WSx(~isnan(WSx)),CL_iso,WSquery,'linear');
        end

        % Compute isolated wing drag coefficient
        AEROdes.(namesdes{i}).(namesa{j}).CD_iso = f.CD(...
            a.(namesa{j}).CD0,AEROdes.(namesdes{i}).(namesa{j}).CL_iso,...
            a.AR,a.(namesa{j}).e);

        % Disk loading for constraint flight condition
        D2Wdes.(namesdes{i}) = f.D2W(WSdes.(namesdes{i}),...
            p.b_dp,p.N,p.dy,a.AR);
        D2Wcon = D2Wdes.(namesdes{i})/m.(namesa{j}).f;

        % Compute thrust loading at design point per constraint. This is
        % NOT the design/sizing thrust loading, but the thrust required per
        % constraint (which should be equal to or smaller than the limiting
        % thrust loading).
        if length(TW.(namesa{j})) == 1
            TW_WSdes.(namesdes{i}).(namesa{j}) = TW.(namesa{j});
        else
            WSx = WS.(namesb{j});
            TWx = TW.(namesa{j});
            WSquery = WSdes.(namesdes{i});
            TW_WSdes.(namesdes{i}).(namesa{j}) = ...
                interp1(WSx(~isnan(WSx)),TWx,WSquery,'linear');
        end
        
        % Correct TW to mass of aircaft during constraint, and consider
        % only the DP fraction
        TW_dp = TW_WSdes.(namesdes{i}).(namesa{j})/m.(namesa{j}).f*...
                p.(namesa{j}).T;

        % Compute thrust coefficient, defined as Tc = T/(rho*v^2*D^2). Only
        % considers the fraction of TW generated by the DP array! Expressed
        % per propulsor.
        AEROdes.(namesdes{i}).(namesa{j}).Tc = ...
            TW_dp/p.N/m.(namesa{j}).rho/...
            AEROdes.(namesdes{i}).(namesa{j}).v^2/D2Wcon;
              
        % Compute prop radius/wing chord ratio
        Rc = 0.5*(D2Wdes.(namesdes{i})*WSdes.(namesdes{i})*a.AR)^0.5;  
        
        % Call WingPropDeltas function to evaluate delta terms
        if strcmp(namesa{j},'Liso')
            dCl = 0;
            dCd0 = 0;
            dCdi = 0;
            detap = 0;
        else
            [dCl,dCd0,dCdi,detap] = WingPropDeltas(...
                AEROdes.(namesdes{i}).(namesa{j}).Tc,Rc,p.xp,a.AR,...
                AEROdes.(namesdes{i}).(namesa{j}).CL_iso,...
                AEROdes.(namesdes{i}).(namesa{j}).M,a.Lambda,...
                p.(namesa{j}).Gamma,p.(namesa{j}).etap,0);
        end

        % Save delta terms relative to full wing span
        AEROdes.(namesdes{i}).(namesa{j}).dCL = p.b_dp*dCl;
        AEROdes.(namesdes{i}).(namesa{j}).dCD0 = p.b_dp*dCd0;
        AEROdes.(namesdes{i}).(namesa{j}).dCDi = p.b_dp*dCdi;

        % Save propulsive efficiency of both powertrains
        if strcmp(namesa{j},'Liso')
            AEROdes.(namesdes{i}).(namesa{j}).detap = NaN;
            AEROdes.(namesdes{i}).(namesa{j}).etap1_iso = NaN;
            AEROdes.(namesdes{i}).(namesa{j}).etap2_iso = NaN;
            AEROdes.(namesdes{i}).(namesa{j}).etap1 = NaN;
            AEROdes.(namesdes{i}).(namesa{j}).etap2 = NaN;
        else
            AEROdes.(namesdes{i}).(namesa{j}).detap = detap;
            AEROdes.(namesdes{i}).(namesa{j}).etap1_iso = ...
                p.(namesa{j}).etap1;
            AEROdes.(namesdes{i}).(namesa{j}).etap2_iso = ...
                p.(namesa{j}).etap2;
            if p.DP == 1
                AEROdes.(namesdes{i}).(namesa{j}).etap1 = ...
                    p.(namesa{j}).etap1 + detap;
                AEROdes.(namesdes{i}).(namesa{j}).etap2 = ...
                    p.(namesa{j}).etap2;
            elseif p.DP == 2
                AEROdes.(namesdes{i}).(namesa{j}).etap1 = ...
                    p.(namesa{j}).etap1;
                AEROdes.(namesdes{i}).(namesa{j}).etap2 = ...
                    p.(namesa{j}).etap2 + detap;
            else
                AEROdes.(namesdes{i}).(namesa{j}).etap1 = ...
                    p.(namesa{j}).etap1;
                AEROdes.(namesdes{i}).(namesa{j}).etap2 = ...
                    p.(namesa{j}).etap2;
            end
        end
        
        % Compute effective CL increase due to CV
        if strcmp(namesa{j},'Liso')
            AEROdes.(namesdes{i}).(namesa{j}).dCL_TV = 0;
        else
            AEROdes.(namesdes{i}).(namesa{j}).dCL_TV = ...
                TW_dp*WSdes.(namesdes{i})*sind(p.(namesa{j}).Gamma)...
                /0.5/m.(namesa{j}).rho/...
                AEROdes.(namesdes{i}).(namesa{j}).v^2;
        end

        % Powered wing lift coefficient and drag coefficient (excl. TV!)
        AEROdes.(namesdes{i}).(namesa{j}).CL = ...
            AEROdes.(namesdes{i}).(namesa{j}).CL_iso + ...
            AEROdes.(namesdes{i}).(namesa{j}).dCL;
        AEROdes.(namesdes{i}).(namesa{j}).CD = ...
            AEROdes.(namesdes{i}).(namesa{j}).CD_iso + ...
            AEROdes.(namesdes{i}).(namesa{j}).dCD0 + ...
            AEROdes.(namesdes{i}).(namesa{j}).dCDi;

        % Compute lift-to-drag-ratio of isolated wing and powered wing
        AEROdes.(namesdes{i}).(namesa{j}).LD_iso = ...
            AEROdes.(namesdes{i}).(namesa{j}).CL_iso/...
            AEROdes.(namesdes{i}).(namesa{j}).CD_iso;
        AEROdes.(namesdes{i}).(namesa{j}).LD = ...
            AEROdes.(namesdes{i}).(namesa{j}).CL/...
            AEROdes.(namesdes{i}).(namesa{j}).CD;
    end
end

% Correction: the WP and AERO values are not known for the vertical
% constraints, unless it is the vertical constraint itself that is
% limiting. Assign NaNs to avoid someone from using the wrong value
for i = 1:nd
    if ~strcmp(namesdes{i},'minWS')
        for j = 1:nv
            TW_WSdes.(namesdes{i}).(namesv{j}) = NaN;
            namespar = fieldnames(AEROdes.(namesdes{i}).(namesv{j}));
            for k = 1:size(namespar,1)
                AEROdes.(namesdes{i}).(namesv{j}).(namespar{k}) = NaN;
            end
        end
    
    % If sized for min WS, we only know the AERO values of the limiting
    % constraint
    else
        for j = 1:nv
            if ~strcmp(namesv{j},AC.minWS.WS)
                TW_WSdes.(namesdes{i}).(namesv{j}) = NaN;
                namespar = fieldnames(AEROdes.(namesdes{i}).(namesv{j}));
                for k = 1:size(namespar,1)
                    AEROdes.(namesdes{i}).(namesv{j}).(namespar{k}) = NaN;
                end
            end
        end
    end
end


%% Generate plots
if s.plotWPWS == 1
    
    % Generate labels for powertrain components
    for i = 1:np
        switch namesp{i}
            case 'p';   labelsp{i} = 'Total propulsive';
            case 'p1';  labelsp{i} = 'Primary propulsive';
            case 'p2';  labelsp{i} = 'Secondary propulsive';
            case 's1';  labelsp{i} = 'Primary shaft';
            case 's2';  labelsp{i} = 'Secondary shaft';
            case 'e2';  labelsp{i} = 'Secondary electric machine';
            case 'f';   labelsp{i} = 'Fuel';
            case 'bat'; labelsp{i} = 'Battery';
            case 'gt';  labelsp{i} = 'Gas turbine';
            case 'gtm'; labelsp{i} = 'Corrected gas turbine';
            case 'e1';  labelsp{i} = 'Primary electric machine';
            case 'gb';  labelsp{i} = 'Gearbox';
            case 'GT';  labelsp{i} = 'Gas turbine comp.';
            case 'GTM'; labelsp{i} = 'Corrected gas turbine comp.';
            case 'GB';  labelsp{i} = 'Gearbox comp.';
            case 'PM';  labelsp{i} = 'PMAD comp.';
            case 'P1';  labelsp{i} = 'Primary propulsor comp.';
            case 'P2';  labelsp{i} = 'Secondary propulsor comp.';
            case 'EM1'; labelsp{i} = 'Primary elec. machine comp.';
            case 'EM2'; labelsp{i} = 'Secondary elec. machine comp.';
            case 'EM1M'; labelsp{i} = 'Corrected primary EM comp.';
            case 'EM2M'; labelsp{i} = 'Corrected secondary EM comp.';
            otherwise;  labelsp{i} = namesp{i};
        end
    end
    
    % Generate labels for constraints
    for i = 1:nv
        switch namesv{i}
            case 'Liso';labelsv{i} = 'Clean wing V_s';
            case 'L';   labelsv{i} = 'Powered wing V_s';
            otherwise;  labelsv{i} = namesv{i};
        end
    end
    for i = 1:nh
        switch namesh{i}
            case 'cr';  labelsh{i} = 'Cruise speed';
            case 'TO';  labelsh{i} = 'Take off (TOP)';
            case 'TO2'; labelsh{i} = 'Take off (Torenbeek 2013)';
            case 'bL1'; labelsh{i} = 'Balked landing, OEI 1';
            case 'bL2'; labelsh{i} = 'Balked landing, OEI 2';
            case 'cI1'; labelsh{i} = 'OEI 1 ceiling';
            case 'cI2'; labelsh{i} = 'OEI 2 ceiling';
            case 'cc';  labelsh{i} = 'Cruise ceiling';
            case 'cl';  labelsh{i} = 'Start-of-climb';
            case 'ct';  labelsh{i} = 'Top-of-climb';
            otherwise;  labelsh{i} = namesh{i};
        end
    end
    
    % Colors
    cols = [1 0 0;
        0 0 1;
        0 1 0;
        1 0 1;
        0 1 1;
        0.8 0.3 0;
        0.3 0 0.8;
        0.8 0 0.3;
        1 0.65 0;
        0 0.3 0.6;
        0 0.8 0.3;
        1 1 0
        0.8 0.8 0.8
        0.2 0.2 0.2];
    
    % Markers
    markers = {'o','s','d','^','v','>','<','+','*','x','p','h','o','s'};
    
    % Change plotting styles, in case figure exists already (for comparison)
    g = groot;
    if ~isempty(g.Children);
        for i = 1:length(g.Children)
            hg = g.Children(i);
            if hg.Number > s.figStart && hg.Number < s.figStart+np
                linesPlot = '--';
                transparency = 0.3;
                areacolor = [0.8 1 1];
                break
            else
                linesPlot = '-';
                transparency = 1;
                areacolor = [0.9 0.9 0.9];
                break
            end
        end
    else
        linesPlot = '-';
        transparency = 1;
        areacolor = [0.9 0.9 0.9];
    end
    
    % Generate one figure per powertrain element. Use ABSOLUTE VALUES
    figureIdx = 0;
    for i = 1:np
        
        % Only if this element exists in the powertrain configuration 
        % selected. Check all constraints, it might be that the power
        % requirement of a given component for a given constraint is zero
        % (i.e. WP = Inf -> WP = NaN) if e.g. phi is zero in that
        % constraint, but that it still is used in other flight conditions
        toggle = 0;
        for j = 1:nh
            if (any(~isnan(WP.(namesh{j}).(namesp{i}))) == 1) && ...
               (any(WP.(namesh{j}).(namesp{i})~=Inf) == 1) 
                toggle = 1;
                break
            end
        end
        if toggle == 1
            
            % Update figure index
            figureIdx = figureIdx + 1;
            
            % Create figure
            figwp(i) = figure(s.figStart-1+figureIdx);
            figwp(i).Name =[labelsp{i} ' wing-loading/power-loading diagram'];
            figwp(i).Color = [1 1 1];
            hold on;
            
            % Fill design space
            ha = area(WSintp,WP_limit.(namesp{i}),'edgecolor','none',...
                'facecolor',areacolor,'facealpha',transparency);
            strings1{1} = 'Feasible design space';
            
            % Plot vertical constraints
            for j = 1:nv
                h1(j) = plot([WSv.(namesv{j}) WSv.(namesv{j})],...
                    [0 s.WPmax],'color',cols(j,:),'linestyle',linesPlot);
                strings1{j+1} = labelsv{j};
            end
            
            % Plot curved/horizontal constraints
            for j = 1:nh
                h1(j+nv) = plot(WSh.(namesh{j}),...
                    WPh.(namesh{j}).(namesp{i}),...
                    'color',cols(j+nv,:),'linestyle',linesPlot);
                strings1{j+nv+1} = labelsh{j};
            end
            
            % Add design point for max wing loading
            h2(1) = plot(WSdes.minWS,WPdes.minWS.(namesp{i}),...
                'marker',markers{1},'markeredgecolor','k',...
                'markerfacecolor',cols(1,:),'color','none');
            strings1{nv+nh+2} = 'Design point for minimum wing size';
            
            % Add design points for max powertrain loading
            for j = 1:np
                h2(j) = plot(WSdes.(namesdes{j+1}),...
                    WPdes.(namesdes{j+1}).(namesp{i}),...
                    'marker',markers{j+1},'markeredgecolor',cols(j+1,:),...
                    'markerfacecolor','none','color','none');
                strings1{j+nv+nh+2} = ['Design point for min. ' labelsp{j} ...
                    ' power'];
            end
            
            % Axis settings
            ax = gca;
            grid on; box on;
            ax.XTick = 0:1000:s.WSmax;
            ax.YTick = 0:0.05:s.WPmax;
            ax.YAxis.TickLabelFormat = '%.2f';
            axis([0 s.WSmax 0 s.WPmax])
            ax.Layer = 'top';
            ax.XLabel.String = ('Wing loading \it{W_{\rm{TO}}\it/S}\rm [N/m^2]');
            ax.YLabel.String = ([labelsp{i} ...
                ' power loading \it{W_{\rm{TO}}\it/P_{\rm{' namesp{i} '}}\rm}\rm [N/W]']);
            legend(strings1,'Location','eastoutside');
        end
    end
else
    figwp = gobjects(1);
end


