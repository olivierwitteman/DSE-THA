function [WP_path,WP_comp,WP_loss] = SizePowertrain(WP_in,con,OEI,m,p,f,s)
%%% Description
%
% This function takes a series of total propulsive power loading values 
% for a given flight condition and computes the power loading required from
% each component of the powertrain, which depends on the power control
% parameters. The batteries are assumed to be sized by energy, and since
% their performance depends on the integral behaviour along the mission,
% they are not explicitly sized here (i.e. the state-of-charge or discharge
% efficiency are not taken into accoun). The components can be sized taking
% into account failure of any one component in the powertrain, except the
% PMAD, which is assumed not to fail (due to redundant wiring).
%
% Input:
%   - WP_in: array containing a series of power-loading values [N/W]
%   - con: string indicating the constraint, e.g. 'cr','TO', etc.
%   - OEI: integer value, indicating whether a component has failed
%       - 0: nominal operation, no failure
%       - 1: primary powertrain component failure
%       - 2: secondary powertrain component failure
%   - m,p,f,s: structures containing mission parameters, powertrain
%       parameters, anonymous functions, and program settings,
%       respectively (see WP_WS_diagram.m).
%
% Output:
%   - WP_path: structure array containing the power loading corresponding
%       to each of the 10 paths along the powertrain. Sizes the 
%       transmission elements (e.g. electrical cables)
%   - WP_comp: structure array containing the power loading that each of
%       the 7 components has to be able to provide in this flight
%       condition, including the gas turbine power corrected to max thrust
%       and sea-level conditions ('gtm'). Sizes the components themselves.
%   - WP_loss: structure array containing the power (loading) losses
%       produced at each element. Sizes the cooling systems associated to
%       each component (except for the propulsors, where this figure is
%       indicative of the aerodynamic losses).
%
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 30-01-18
%%% Last modified: 11-02-18


%% Rearrange input variables

% Use power per unit weight [W/N]
PW = 1./WP_in;

% Density
rho = m.(con).rho;

% Powertrain component efficiencies
etas.GT = p.eta_GT;     % Gas turbine
etas.GB = p.eta_GB;     % Gearbox
etas.EM1 = p.eta_EM1;   % Primary electrical machine
etas.PM = p.eta_PM;     % PMAD (power management and distribution)
etas.EM2 = p.eta_EM2;   % Secondary electrical machine

% Control parameters
phi = m.(con).phi;
Phi = m.(con).Phi;
xi = m.(con).t;         % GT/EM throttle (note: should be renamed,
                        %   originally "t" for throttle, now "xi" like
                        %   in the paper.
% Check throttle is defined
if isnan(xi) || isempty(xi)
    error('Throttle must be specified to size the powertrain')
end                        
                        
% Check etas have correct values
compnames = fieldnames(etas);
for i = 1:size(etas,1)
    if etas.(compnames{i}) > 1 || etas.(compnames{i}) < 0
        error('Component efficiencies must have values between 0 and 1.')
    end
end
                        

%% Compute power paths, component powers, and power losses

% Loop over all wing/power loading values
for i = 1:length(PW)
    
    % Correct propulsive efficiency of DP propulsors
    if p.DP == 1
        etas.P1 = p.(con).etap1 + p.(con).detap(i);
        etas.P2 = p.(con).etap2;
    elseif p.DP == 2
        etas.P1 = p.(con).etap1;
        etas.P2 = p.(con).etap2 + p.(con).detap(i);
    else 
        etas.P1 = p.(con).etap1;
        etas.P2 = p.(con).etap2;
    end
    
    % Verify no negative propulsive power has been obtained
    if (etas.P1 > 1) || (etas.P1 < 0) || (etas.P2 > 1) || (etas.P2 < 0) 
        error('Propulsive efficiencies must have values between 0 and 1.')
    end
    
    % Call power transmission function, neglecting throttle!
    [PW_out,~,phi_out,Phi_out,solution] = ...
        PowerTransmissionComputation_v2...
        (p.config,etas,phi,Phi,[],PW(i),[],[],[]);

    % Initialize variables
    if i == 1; 
        pathnames = fieldnames(PW_out);
        string1 = [];
        for j = 1:size(pathnames,1)
            if j == size(pathnames,1)
                string1 = strcat(string1,'''', pathnames{j},''',NaN(size(WP_in))');
            else
                string1 = strcat(string1,'''', pathnames{j},''',NaN(size(WP_in)),');
            end
        end
        eval(['PW_path = struct(' string1 ');'])
    end
    
    % Store in array
    for j = 1:size(pathnames,1)
        PW_path.(pathnames{j})(i) = PW_out.(pathnames{j});
    end
    
    % Issue warning if negative power ratios have been found
    if (~(phi<0) && phi_out<0) || (~(phi>0) && phi_out>1)
        disp([s.levelString '  > Warning: phi = ' num2str(phi_out) ...
              ' for ''' con ''' condition (battery is charging).'...
              ' Solution array = ' num2str(solution')]);
    end
    if (~(Phi<0) && Phi_out<0) || (~(Phi>0) && Phi_out>1)
        disp([s.levelString '  > Warning: Phi = ' num2str(Phi_out) ...
              ' for ''' con ''' condition (propulsors are windmilling.'...
              ' Solution array = ' num2str(solution')]);
    end

    % Compute losses and sizing condition of each component
    %
    % The gas turbine is sized in terms of the shaft power it provides, since
    % that is the convention, and the power cannot flow in opposite direction.
    % For the other components (GB, EG, PM, EM, P1, P2), the "sizing" power is
    % assumed to be the summation of powers flowing into the element, which
    % will always be greater than the summation of powers flowing out of the
    % element. The difference between the ones going in and out is the loss of
    % that component, which will have to be accounted for with e.g. cooling.
    PW_comp.GT(i) = PW_path.gt(i);
    PW_comp.EM1(i) = max([abs(PW_path.gb(i)) abs(PW_path.e1(i))]);
    PW_comp.EM2(i) = max([abs(PW_path.e2(i)) abs(PW_path.s2(i))]);
    PW_comp.P1(i) = max([abs(PW_path.p1(i)) abs(PW_path.s1(i))]);
    PW_comp.P2(i) = max([abs(PW_path.p2(i)) abs(PW_path.s2(i))]);
    
    % For the GB and PMAD nodes, take all powers flowing into the node as
    % positive. Then add the ones with positive values.
    P1_GB = PW_path.gt(i); P2_GB = -PW_path.s1(i); P3_GB = -PW_path.gb(i);
    P1_PM = PW_path.e1(i); P2_PM = PW_path.bat(i); P3_PM = -PW_path.e2(i);
    P_GB = [P1_GB P2_GB P3_GB];
    P_PM = [P1_PM P2_PM P3_PM];
    if sum(isnan(P_GB)) == 3
        PW_comp.GB(i) = NaN;
    else
        PW_comp.GB(i) = sum(P_GB(P_GB>0));
    end
    if sum(isnan(P_PM)) == 3
        PW_comp.PM(i) = NaN;
    else
        PW_comp.PM(i) = sum(P_PM(P_PM>0));
    end
    
    % Compute losses as difference between inflow and outflow at each component
    PW_loss.GT(i) = abs(abs(PW_path.gt(i))-abs(PW_path.f(i)));
    PW_loss.EM1(i) = abs(abs(PW_path.gb(i))-abs(PW_path.e1(i)));
    PW_loss.EM2(i) = abs(abs(PW_path.e2(i))-abs(PW_path.s2(i)));
    PW_loss.P1(i) = abs(abs(PW_path.p1(i))-abs(PW_path.s1(i)));
    PW_loss.P2(i) = abs(abs(PW_path.p2(i))-abs(PW_path.s2(i)));
    if sum(isnan(P_GB)) == 3
        PW_loss.GB(i) = NaN;
    else
        PW_loss.GB(i) = sum(P_GB(P_GB>0)) + sum(P_GB(P_GB<0));
    end
    if sum(isnan(P_PM)) == 3
        PW_loss.PM(i) = NaN;
    else
        PW_loss.PM(i) = sum(P_PM(P_PM>0)) + sum(P_PM(P_PM<0));
    end    
    
    % For consistency, check that no negative losses are found
    for j = 1:size(compnames,1)
        if PW_loss.(compnames{j}) < 0
            error(['Negative power losses found in: ' compnames{j}])
        end
    end
end


%% Correct for OEI

% For powertrains with only primary components, overwrite N2 with N1. The
% opposite for powertrains with only secondary components. 
if strcmp(p.config,'parallel') || strcmp(p.config,'conventional')
    N2 = p.N1;
    N1 = p.N1;
elseif strcmp(p.config,'e-2')
    N2 = p.N2;
    N1 = p.N1;
else
    N2 = p.N2;
    N1 = p.N1;
end

% Select correction factors for OEI depending on whether the primary or
% secondary powertrain has failed
if OEI == 0
    k1 = 1;
    k2 = 1;
elseif OEI == 1
    k1 = N1/(N1-1);
    k2 = 1;
elseif OEI == 2
    k1 = 1;
    k2 = N2/(N2-1);
else
    error('Incorrect OEI condition specified')
end

% Correct power required from each component in case of component failure.
% It is assumed that the PMAD does not fail (requires redundant wiring).
% The losses are also affected by the component failure (although in
% practice it may be viable to have an under-sized cooling system since the
% aircraft will not finish the complete mission).
PW_comp.GT = PW_comp.GT*k1; PW_loss.GT = PW_loss.GT*k1; 
PW_comp.EM1 = PW_comp.EM1*k1; PW_loss.EM1 = PW_loss.EM1*k1;
PW_comp.EM2 = PW_comp.EM2*k2; PW_loss.EM2 = PW_loss.EM2*k2;
PW_comp.P1 = PW_comp.P1*k1; PW_loss.P1 = PW_loss.P1*k1;
PW_comp.P2 = PW_comp.P2*k2; PW_loss.P2 = PW_loss.P2*k2;
PW_comp.GB = PW_comp.GB*k1; PW_loss.GB = PW_loss.GB*k1;
PW_path.gt = PW_path.gt*k1; 
PW_path.e1 = PW_path.e1*k1; 
PW_path.e2 = PW_path.e2*k2;
PW_path.s1 = PW_path.s1*k1; 
PW_path.s2 = PW_path.s2*k2;
PW_path.gb = PW_path.gb*k1; 
PW_path.p1 = PW_path.p1*k1; 
PW_path.p2 = PW_path.p2*k2; % Note that p1 and p2 are also corrected for 
                            % OEI, since they may size the fans in
                            % windmilling conditions

% Battery failure is not considered, since this depends on the mission
% profile (batteries are assumed to be sized for energy, while the power
% requirement can be obtained by placing the battery cells in series).


%% Compute powers corrected to SL, max power conditions

% For fully-electrical configurations, only correct for throttle
% Electrical configurations where throttle affects EM1
if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
    
    % Correct electrical-machine installed power for throttle
    PW_comp.EM1M = PW_comp.EM1/xi;
    PW_loss.EM1M = PW_comp.EM1M*(1/etas.EM1-1);
    
    % Assign NaNs to other motors
    PW_comp.GTM = NaN(size(PW_comp.GT));
    PW_loss.GTM = NaN(size(PW_comp.GT));
    PW_comp.EM2M = NaN(size(PW_comp.GT));
    PW_loss.EM2M = NaN(size(PW_comp.GT));
    
% Electrical configurations where throttle affects EM2    
elseif strcmp(p.config,'e-2')
    
    % Correct electrical-machine installed power for throttle
    PW_comp.EM2M = PW_comp.EM2/xi;
    PW_loss.EM2M = PW_comp.EM2M*(1/etas.EM1-1);
    
    % Assign NaNs to other motors
    PW_comp.GTM = NaN(size(PW_comp.GT));
    PW_loss.GTM = NaN(size(PW_comp.GT));
    PW_comp.EM1M = NaN(size(PW_comp.GT));
    PW_loss.EM1M = NaN(size(PW_comp.GT));
    
% For powertrains containing gas turbines    
else
    % Correct to max thust settings in given flight condition
    PW_comp.GTM = PW_comp.GT/xi;
    
    % Correct to SL conditions. Not corrected for MTOW, since this is done
    % outside this function
    PW_comp.GTM = PW_comp.GTM/f.Alpha(rho);
       
    % Include losses in max-thrust setting for consistency of the arrays
    PW_loss.GTM = PW_comp.GTM*(1/etas.GT-1);
    
    % Assign NaNs to other motors
    PW_comp.EM1M = NaN(size(PW_comp.GT));
    PW_loss.EM1M = NaN(size(PW_comp.GT));
    PW_comp.EM2M = NaN(size(PW_comp.GT));
    PW_loss.EM2M = NaN(size(PW_comp.GT));
end


%% Convert back to power loading

% Update fieldnames to include gtm
compnames = fieldnames(PW_comp);
pathnames = fieldnames(PW_path);
lossnames = fieldnames(PW_loss);

% Component power loadings, and associated losses
for i = 1:size(compnames,1)
    WP_comp.(compnames{i}) = 1./PW_comp.(compnames{i});
end

for i = 1:size(lossnames,1)
    WP_loss.(compnames{i}) = 1./PW_loss.(compnames{i});
end

% Path power loadings
for i = 1:size(pathnames,1)
    WP_path.(pathnames{i}) = 1./PW_path.(pathnames{i});
end







