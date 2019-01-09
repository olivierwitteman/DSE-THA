function [M,MA,DOHmiss,DOHtot] = ComputeWeights(a,m,p,f,s,c,aircraft,PL,MA)
%%% Description
% This function estimates the "Class-1.5" weight breakdown of the aircraft
% for a given output of the WP-WS diagram and a determined battery and fuel
% energy required to complete the mission. A reference OEM is estimated
% as a function of TOM using Roskam's method. The weight of a reference
% (conventional aircraft) wing and gas turbine is subtracted to obtain an
% OEM which excludes wing and powerplant. The actual wing weight and
% powerplant weight are then computed using Roskam/Torenbeek methods.
%
% Input:
%   - a,m,p,f,s,c: structures containing aircraft and program data (see
%       main input file)
%   - aircraft: structure containing:
%       - Pdes: installed power of the powerplant components on the
%           aircraft, as obtained from WPdes in the WP-WS diagram [W]
%       - TOM: initial guess of the take-off mass of the aircraft [kg]
%       - Sw: wing area [m2]
%   - PL: payload [kg]
%   - MA: structure generated in MissionAnalysis.m
%
% Output:
%   - M: structure containing masses in [kg] of:
%       - f: fuel mass
%       - bat: installed battery mass
%       - bat_miss: battery mass required for nominal mission
%       - bat_E: battery mass required for total mission energy(incl. div.)
%       - bat_P: battery mass required to meet power requirements
%       - EM1: primary electrical machines (total) mass
%       - EM2: secondary electrical machines (total) mass
%       - GT: gas turbines (total) mass
%       - w: wing mass
%       - OEM: operative empty mass EXCLUDING wing and powertrain
%       - TOM: total take-off mass
%       - PL: payload (same as input)
%   - MA: update MA structure, containing the corrected mission energy
%       profiles.
%   - DOHmiss: degree-of-hybridization of nominal mission, i.e. ratio
%       between battery energy spent during nominal mission and the total
%       energy spent during the nominal mission [-]
%   - DOHtot: degree-of-hybridization of complete mission, i.e. ratio
%       between battery energy spent during complete mission and the total
%       energy spent during the complete (= incl. diversion) mission [-]
%
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 04-04-18
%%% Last modified: 16-04-18


%% Rearrange Input

% Structure containing installed power of GT and EM's
P = aircraft.Pdes;

% Initial guess for TOM
TOM0 = aircraft.TOM;

% Check if powertrain is fully electric
if strcmp(p.config,'e-1') || strcmp(p.config,'e-2') ...
                          || strcmp(p.config,'dual-e')
    electric = 1;
else
    electric = 0;
end


%% Compute battery weight
            
% Shift battery energy consumption line during descent to get a
% continuous curve along the missions
MA.Dcr.Ebat = MA.Dcr.Ebat + (MA.Dde.Ebat(1) - MA.Dcr.Ebat(end));
MA.Dcl.Ebat = MA.Dcl.Ebat + (MA.Dcr.Ebat(1) - MA.Dcl.Ebat(end));
MA.de.Ebat = MA.de.Ebat + (MA.Dcl.Ebat(1) - MA.de.Ebat(end));
MA.cr.Ebat = MA.cr.Ebat + (MA.de.Ebat(1) - MA.cr.Ebat(end));
MA.cl.Ebat = MA.cl.Ebat + (MA.cr.Ebat(1) - MA.cl.Ebat(end));

% Calculate energy-delta between minimum during nominal mission and the
% global maximum (since the min SOC limit is a fraction of the global
% maximum)
EbatTot = [MA.cl.Ebat MA.cr.Ebat MA.de.Ebat ...
                MA.Dcl.Ebat MA.Dcr.Ebat MA.Dde.Ebat];
EbatMiss = [MA.cl.Ebat MA.cr.Ebat MA.de.Ebat];
DeltaMiss = max(EbatTot) - min(EbatMiss);

% Calculate what the minimum energy at any point along the mission would be
% if the minimum value is such that the minSOC is not exceeded during the
% nominal mission
EMissMin = p.minSOC_miss/(1-p.minSOC_miss)*DeltaMiss;

% Calculate energy-delta between minimum during diversion mission and the
% global maximum
EbatDiv = [MA.Dcl.Ebat MA.Dcr.Ebat MA.Dde.Ebat];
DeltaDiv = max(EbatTot) - min(EbatDiv);

% Calculate what the minimum energy at any point along the mission would be
% if the minimum value is such that the minSOC is not exceeded during the
% diversion mission
EDivMin = p.minSOC_tot/(1-p.minSOC_tot)*DeltaDiv;

% Selected the maximum minimum-energy value, including zero, in order to
% meet the energy requirements for both conditions
EMin = max([EDivMin EMissMin 0]);

% Offset all mission segments such that EMin is not exceeded
names = fieldnames(MA);
Ncon = size(names,1);
for i = 1:Ncon
    MA.(names{i}).Ebat = MA.(names{i}).Ebat - min(EbatTot) + EMin;
end
EbatTot1 = [MA.cl.Ebat MA.cr.Ebat MA.de.Ebat ...
                MA.Dcl.Ebat MA.Dcr.Ebat MA.Dde.Ebat];
    
% Compute battery mass required to provide maximum energy along mission
M_bat_E = max(EbatTot1)/p.SE.bat;

% Compute battery mass required to provide maximum power
M_bat_P = aircraft.Pdes.bat/p.SP.bat;

% Select maximum
M_bat = max([M_bat_E M_bat_P]);


%% Calculate remaining masses which do not depend on TOM
% Note: TOM is updated in a loop outside this function!

% Fuel mass [kg]
if electric == 1
    M_f_miss = 0;
    M_f_tot = 0;
else
    
    % Shift fuel energy consumption line during descent to get a
    % continuous curve along the missions
    MA.Dcr.Ef = MA.Dcr.Ef + (MA.Dde.Ef(1) - MA.Dcr.Ef(end));
    MA.Dcl.Ef = MA.Dcl.Ef + (MA.Dcr.Ef(1) - MA.Dcl.Ef(end));
    MA.de.Ef = MA.de.Ef + (MA.Dcl.Ef(1) - MA.de.Ef(end));
    MA.cr.Ef = MA.cr.Ef + (MA.de.Ef(1) - MA.cr.Ef(end));
    MA.cl.Ef = MA.cl.Ef + (MA.cr.Ef(1) - MA.cl.Ef(end));
    
    Ef_miss = MA.cl.Ef(1)-MA.de.Ef(end);
    Ef_tot = Ef_miss + MA.Dcl.Ef(1) - MA.Dde.Ef(end);
    M_f_tot = Ef_tot/p.SE.f;
    M_f_miss = Ef_miss/p.SE.f;
end

% Mass of reference gas turbines [kg]
if electric == 1; N = 2; else N = p.N1; end;
if isnan(P.s1) && ~isnan(P.s2)
    P_ref = P.s2/N;
elseif ~isnan(P.s1) && isnan(P.s2)
    P_ref = P.s1/N;
else
    P_ref = (P.s1 + P.s2)/N;
end
M_GT_ref = N*f.W.GT(P_ref);

% Actual mass of gas turbines (using turboshaft weight correlation) [kg]
if electric == 1
    M_GT = 0;
else
    M_GT = p.N1*f.W.GT(P.GTM/p.N1);
end

% Weight of electrical machines [kg]
if ~isnan(P.EM1); M_EM1 = p.N1*f.W.EM(P.EM1/p.N1); else M_EM1 = 0; end
if ~isnan(P.EM2); M_EM2 = p.N2*f.W.EM(P.EM2/p.N2); else M_EM2 = 0; end

% Estimate DOH of nominal mission and total mission (to save in output MA
% structure)
DOHmiss = DeltaMiss/(DeltaMiss + Ef_miss);
DOHtot = (max(EbatTot1)-min(EbatTot1))/...
    ((max(EbatTot1)-min(EbatTot1)) + Ef_tot);


%% Loop till TOM convergence
err = 1;
iter = 0;
while err > s.errmax
    
    % Update counter
    iter = iter + 1;
    
    % Estimate operative empty mass (incl. wing and conventional power
    % plant)
    OEM_ref = f.W.OEM(TOM0);
    
    % Compute zero fuel weight
    MZFW = TOM0-M_f_tot;
    
    % Compute reference aircraft (no-DP) wing weight. For this assume the
    % wing loading to be determined by stall speed in unpowered landing
    % conditions
    % Estimate wing geometry of reference "conventional" aircraft
    q = 0.5*m.L.rho*m.L.vs^2;                   % Dynamic pressure [Pa]
    WS_ref = a.L.CLmax*q/m.L.f;                 % Wing loading [N/m2]
    Sw_ref = TOM0*c.g/WS_ref;                   % Wing area [m2]
    b_ref = (Sw_ref*a.AR)^0.5;                  % Wing span [m]
    tr_ref = a.tc*(Sw_ref/a.AR)^0.5*2/(1+a.TR); % Root thickness [m]
    bs_ref = b_ref/cosd(a.Lambda);              % Corrected span [m]
    
    % Reference wing mass (Torenbeek S. 8.4.1b, p.280). Ignoring weight
    % reduction due to wing-mounted engines (would require specifying
    % whether engines are or not on wing, and probably not accurate for 
    % e.g. tip-mounted engines and DP).
    M_w_ref = MZFW*6.67e-3*bs_ref^0.75*(1+(1.905/bs_ref)^0.5)*...
                a.nUlt^0.55*((bs_ref/tr_ref)/(MZFW/Sw_ref))^0.3;
    
    % OEM excl. wing and conventional powerplant
    OEM = OEM_ref - M_w_ref - M_GT_ref;
    
    % Actual wing mass 
    bs = (aircraft.Sw*a.AR)^0.5/cosd(a.Lambda);
    tr = a.tc*(aircraft.Sw/a.AR)^0.5*2/(1+a.TR);
    M_w = MZFW*6.67e-3*bs^0.75*(1+(1.905/bs)^0.5)*...
                a.nUlt^0.55*((bs/tr)/(MZFW/aircraft.Sw))^0.3;

    % Recompute TOM
    TOM1 = OEM + PL + M_bat + M_f_tot + M_w + M_GT + M_EM1 + M_EM2;
    
    % Check error and update
    err = abs((TOM1-TOM0)/TOM0);
    TOM0 = TOM1;
    
    if iter > s.itermax
        disp([s.levelString '    > TOM did not converge!'])
        TOM1 = NaN;
        break
    end
end


%% Output

% Collect function output structure
M.f = M_f_tot;
M.f_miss = M_f_miss;
M.bat = M_bat;
M.bat_miss = DeltaMiss/p.SE.bat;
M.bat_E = M_bat_E;
M.bat_P = M_bat_P;
M.EM1 = M_EM1;
M.EM2 = M_EM2;
M.GT = M_GT;
M.w = M_w;
M.OEM = OEM;
M.TOM = TOM1;
M.PL = PL;



