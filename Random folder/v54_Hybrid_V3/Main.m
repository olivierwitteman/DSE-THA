%%% Description
%
% This script combines the wing-loading power-loading diagrams and mission
% analysis to carry out the complete "Class-1.5" sizing of an HEP aircraft
% including aero-propulsive interaction effects. A description of the
% method can be found in the paper by de Vries, Brown and Vos (2018). The
% way in which aero-propulsive interaction effects are accounted for is
% briefly described in "WP_WS_diagram.m". The aero-propulsive models are
% specified in "WingPropDeltas.m". The powertrain model is described in
% "PowerTransmissionComputation_v2". The input parameters of the code are 
% described and specified in "Input.m". 
% 
% After running, the following variables (structures) should appear in the 
% workspace:
%
%   - a: contains aerodynamic properties such as the assumed CLmax or CD0
%       in the different flight conditions evaluated in the WS-WP diagram,
%       as well as wing geometries such as AR or TR. Defined by user in
%       "Input.m".
%   - AC: indicates which constraint is active/limiting for each powertrain
%       component's power loading, as well as wing loading, depending on
%       the design criteria chosen. For example, "AC.minGT.bat" gives the
%       name of the constraint which limits the minimum battery size when
%       selecting the design wing-loading corresponding to minimum gas
%       turbine size. AC is computed in "ComputeDesignPoint.m".
%   - AEROdes: provides the aerodynamic (such as CL, L/D,...) and
%       operational (such as v, M,...) properties of the aircraft in a
%       determined flight condition (= constraint) for a given design
%       point. For example, "AEROdes.minWS.cr.v" returns the velocity
%       required during the cruise constraint at the design wing loading
%       corresponding to minimum wing area (maximum W/S). AEROdes is
%       obtained in "ComputeDesignPoint.m".
%   - c: constants specified by the user in "Input.m", such as gravity or
%       sea-level atmospheric density.
%   - f: contains anonymous functions such as the power lapse of the
%       engine, atmospheric density lapse, or the weight of a powertrain
%       component as a function of installed power. Specified by user in
%       "Input.m".
%   - m: specifies operational/mission parameters, divided into flight
%       conditions (= constraints). For exmaple, "m.TO.h" is the altitude
%       at which take-off is performed. Some variables are specified by the
%       user in "Input.m"; others are computed along the way and added to
%       the structure.
%   - M: mass breakdown of the aircraft, as obtained from the mission
%       analysis. It is obtained from the "ComputeWeights.m" function.
%   - MA: collects results of the mission analysis ("MissionAnalysis.m",
%       unsuprisingly). One field exists per mission segment. For each
%       mission segment, several parameters are specified as an array, with
%       each element of the array corresponding to a timestep along the
%       mission segment. For example, "MA.Dcl.Ef" refers to the fuel energy
%       remaining on the aircraft (as a function of time) during the climb
%       phase of the diversion mission. MA also collects the resulting fuel
%       fractions and degree-of-hybridization of the aircraft.
%   - MA_in: input parameters for the mission analysis specified by the
%       user in "Input.m", such as Range, diversion altitude, and initial
%       guesses for OEM or DOH. Additional parameters are added per mission
%       segment throughout "MissionAnalysis.m".
%   - p: powertrain parameters specified by user in "Input.m", including
%       several which are flight-condition dependent. For example,
%       "p.L.etap1" refers to the propulsive efficiency of the primary
%       propulsors in landing conditions.
%   - s: program settings, such as convergence tolerances, figure numbers,
%       etc. Specified by user in "Input.m". The handles of the figures
%       generated are added to "s.figs".
%   - TW: Thrust-to-weight ratio of the different constraints evaluated in
%       the WP-WS diagram. Created in the different constraint functions
%       called in "WP_WS_diagram.m". The wing loading values at which the
%       TW arrays are specified are given by the structure "WS".
%   - TW_WSdes: Thrust-to-weight ratio evaluated at the different design
%       wing-loadings (specified in "WSdes") for each constraint. For
%       example, "TW_WSdes.minp.cr" refers to the thrust-to-weight ratio
%       obtained during cruise when selecting the wing loading
%       corresponding to minimum propulsive power as design point. TW_WSdes
%       is computed in "ComputeDesignPoint.m".
%   - WP_comp: contains the sizing power-loading values of each component 
%       of the powertrain, for each constraint and as a function of "WS".
%       For example, "WP_comp.cr.GB" gives, as a function of wing loading,
%       the sizing power-loading value (i.e. the maximum value the
%       component has to be able to produce/absorb) of the gearbox in
%       cruise conditions. This structure is generated as output of the
%       constraint functions in "WP_WS_diagram.m".
%   - WP_loss: contains the power-loading losses of each component 
%       of the powertrain, for each constraint and as a function of "WS".
%       For example, "WP_loss.cr.GB" gives, as a function of wing loading,
%       the amount of power lost due to e.g. heat dissipation (expressed
%       as a power-lodaing) in the gearbox in cruise conditions. This 
%       structure is generated as output of the constraint functions in 
%       "WP_WS_diagram.m".
%   - WP_path: contains the power-loading values of the different paths
%       that link the different components of the powertrain, for each 
%       constraint and as a function of "WS". For example, "WP_path.cr.s1" 
%       gives, as a function of wing loading, the power transmitted through
%       the primary shaft, from the gearbox to the primary propulsor 
%       (expressed as a power-loading). This structure is generated as 
%       output of the constraint functions in "WP_WS_diagram.m".
%   - WP_select: a series of WP-arrays obtained from WP_comp, WP_loss, or 
%       WP_path, which are manually selected at the end of 
%       "WP_WS_diagram.m" and passed on to "ComputeDesignPoint.m" in order
%       to evaluate the design points which lead to maximum power-loading
%       values of the fields specified (per constraint) in WP_select.
%   - WPdes: design power-loadings obtained from "ComputeDesignPoint.m" for
%       different design criteria. For example, "WPdes.minWS.GTM" refers to
%       the required power loading of the gas turbine (corrected to
%       TO/SL/max thrust conditions) when selecting the design point
%       corresponding to minimum wing area. The results in this structure
%       determine the installed power of the different powertrain 
%       components during the mission analysis, once the aircraft weight is 
%       known.
%   - WS: contains the wing-loading values at which TW and WP are sampled
%       for each constraint in the WP-WS diagram. These arrays are
%       generated in the constraint functions in "WP_WS_diagram.m".
%   - WSdes: wing-loading values of the design points selected. For
%       example, "WSdes.minbat" gives the wing loading value at which the
%       required battery power is minimized. Obtained from
%       "ComputeDesignPoint.m".
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 04-04-18
%%% Last modified: 16-04-18


%% Initialize
clc
clear all
close all
tic


%% Main body

% Load input parameters
Input;

% Check input is consistent
disp([s.levelString '> Checking powertrain input settings'])
CheckInput;

% Generate wing-loading power-loading diagrams
disp([s.levelString '> Evaluating W/S-W/P diagram'])
WP_WS_diagram;

% Run mission analysis
disp([s.levelString '> Starting Mission Analysis'])
MissionAnalysis;


%% Generate additional plots if desired

% I. Powertrain plots
% (Note: the power paths used here come from the
% constraints, not the mission analysis! So the powers shown do not
% necessarily coincide with a given point on the MA mission profile)
if s.plotPowertrain == 1
    disp([s.levelString '> Generating powertrain plots'])
    P_path.cr = structfun(@(x) 1./x/m.cr.f*M.TOM*c.g/1e6,...
        WP_path.cr,'UniformOutput',0);
    P_path.TO = structfun(@(x) 1./x/m.TO.f*M.TOM*c.g/1e6,...
        WP_path.TO,'UniformOutput',0);
    P_path.L = structfun(@(x) 1./x/m.L.f*M.TOM*c.g/1e6,...
        WP_path.L,'UniformOutput',0);
    [s.figs(end+1)] = PlotPowertrain(P_path.cr,WS.cr,...
        WSdes.(s.SelDes),s.figStart+size(s.figs,2),'%6.2f',...
        'Cruise constraint [MW]');
    [s.figs(end+1)] = PlotPowertrain(P_path.TO,WS.TO,...
        WSdes.(s.SelDes),s.figStart+size(s.figs,2),'%6.2f',...
        'Take-off constraint [MW]');
    [s.figs(end+1)] = PlotPowertrain(P_path.L,NaN,...
        NaN,s.figStart+size(s.figs,2),'%6.2f',...
        'Landing constraint [MW]');
    clear('P_path')
end

% II. Aerodynamic polar in cruise conditions
if s.Polar.plot == 1
    disp([s.levelString '> Generating cruise polar'])
    [~,~,~,~,~] = CreateLiftDragPolarMap(WSdes.(s.SelDes),'cr',a,p,f,s,1);
end

% III. Power-control envelopes
if s.Env.plot == 1
    disp([s.levelString '> Generating power-control envelope(s)'])
    [s] = CreatePowerControlEnvelope(p,s,f,c,WPdes,AEROdes,MA,M);
end

% IV. Add drag requirements to stall constraint
if s.LandingCheck.plot == 1
    disp([s.levelString '> Adding landing drag requirements'])
    [s] = CreateLandingDragRequirements_v2(a,m,p,f,s,c,WPdes);
end


%% End
disp([s.levelString '> Completed. Run time: ' num2str(toc) ' seconds'])

figures_to_close = [18, 30, 29, 27, 25, 24, 23, 21, 15, 17, 14, 13, 12, 10]

close Figure 18
close Figure 30
close Figure 29
close Figure 27
close Figure 25
close Figure 24
close Figure 23
close Figure 21
close Figure 15
close Figure 17
close Figure 14
close Figure 13
close Figure 12
close Figure 10

save Figure

close all


