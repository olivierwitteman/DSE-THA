function [s] = CreatePowerControlEnvelope(p,s,f,c,WPdes,AEROdes,MA,M)
%%% Description
% Function used to generate the power-control envelope, which establishes
% the bounds that the power-control parameters cannot surpass. For a SPPH
% powertrain, the envelope is generated for a series of predefined phi or
% Phi values, since the envelope is actually a 4D space: limit =
% f(available power, throttle, phi, Phi). For all other powertrains, it
% presents just one envelope. Each "envelope" shows the maximum and minimum
% allowable phi (or phi) (Z-axis) as a function of throttle (X-axis), for a
% series of altitude values (Y-axis), which determine the power available
% from the gas turbine. This code currently only considers phi and Phi 
% values in the interval [0,1].
%
% N.B.: this code is not perfected. A more "clean" and rigourous method is
% required to consistently define powertrain control limits throughout the
% mission analysis.
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 06-04-18
%%% Last modified: 16-04-18


%% Input 

% Installed power
Pdes = structfun(@(x) M.TOM*c.g./x, WPdes.(s.SelDes),'UniformOutput',0); 

% Powertrain component efficiencies
etas.GT = p.eta_GT;     % Gas turbine
etas.GB = p.eta_GB;     % Gearbox
etas.EM1 = p.eta_EM1;   % Primary electrical machine
etas.PM = p.eta_PM;     % PMAD (power management and distribution)
etas.EM2 = p.eta_EM2;   % Secondary electrical machine
etas.P1 = AEROdes.(s.SelDes).(s.Env.con).etap1;% Includes delta
etas.P2 = AEROdes.(s.SelDes).(s.Env.con).etap2;

% Powertrain configuration
switch p.config
    case 'conventional'; AnalysisCase = 0;
    case 'turboelectric'; AnalysisCase = 0;
    case 'serial'; AnalysisCase = 1;
    case 'parallel'; AnalysisCase = 1;
    case 'PTE'; AnalysisCase = 2;
    case 'SPPH'; AnalysisCase = 3;
    case 'e-1'; AnalysisCase = 0;
    case 'e-2'; AnalysisCase = 0;
    case 'dual-e'; AnalysisCase = 2;
end

% Amount of plots generated/conditions analyzed
if AnalysisCase == 3
    if length(s.Env.SPPH_phis) ~= length(s.Env.SPPH_Phis)
        error('SPPH_phis and SPPH_Phis must have the same length')
    end  
    Ncons = length(s.Env.SPPH_phis);
else
    Ncons = 1;
end
    

%% Determine altitude values to be sampled

% Mission segments
names = {'cl','cr','de','Dcl','Dcr','Dde'};
Nseg = size(names,2);

% Gather all maximum/minimum values in one array
hsMax = [];
hsMin = [];
for k = 1:Nseg
    hMax = max(MA.(names{k}).h);
    hMin = min(MA.(names{k}).h);
    hsMax = [hsMax hMax];
    hsMin = [hsMin hMin];
end

% Generate array of altitude values
h_array = linspace(min(hsMin),max(hsMax),s.Env.Nh);


%% Evaluate envelope for all rho values

% Initialize variables
X_phiMin = cell(1,Ncons);
Y_phiMin = cell(1,Ncons);
Z_phiMin = cell(1,Ncons);
X_phiMax = cell(1,Ncons);
Y_phiMax = cell(1,Ncons);
Z_phiMax = cell(1,Ncons);
X_xiMax = cell(1,Ncons);
Y_xiMax = cell(1,Ncons);
Z_xiMax = cell(1,Ncons);

% For configurations with just one control parameter, analyze phi or phi to
% generate a 3D plot
if AnalysisCase == 1 || AnalysisCase == 2
    
    % Initialize arrays
    X_phiMin{1} = NaN(s.Env.Nh,s.Env.Nxi);
    Y_phiMin{1} = NaN(s.Env.Nh,s.Env.Nxi);
    Z_phiMin{1} = NaN(s.Env.Nh,s.Env.Nxi);
    X_phiMax{1} = NaN(s.Env.Nh,s.Env.Nxi);
    Y_phiMax{1} = NaN(s.Env.Nh,s.Env.Nxi);
    Z_phiMax{1} = NaN(s.Env.Nh,s.Env.Nxi);
    X_xiMax{1} = NaN(s.Env.Nh,2);
    Y_xiMax{1} = NaN(s.Env.Nh,2);
    Z_xiMax{1} = NaN(s.Env.Nh,2);
    
    % Loop over all density values encountered during the MA    
    for i = 1:s.Env.Nh
        
        % Display progress
        disp([s.levelString '  > Evaluating altitude value ' num2str(i) ...
            ' of ' num2str(s.Env.Nh)]);
        
        % Select density value
        h = h_array(i);
        
        % Call function to generate envelope for a given rho
        [xi_phiMin,phi_phiMin,...
            xi_phiMax,phi_phiMax,...
            xi_xiMax,phi_xiMax] = ...
            ComputePowerControlLimits(Pdes,h,etas,[],[],p,s,f);
        
        % Store results in a meshgrid. Note that in the case of PTE or dual-e
        % configurations, "phi" actually refers to "Phi"! (see
        % ComputePowerControlLimits.m).
        X_phiMin{1}(i,:) = xi_phiMin;
        Y_phiMin{1}(i,:) = h*ones(size(xi_phiMin));
        Z_phiMin{1}(i,:) = phi_phiMin;
        X_phiMax{1}(i,:) = xi_phiMax;
        Y_phiMax{1}(i,:) = h*ones(size(xi_phiMax));
        Z_phiMax{1}(i,:) = phi_phiMax;
        X_xiMax{1}(i,:) = xi_xiMax;
        Y_xiMax{1}(i,:) = h*ones(size(xi_xiMax));
        Z_xiMax{1}(i,:) = phi_xiMax;
    end
    
% For SPPH, additionally loop over several phi or Phi values specified to 
% generate multiple 3D plots
elseif AnalysisCase == 3
    
    % Loop over conditions
    for j = 1:Ncons
        
        % Initialize arrays
        X_phiMin{j} = NaN(s.Env.Nh,s.Env.Nxi);
        Y_phiMin{j} = NaN(s.Env.Nh,s.Env.Nxi);
        Z_phiMin{j} = NaN(s.Env.Nh,s.Env.Nxi);
        X_phiMax{j} = NaN(s.Env.Nh,s.Env.Nxi);
        Y_phiMax{j} = NaN(s.Env.Nh,s.Env.Nxi);
        Z_phiMax{j} = NaN(s.Env.Nh,s.Env.Nxi);
        X_xiMax{j} = NaN(s.Env.Nh,2);
        Y_xiMax{j} = NaN(s.Env.Nh,2);
        Z_xiMax{j} = NaN(s.Env.Nh,2);
        
        % Select corresponding specified 
        phiSpec = s.Env.SPPH_phis(j);
        PhiSpec = s.Env.SPPH_Phis(j);
        
        % Loop over altitude values
        for i = 1:s.Env.Nh
            
            % Display progress
            disp([s.levelString '  > Evaluating altitude value ' ...
                    num2str(i) ' of ' num2str(s.Env.Nh) ', condition ' ...
                    num2str(j) ' of ' num2str(Ncons)]);
            
            % Select density value
            h = h_array(i);
            
            % Call function to generate envelope for a given rho
            [xi_phiMin,phi_phiMin,...
                xi_phiMax,phi_phiMax,...
                xi_xiMax,phi_xiMax] = ComputePowerControlLimits...
                (Pdes,h,etas,phiSpec,PhiSpec,p,s,f);
            
            % Store results in a meshgrid. 
            X_phiMin{j}(i,:) = xi_phiMin;
            Y_phiMin{j}(i,:) = h*ones(size(xi_phiMin));
            Z_phiMin{j}(i,:) = phi_phiMin;
            X_phiMax{j}(i,:) = xi_phiMax;
            Y_phiMax{j}(i,:) = h*ones(size(xi_phiMax));
            Z_phiMax{j}(i,:) = phi_phiMax;
            X_xiMax{j}(i,:) = xi_xiMax;
            Y_xiMax{j}(i,:) = h*ones(size(xi_xiMax));
            Z_xiMax{j}(i,:) = phi_xiMax;
        end
    end  
end


%% Generate plots

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'Power-control envelope(s)';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% For throttle-only envelopes
if AnalysisCase == 0
    hold on; grid on; box on
    xlabel('Altitude [m]');    
    ylabel('GT throttle \xi [-]');
    H0 = plot([min(h_array) max(h_array)],[1 1],'g');
    H1 = plot(MA.cl.h,MA.cl.P.xi,'vk','markersize',6,'markerfacecolor','y');
    H2 = plot(MA.cr.h,MA.cr.P.xi,'ok','markersize',6,'markerfacecolor','y');
    H3 = plot(MA.de.h,MA.de.P.xi,'sk','markersize',6,'markerfacecolor','y');
    H4 = plot(MA.Dcl.h,MA.Dcl.P.xi,'^k','markersize',6,'markerfacecolor','y');
    H5 = plot(MA.Dcr.h,MA.Dcr.P.xi,'dk','markersize',6,'markerfacecolor','y');
    H6 = plot(MA.Dde.h,MA.Dde.P.xi,'pk','markersize',8,'markerfacecolor','y');
    legend([H0 H1 H2 H3 H4 H5 H6],'\xi_{max}',...
        'Climb','Cruise','Descent','Div. climb',...
        'Div. cruise','Div. descent')
    
% For phi/Phi cases
elseif AnalysisCase == 1 || AnalysisCase == 2
    hold on; grid on; box on
    xlabel('GT throttle \xi [-]');
    ylabel('Altitude [m]');
    if AnalysisCase == 1
        zlabel('Supplied power ratio \phi  [-]')
    elseif AnalysisCase == 2
        zlabel('Shaft power ratio \Psi  [-]')
    end
    
    % Envelope
    h1 = surf(X_phiMin{1},Y_phiMin{1},Z_phiMin{1},'facecolor','b','edgecolor',...
        'none','FaceAlpha',0.6);
    h2 = surf(X_phiMax{1},Y_phiMax{1},Z_phiMax{1},'facecolor','r','edgecolor',...
        'none','FaceAlpha',0.6);
    h3 = surf(X_xiMax{1},Y_xiMax{1},Z_xiMax{1},'facecolor','g','edgecolor',...
        'none','FaceAlpha',0.6);
    
    % Mission profiles
    if AnalysisCase == 0
        disp([s.levelString '  > Power control envelope is not '...
            'applicable to a ' p.config ' configuration.'])
        varplot = 'phi';
    elseif AnalysisCase == 1
        varplot = 'phi';
    elseif AnalysisCase == 2
        varplot = 'Phi';
    end
    H1 = plot3(MA.cl.P.xi,MA.cl.h,MA.cl.P.(varplot),'vk','markersize',6,'markerfacecolor','y');
    H2 = plot3(MA.cr.P.xi,MA.cr.h,MA.cr.P.(varplot),'ok','markersize',6,'markerfacecolor','y');
    H3 = plot3(MA.de.P.xi,MA.de.h,MA.de.P.(varplot),'sk','markersize',6,'markerfacecolor','y');
    H4 = plot3(MA.Dcl.P.xi,MA.Dcl.h,MA.Dcl.P.(varplot),'^k','markersize',6,'markerfacecolor','y');
    H5 = plot3(MA.Dcr.P.xi,MA.Dcr.h,MA.Dcr.P.(varplot),'dk','markersize',6,'markerfacecolor','y');
    H6 = plot3(MA.Dde.P.xi,MA.Dde.h,MA.Dde.P.(varplot),'pk','markersize',8,'markerfacecolor','y');
    
    % Add legend
    if AnalysisCase == 1
        legend([h1 h2 h3 H1 H2 H3 H4 H5 H6],'\phi_{min}','\phi_{max}',...
            '\xi_{max}','Climb','Cruise','Descent','Div. climb',...
            'Div. cruise','Div. descent')
    elseif AnalysisCase == 2
        legend([h1 h2 h3 H1 H2 H3 H4 H5 H6],'\Psi_{min}','\Psi_{max}',...
            '\xi_{max}','Climb','Cruise','Descent','Div. climb',...
            'Div. cruise','Div. descent')
    end
    
% SPPH plots
elseif AnalysisCase == 3
    
    % Determine subplots layout
    if Ncons > 2
        Nrows = 2;
        Ncols = ceil(Ncons/2);
    else
        Nrows = 1;
        Ncols = Ncons;
    end
    
    % Loop over different conditions
    for j = 1:Ncons
        
        % Determine whether this subplot is at constant phi or at constant
        % Phi
        if ~isnan(s.Env.SPPH_phis(j)) && isnan(s.Env.SPPH_Phis(j))
            plotCase = 2;   % If phi is fixed
            varplot = 'Phi';
            OppVarplot = 'phi';
        elseif isnan(s.Env.SPPH_phis(j)) && ~isnan(s.Env.SPPH_Phis(j))
            plotCase = 1;   % If Phi is fixed
            varplot = 'phi';
            OppVarplot = 'Phi';
        else
            error('Incorrect s.Env.SPPH_phis/Phis specified')
        end
        
        % New subplot
        subplot(Nrows,Ncols,j)
        hold on; grid on; box on
        xlabel('GT throttle \xi [-]');
        ylabel('Altitude [m]');
        if plotCase == 1
            zlabel('Supplied power ratio \phi  [-]')
        elseif plotCase == 2
            zlabel('Shaft power ratio \Psi  [-]')
        end
        
        % Envelope
        h1 = surf(X_phiMin{j},Y_phiMin{j},Z_phiMin{j},'facecolor','b','edgecolor',...
            'none','FaceAlpha',0.6);
        h2 = surf(X_phiMax{j},Y_phiMax{j},Z_phiMax{j},'facecolor','r','edgecolor',...
            'none','FaceAlpha',0.6);
        h3 = surf(X_xiMax{j},Y_xiMax{j},Z_xiMax{j},'facecolor','g','edgecolor',...
            'none','FaceAlpha',0.6);
        
        % Mission profiles
        H1 = plot3(MA.cl.P.xi,MA.cl.h,MA.cl.P.(varplot),'vk','markersize',6,'markerfacecolor','y');
        H2 = plot3(MA.cr.P.xi,MA.cr.h,MA.cr.P.(varplot),'ok','markersize',6,'markerfacecolor','y');
        H3 = plot3(MA.de.P.xi,MA.de.h,MA.de.P.(varplot),'sk','markersize',6,'markerfacecolor','y');
        H4 = plot3(MA.Dcl.P.xi,MA.Dcl.h,MA.Dcl.P.(varplot),'^k','markersize',6,'markerfacecolor','y');
        H5 = plot3(MA.Dcr.P.xi,MA.Dcr.h,MA.Dcr.P.(varplot),'dk','markersize',6,'markerfacecolor','y');
        H6 = plot3(MA.Dde.P.xi,MA.Dde.h,MA.Dde.P.(varplot),'pk','markersize',8,'markerfacecolor','y');
        
        % Add text labels indicating at which phi/psi mission segments are
        % carried out, since the values plotted in the chart are not per se
        % the ones dipslayed in the title
        for k = 1:Nseg
            if plotCase == 1
                text(MA.(names{k}).P.xi(1),MA.(names{k}).h(1),MA.(names{k}).P.(varplot)(1),...
                    ['  \Psi = ' num2str(MA.(names{k}).P.(OppVarplot)(1),'%.2f')],...
                    'HorizontalAlignment','Left');
                text(MA.(names{k}).P.xi(end-1),MA.(names{k}).h(end-1),MA.(names{k}).P.(varplot)(end-1),...
                    ['\Psi = ' num2str(MA.(names{k}).P.(OppVarplot)(end-1),'%.2f') '  '],...
                    'HorizontalAlignment','Right');
            elseif plotCase == 2
                text(MA.(names{k}).P.xi(1),MA.(names{k}).h(1),MA.(names{k}).P.(varplot)(1),...
                    ['  \phi = ' num2str(MA.(names{k}).P.(OppVarplot)(1),'%.2f')],...
                    'HorizontalAlignment','Left');
                text(MA.(names{k}).P.xi(end-1),MA.(names{k}).h(end-1),MA.(names{k}).P.(varplot)(end-1),...
                    ['\phi = ' num2str(MA.(names{k}).P.(OppVarplot)(end-1),'%.2f') '  '],...
                    'HorizontalAlignment','Right');
            end
        end
        
        % Add tile
        if plotCase == 1
            title(['Envelope at constant \Psi = ' num2str(s.Env.SPPH_Phis(j),'%.2f')])
        elseif plotCase == 2
            title(['Envelope at constant \phi = ' num2str(s.Env.SPPH_phis(j),'%.2f')])        
        end
    end
    
    % Add legend
    legend([h1 h2 h3 H1 H2 H3 H4 H5 H6],'\phi_{min} or \Psi_{min}',...
        '\phi_{max} or \Psi_{max}',...
        '\xi_{max}','Climb','Cruise','Descent','Div. climb',...
        'Div. cruise','Div. descent')
end






