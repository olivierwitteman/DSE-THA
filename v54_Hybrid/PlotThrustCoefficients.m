function [fig,warnings] = PlotThrustCoefficients(WS,WP,WSdes,WPdes,a,m,p,f,s)
%%% Description
%
% This function computes the thrust coefficients required from the DP
% propulsors ONLY (expressed per 1 propulsor) for the different flight
% conditions analyzed, and checks whether the design points identified in
% the WP-WS diagram lead to excessively high thrust coefficient values.
%
% Input:
%   - WS, WP: wing and power loading, see WP_WS_diagram_DP.m
%   - WSdes, WPdes: wing and power loading values at the design points
%       identified in ComputeDesignPoint.m
%   - a,m,p,f,s: structures containing aircraft, mission and powertrain
%       parameters, anonymous functions, and program settings,
%       respectively. See WP_WS_diagram_DP.m. 
%   - plotOption: a figure showing the Tc maps and constraints is generated
%       if set to 1.
%
% Output:
%   - fig: figure handle
%   - warnings: 1D array containing ones for the
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 06-02-18
%%% Last modified: 08-02-18


%% Input settings

% Resolution increase in plots, to generate smooth fills
resfactor = 2;

% Maximum thrust coefficient (defined as Tc = T/rho/v^2/D^2)
Tcmax = s.Tcmax;

% Number of constraints ignoring the OEI1/OEI2 cases and unpowered
% constraints
if isfield(WP,'Liso'); WP = rmfield(WP,'Liso'); end
namesFields = fieldnames(WP);
counter = 0;
for k = 1:size(namesFields,1)
    toggle1 = find(namesFields{k}=='1',1);
    toggle2 = find(namesFields{k}=='2',1);
    trig = 0;
    if ~isempty(toggle1)
        namesFields{k} = namesFields{k}(1:toggle1-1);
        trig = 1;
    end
    if isempty(toggle2)
        counter = counter+1;
        namesa{counter,1} = namesFields{k};
        if trig == 1; OEIcase(counter) = 1; else OEIcase(counter) = 0; end;
    end
end
na = size(namesa,1);


%% Loop over constraint being evaluated
WP_array = linspace(0,s.WPmax,resfactor*s.n);
WS_array = linspace(0,s.WSmax,resfactor*s.n);
for x = 1:na

    % Select constraint
    con = namesa{x};
    
    % AEO constraints: extract and translate into flight conditions
    if OEIcase(x) == 0
        WP_con{x} = WP.(namesa{x}).p*m.(namesa{x}).f/p.(namesa{x}).T;
        WS_con{x} = WS.(namesa{x})*m.(namesa{x}).f;
        
        % OEI constraints
    elseif OEIcase(x) == 1
        WS_con{x} = WS.([namesa{x} '1'])*m.(namesa{x}).f;
        WP_con{x} = WP.([namesa{x} '1']).p*m.(namesa{x}).f/p.(namesa{x}).T;
    end
    
    % Initialize arrays for Tc limit curve and map
    Tc_map{x} = NaN(resfactor*s.n,resfactor*s.n);
    WP_lim{x} = NaN(1,resfactor*s.n);
    
    % Loop over WS (i) and WP(j) values
    for i = 1:length(WS_array)
        
        % Use velocity corresponding to current wing loading value
        if length(m.(con).v) > 1
            v = interp1(WS_con{x},m.(con).v,WS_array(i),'linear');%(m.(con).v(i);
        else
            v = m.(con).v;
        end
        
        % Limiting curve
        D2W = f.D2W(WS_array(i),p.b_dp,p.N,p.dy,a.AR);
        TW_aux = Tcmax/(p.(con).T/p.N/m.(con).rho/v^2/D2W);
        WP_lim{x}(i) = 1/TW_aux/v;
        
        % Contour map
        for j = 1:length(WP_array)
            TW_aux = 1/WP_array(j)/v;
            Tc_map{x}(i,j) = p.(con).T*TW_aux/p.N/m.(con).rho/v^2/D2W;
        end
    end
end


%% Establish values of mesh grid where design is valid.
% Distinguish between vertical constraints and inclined
% constraints

% Identify design points' names
namesdes = fieldnames(WSdes);

% Generate matrix indicating feasible design region
[X,Y] = meshgrid(WS_array,WP_array);
A = cell(1,na);

% Loop over design constraints
for x = 1:na
    
    % Select constraint
    con = namesa{x};
    A{x} = NaN(size(X));
    
    % Vertical constraint
    if length(WS_con{x}) == 1
        
        % Loop over columns (X-coordinate) in matrix
        for i = 1:length(WS_array)
            
            % Obtain WP_lim value at this wing loading
            WPmin = interp1(WS_array,WP_lim{x},X(1,i),'linear');
            
            % Loop over rows (Y-coordinate)
            for j = 1:length(WP_array)
                        
                % Check if coordinate is in feasible region
                if (X(j,i) <= WS_con{x}) && (Y(j,i) >= WPmin)
                    A{x}(j,i) = 1;
                end
            end
        end
        
    % Inclined constraints
    else
        % Loop over columns (X-coordinate) in matrix
        for i = 1:length(WS_array)
            
            % Obtain WP_lim value and constraint at this wing loading
            WPmin = interp1(WS_array,WP_lim{x},X(1,i),'linear');
            WPmax = interp1(WS_con{x},WP_con{x},X(1,i),'linear');
            
            % Loop over rows (Y-coordinate)
            for j = 1:length(WP_array)
                        
                % Check if coordinate is in feasible region
                if (Y(j,i) <= WPmax) && (Y(j,i) >= WPmin)
                    A{x}(j,i) = 1;
                end
            end
        end    
    end
    
    % Assign NaNs to Tc_map in infeasible region
    Tc_map{x} = Tc_map{x}';
    Tc_map{x}(isnan(A{x})) = NaN;
    
    % Compute design points to include in graph and check if they violate
    % the Tc constraint
    warnings = zeros(size(namesdes));
    for n = 1:size(namesdes,1)
        
        % If design points have been specified in terms of total propulsive
        % power (path)
        if isfield(WPdes.(namesdes{n}),'p')
            WP_des{x}(n) = WPdes.(namesdes{n}).p*m.(namesa{x}).f/p.(namesa{x}).T;
            
        % If design points have been specified only in terms of e.g.
        % component powers
        elseif isfield(WPdes.(namesdes{n}),'P1') || isfield(WPdes.(namesdes{n}),'P1') 
            if p.DP == 1
                WP_des{x}(n) = WPdes.(namesdes{n}).P1*m.(namesa{x}).f/p.(namesa{x}).T;
            elseif p.DP == 2
                WP_des{x}(n) = WPdes.(namesdes{n}).P2*m.(namesa{x}).f/p.(namesa{x}).T;
            else
                WP_des{x}(n) = NaN;
            end
        else
            error(['The total propulsive power must be included in the '...
                'list of power-loadings fed to ComputeDesignPoint.m if '...
                'the thrust coefficient is checked.'])
        end
        WS_des{x}(n) = WSdes.(namesdes{n})*m.(namesa{x}).f;
        
        % Interpolate to find limit Tc value at the WS of the design point
        WPlim = interp1(WS_array,WP_lim{x},WS_des{x}(n),'linear');
        
        % Check if design point violates constraint
        if WP_des{x}(n) < WPlim; warnings(n) = 1; end
    end
    
    % Issue warning if violating design points have been found
    if sum(warnings)>0
        disp([s.levelString '  > Warning: ' num2str(sum(warnings)) ...
              ' design point(s) lead(s) to excessively high thrust '...
              'coefficients (Tc > ' num2str(Tcmax) ') in ''' namesa{x} ...
              ''' conditions. Propulsor disk area must be increased.'])
    end
end


%% Plot
if s.plotTc == 1
    
    % Generate figure
    fig = figure(s.figStart+size(s.figs,2));
    fig.Name = 'Thrust coefficients of DP system';
    fig.Color = [1 1 1];
        
    % Loop over constraints (one per subplot)
    for x = 1:na
        subplot(2,ceil(na/2),x)
        hold on;
        
        % Plot Tc limit and contour
        h0 = surf(WS_array,WP_array,Tc_map{x},'edgecolor','none');
        h1 = plot3(WS_array,WP_lim{x},Tcmax*ones(size(WS_array)),'-b','linewidth',2);
        
        % Vertical constraints
        if length(WS_con{x})==1
            h2 = plot3([WS_con{x} WS_con{x}],[0 s.WPmax],...
                [Tcmax Tcmax],'-g','linewidth',2);
            
            % Inclined constraints
        else
            h2 = plot3(WS_con{x},WP_con{x},Tcmax*...
                ones(size(WP_con{x})),'-g','linewidth',2);
        end
        
        % Add design points
        for n = 1:size(namesdes,1)
            h3 = plot3(WS_des{x}(n),WP_des{x}(n),Tcmax+1,'.',...
                'markeredgecolor','r','markerfacecolor','r');
        end
        
        % Graph settings
        axis([0 s.WSmax 0 s.WPmax 0 Tcmax+1]);
        grid on; box on;
        caxis([0 Tcmax])
        colormap(flip(gray))
        ax =gca; ax.Layer = 'Top';
        ax.YAxis.TickLabelFormat = '%.2f';
        xlabel('Wing loading \it{W/S}\rm [N/m^2]')
        if x == 1 || x == ceil(na/2) + 1
            ylabel('Total propulsive power loading \it{W/P_{\rm{p}}}\rm [N/W]')
        end
        if x == na
            legend([h0 h1 h2 h3],'Feasible design space','Maximum T_c',...
                'Performance constraint','Design points')
        end
        
        % Title
        switch namesa{x}
            case 'cr';  title('Cruise speed');
            case 'TO';  title('Take off (TOP)');
            case 'TO2'; title('Take off (Torenbeek 2013)');
            case 'bL';  title('Balked landing');
            case 'cI';  title('OEI ceiling');
            case 'cc';  title('Cruise ceiling');
            case 'cl';  title('Start-of-climb');
            case 'ct';  title('Top-of-climb');
            case 'L';   title('V_s, landing, powered');
            otherwise; title(namesa{x});
        end
    end
    
% If no figure should be created    
else
    fig = gobjects(1);
end

