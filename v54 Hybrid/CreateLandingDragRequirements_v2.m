function [s] = CreateLandingDragRequirements_v2(a,m,p,f,s,c,WPdes)
%% Landing operating conditions

% Density
m.L.rho = f.rho(m.L.h);

% Velocity/dynamic pressure
m.L.v = m.L.vs;
m.L.M = m.L.v/f.a(m.L.h);
q = 0.5*m.L.rho*m.L.v^2; 
           
% No horizontal acceleration
dvdt = 0;               

% No climb rate, no climb gradient
climb = [0 NaN];        

% No bank angle/load factor req./turn radius
maneuver = [0 NaN NaN]; 

% CL/CD0 values
CL_map = s.LandingCheck.CL_map;
CD0_map = s.LandingCheck.CD0_map;


%% Loop over CL and CD0 values to generate map

% Initialize arrays
WS = linspace(s.WSmin,s.WSmax,s.n);
WS_map = NaN(length(CL_map),length(CD0_map));
WP_map = NaN(length(CL_map),length(CD0_map));
WS_CLmax = NaN(1,length(CD0_map));
WP_CLmax = NaN(1,length(CD0_map));
WS_CD0max = NaN(length(CL_map),1);
WP_CD0max = NaN(length(CL_map),1);
initial_CD0 = a.L.CD0;

% Loop over CD0 values
for j = 1:length(CD0_map)
    
    % Update CD0 input value
    a.L.CD0 = CD0_map(j);
    
    % Loop over WS values, each of which corresponds to a determined CL
    for i = 1:length(WS)
        
        % Wing loading in flight condition
        WS_in = WS(i)*m.L.f;
        
        % Initial guess to speed up convergence
        TW_in = q*a.L.CD0./WS_in + WS_in/pi/a.AR/a.L.e/q;
        
        % Compute thrust and power loading
        [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc] = ComputeThrustLoading_vKnown(...
            'L',TW_in,WS_in,climb,maneuver,dvdt,a,m,p,f,s,c);
        
        % Correct to MTOW
        TW(i) = TW_out*m.L.f;
        WP(i) = WP_out/m.L.f;
        
        % If the wing loading is unrealistically high and it is starting to
        % diverge, stop
        if TW(i) > 10e2
            TW(i) = NaN;
            WP(i) = NaN;
            break
        end
        
        % Save aerodynamic variables
        a.L.CL(i) = CL;
        a.L.dCL(i) = dCL;
        a.L.dCDi(i) = dCDi;
        p.L.Tc(i) = Tc;
        p.L.detap(i) = detap;
        
    end
    
    % For unrealistically high WS values the solution will not converge,
    % remove indexes where any of the results has NaN value
    NaNarray = sum([TW; WP; a.L.CL; a.L.dCL; a.L.dCDi; p.L.Tc],1);
    TW = TW(~isnan(NaNarray));
    WP = WP(~isnan(NaNarray));
    WS = WS(~isnan(NaNarray));
    a.L.CL = a.L.CL(~isnan(NaNarray));
    a.L.dCL = a.L.dCL(~isnan(NaNarray));
    a.L.dCDi = a.L.dCDi(~isnan(NaNarray));
    p.L.Tc = p.L.Tc(~isnan(NaNarray));
    
    % Interpolate CL values to common CL map, store in columns
    WS_map(:,j) = interp1(a.L.CL,WS,CL_map,'linear',NaN);
    WP_map(:,j) = interp1(a.L.CL,WP,CL_map,'linear',NaN);
    
    % Obtain WP/WS curve at CL max
    WS_CLmax(j) = interp1(CL_map,WS_map(:,j),a.L.CLmax,'linear');
    WP_CLmax(j) = interp1(CL_map,WP_map(:,j),a.L.CLmax,'linear');
end

% Obtain WP/WS curve at CD0 max
for i = 1:length(CL_map)
    WS_CD0max(i) = interp1(CD0_map,WS_map(i,:),initial_CD0,'linear');
    WP_CD0max(i) = interp1(CD0_map,WP_map(i,:),initial_CD0,'linear');
end


%% Add plot to WP-WS diagram

% Find propulsive WP-WS plot
index = 0;
figName = 'Total propulsive wing-loading/power-loading diagram';
for k = 1:size(s.figs,2)
    figHandle = s.figs(k);
    try
        if strcmp(figHandle.Name,figName)
            index = k;
            break
        end
    catch
        % Do nothing if an error is returned, which is normally due to
        % empty graphics holders in the graphics array (e.g. if a figure
        % has been closed)
    end
    
end

% If no figure was found, issue warning and do nothing else
if index == 0
    disp([s.levelString '  > No propulsive power loading plot was '...
        'found to show the landing drag requirements'])
else
    
    % Select figure
    figure(s.figs(index))
    
    % Plot lines of constant CL
    for i = 1:length(CL_map)
        plot(WS_map(i,:),WP_map(i,:),'-','color',[0.5 0.5 0.5])
        [WSlabel,idx] = min(WS_map(i,:));
        WPlabel = WP_map(i,idx);
        text(WSlabel,WPlabel,['C_L = ' num2str(CL_map(i),'%.2f') '  '],...
            'HorizontalAlignment','right','Rotation',-90)
    end
    
    % Plot lines of constant CD0
    for j = 1:length(CD0_map)
        plot(WS_map(:,j),WP_map(:,j),'-','color',[0.5 0.5 0.5])
        [WSlabel,idx] = max(WS_map(:,j));
        WPlabel = WP_map(idx,j);
        text(WSlabel,WPlabel,['  C_{D0} = ' num2str(CD0_map(j),'%.2f')])
    end
    
    % Plot limits
    plot(WS_CLmax,WP_CLmax,'-','color',[0.5 0.5 0.5],'linewidth',2)
    plot(WS_CD0max,WP_CD0max,'-','color',[0.5 0.5 0.5],'linewidth',2)
end


%% Verify throttle setting required to perform landing
% Note: this throttle is not equal to the throttle used in the input file
% during the landing constraint. That one simply determines the GT size
% necessary to perform this constraint at the given throttle setting. The
% one computed here assumes the design power loading has already been
% selected based on all the constraints evaluated previously. Since the
% power loading obtained from the other constraints is typically much more 
% restrictive than the one during the approach phase (which is not 
% explicitly shown in the WP-WS diagrams since it is assumed not be 
% limiting), a throttle value smaller than 1 will be required to perform
% the landing condition. 

% Interpolate limiting curves to common grid
% When no aero-propulsive effects are present, the constant CL line is
% vertical. In that case the intersection must be computed differently
aux = WS_CLmax(~isnan(WS_CLmax));
if all(aux == aux(1))
    WS_int = aux(1);
    WP_int = interp1(WS_CD0max(~isnan(WS_CD0max)),...
        WP_CD0max(~isnan(WS_CD0max)),WS_int,'linear',NaN);
else
    WP_CLmax_intp = interp1(WS_CLmax(~isnan(WS_CLmax)),...
        WP_CLmax(~isnan(WS_CLmax)),WS,'linear',NaN);
    WP_CD0max_intp = interp1(WS_CD0max(~isnan(WS_CD0max)),...
        WP_CD0max(~isnan(WS_CD0max)),WS,'linear',NaN);
    
    % Compute propulsive power loading in stall conditions
    for n = 2:length(WP_CLmax_intp)
        
        % Check when lines cross each other
        if (WP_CLmax_intp(n) >= WP_CD0max_intp(n) && ...
                WP_CLmax_intp(n-1) <= WP_CD0max_intp(n-1)) || ...
                (WP_CLmax_intp(n) <= WP_CD0max_intp(n) && ...
                WP_CLmax_intp(n-1) >= WP_CD0max_intp(n-1))
            D1 = WP_CLmax_intp(n-1) - WP_CD0max_intp(n-1);
            D2 = WP_CD0max_intp(n) - WP_CLmax_intp(n);
            D3 = WP_CD0max_intp(n) - WP_CD0max_intp(n-1);
            WP_int = WP_CD0max_intp(n-1) + D3*D1/(D1+D2);
            WS_int = WS(n-1) + (WS(n) - WS(n-1))*D1/(D1+D2);
            break
        end
    end
end

% Add point to plot
plot(WS_int,WP_int,'ok','markerfacecolor','k')

% Calculate detap in landing condition
TW_int = 1/WP_int/m.L.v;
D2W = f.D2W(WS_int,p.b_dp,p.N,p.dy,a.AR);
Tc = p.L.T*TW_int/p.N/m.L.rho/m.L.v^2/D2W;
Rc = 0.5*(D2W*WS_int*a.AR)^0.5;
[~,~,~,detap] = WingPropDeltas(Tc,Rc,p.xp,a.AR,a.L.CLmax,m.L.M,a.Lambda,...
                                            p.L.Gamma,p.L.etap,0);
    
% Collect powertrain component efficiencies
etas.GT = p.eta_GT;     % Gas turbine
etas.GB = p.eta_GB;     % Gearbox
etas.EM1 = p.eta_EM1;   % Primary electrical machine
etas.PM = p.eta_PM;     % PMAD (power management and distribution)
etas.EM2 = p.eta_EM2;   % Secondary electrical machine
if p.DP == 1
    etas.P1 = p.L.etap1 + detap;
    etas.P2 = p.L.etap2;
elseif p.DP == 2
    etas.P1 = p.L.etap1;
    etas.P2 = p.L.etap2 + detap;
else
    etas.P1 = p.L.etap1;
    etas.P2 = p.L.etap2;
end

% Call function to compute power required from GT (or EM)
PW_int = 1/WP_int;  % Not correcting for mass fraction here, since the 
                    %   available power computed below is also not
                    %   corrected, so the resulting throttle value is the
                    %   same.
[PW_out,~,~,~,~] = PowerTransmissionComputation_v2(p.config,etas,...
                    m.L.phi,m.L.Phi,[],PW_int,[],[],[]);
                
% Distinguish between powertrains with and without GT's
if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
    PW_req = max([abs(PW_out.gb) abs(PW_out.e1)]);
    PW_avail = 1/WPdes.(s.SelDes).EM1M;
elseif strcmp(p.config,'e-2')
    PW_req = max([abs(PW_out.e2) abs(PW_out.s2)]);
    PW_avail = 1/WPdes.(s.SelDes).EM2M;
else
    PW_req = PW_out.gt;
    PW_avail = 1/WPdes.(s.SelDes).GTM;
end

% throttle required for equilibrium flight point
xi_int = PW_req/PW_avail;

% Add label to plot
text(WS_int,WP_int,['  \xi = ' num2str(xi_int,'%.2f')],...
        'VerticalAlignment','Bottom');

% Add warning if the powers on one of the other branches have been exceeded
if PW_out.gt > 1/WPdes.(s.SelDes).GTM
   disp([s.levelString '    > Warning: the design power loading of '...
            'the gas turbine is exceeded during '...
            'landing conditions'])
end
if max([abs(PW_out.gb) abs(PW_out.e1)]) > 1/WPdes.(s.SelDes).EM1M
    disp([s.levelString '    > Warning: the design power loading of '...
            'the primary electrical machine is exceeded during '...
            'landing conditions'])
end
if max([abs(PW_out.e2) abs(PW_out.s2)]) > 1/WPdes.(s.SelDes).EM2M
    disp([s.levelString '    > Warning: the design power loading of '...
            'the secondary electrical machine is exceeded during '...
            'landing conditions'])
end
if PW_out.bat > 1/WPdes.(s.SelDes).bat
    disp([s.levelString '    > Warning: the design power loading of '...
            'the batteries is exceeded during '...
            'landing conditions'])
end
if PW_out.s1 > 1/WPdes.(s.SelDes).s1
    disp([s.levelString '    > Warning: the design power loading of '...
            'the primary shaft is exceeded during '...
            'landing conditions'])
end


%% Evaluate additional thrust shares/settings
 %{
% chi-values (i.e. T-vaues)
chi_array = 0:0.25:1;

% Make cell arrays; first cell corresonds to previous results, i.e. at the
% specified design T
cWS_map = cell(1,length(chi_array));
cWP_map = cell(1,length(chi_array));
cWS_CLmax = cell(1,length(chi_array));
cWP_CLmax = cell(1,length(chi_array));
cWS_CD0max = cell(1,length(chi_array));
cWP_CD0max = cell(1,length(chi_array));
cWS_map{1} = WS_map;
cWP_map{1} = WP_map;
cWS_CLmax{1} = WS_CLmax;
cWP_CLmax{1} = WP_CLmax;
cWS_CD0max{1} = WS_CD0max;
cWP_CD0max{1} = WP_CD0max;
chi_array = [p.L.T chi_array];

% Loop over cells
for k = 2:length(chi_array)
    
    % Update power-share value
    p.L.T = chi_array(k);
    
    % Initialize arrays
    WS = linspace(s.WSmin,s.WSmax,s.n);
    WS_map = NaN(length(CL_map),length(CD0_map));
    WP_map = NaN(length(CL_map),length(CD0_map));
    WS_CLmax = NaN(1,length(CD0_map));
    WP_CLmax = NaN(1,length(CD0_map));
    WS_CD0max = NaN(length(CL_map),1);
    WP_CD0max = NaN(length(CL_map),1);
    
    % Loop over CD0 values
    for j = 1:length(CD0_map)
        
        % Update CD0 input value
        a.L.CD0 = CD0_map(j);
        
        % Loop over WS values, each of which corresponds to a determined CL
        for i = 1:length(WS)
            
            % Wing loading in flight condition
            WS_in = WS(i)*m.L.f;
            
            % Initial guess to speed up convergence
            TW_in = q*a.L.CD0./WS_in + WS_in/pi/a.AR/a.L.e/q;
            
            % Compute thrust and power loading
            [WP_out,TW_out,CL,dCL,dCDi,~,detap,Tc] = ComputeThrustLoading_vKnown(...
                'L',TW_in,WS_in,climb,maneuver,dvdt,a,m,p,f,s,c);
            
            % Correct to MTOW
            TW(i) = TW_out*m.L.f;
            WP(i) = WP_out/m.L.f;
            
            % If the wing loading is unrealistically high and it is starting to
            % diverge, stop
            if TW(i) > 10e2
                TW(i) = NaN;
                WP(i) = NaN;
                break
            end
            
            % Save aerodynamic variables
            a.L.CL(i) = CL;
            a.L.dCL(i) = dCL;
            a.L.dCDi(i) = dCDi;
            p.L.Tc(i) = Tc;
            p.L.detap(i) = detap;
            
        end
        
        % For unrealistically high WS values the solution will not converge,
        % remove indexes where any of the results has NaN value
        NaNarray = sum([TW; WP; a.L.CL; a.L.dCL; a.L.dCDi; p.L.Tc],1);
        TW = TW(~isnan(NaNarray));
        WP = WP(~isnan(NaNarray));
        WS = WS(~isnan(NaNarray));
        a.L.CL = a.L.CL(~isnan(NaNarray));
        a.L.dCL = a.L.dCL(~isnan(NaNarray));
        a.L.dCDi = a.L.dCDi(~isnan(NaNarray));
        p.L.Tc = p.L.Tc(~isnan(NaNarray));
        
        % Interpolate CL values to common CL map, store in columns
        WS_map(:,j) = interp1(a.L.CL,WS,CL_map,'linear',NaN);
        WP_map(:,j) = interp1(a.L.CL,WP,CL_map,'linear',NaN);
        
        % Obtain WP/WS curve at CL max
        WS_CLmax(j) = interp1(CL_map,WS_map(:,j),a.L.CLmax,'linear');
        WP_CLmax(j) = interp1(CL_map,WP_map(:,j),a.L.CLmax,'linear');
    end
    
    % Obtain WP/WS curve at CD0 max
    for i = 1:length(CL_map)
        WS_CD0max(i) = interp1(CD0_map,WS_map(i,:),initial_CD0,'linear');
        WP_CD0max(i) = interp1(CD0_map,WP_map(i,:),initial_CD0,'linear');
    end
    
    % Save to arrays
    cWS_map{k} = WS_map;
    cWP_map{k} = WP_map;
    cWS_CLmax{k} = WS_CLmax;
    cWP_CLmax{k} = WP_CLmax;
    cWS_CD0max{k} = WS_CD0max;
    cWP_CD0max{k} = WP_CD0max;
end


%% Generate 3D plot

% Select figure
figure(300)
hold on; grid on

% Color schemes
cols = [0.5 0.5 0.5
        1 0 0 
        0 1 0
        0 0 1
        1 0 1
        0 1 1];

% Loop over T-values
for k = 1:length(chi_array)
    
    % Plot lines of constant CL
    for i = 1:length(CL_map)
        plot3(cWS_map{k}(i,:),chi_array(k)*ones(size(cWS_map{k}(i,:))),cWP_map{k}(i,:),'-','color',cols(k,:))
        [WSlabel,idx] = min(cWS_map{k}(i,:));
        WPlabel = cWP_map{k}(i,idx);
        text(WSlabel,chi_array(k),WPlabel,['C_L = ' num2str(CL_map(i),'%.2f') '  '],...
            'HorizontalAlignment','right','Rotation',-90)
    end
    
    % Plot lines of constant CD0
    for j = 1:length(CD0_map)
        plot3(cWS_map{k}(:,j),chi_array(k)*ones(size(cWS_map{k}(:,j))),cWP_map{k}(:,j),'-','color',cols(k,:))
        [WSlabel,idx] = max(cWS_map{k}(:,j));
        WPlabel = cWP_map{k}(idx,j);
        text(WSlabel,chi_array(k),WPlabel,['  C_{D0} = ' num2str(CD0_map(j),'%.2f')])
    end
    
    % Plot limits
    plot3(cWS_CLmax{k},chi_array(k)*ones(size(cWS_CLmax{k})),cWP_CLmax{k},'-','color','k','linewidth',2)
    plot3(cWS_CD0max{k},chi_array(k)*ones(size(cWS_CD0max{k})),cWP_CD0max{k},'-','color','k','linewidth',2)
end
xlabel('W/S')
ylabel('\chi')
zlabel('W/P')
ax=gca;
ax.View = [0 0];
%}




