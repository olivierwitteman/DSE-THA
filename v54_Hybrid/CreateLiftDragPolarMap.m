function [Tc_array,M_array,CLiso_opt,LD_opt,fig] = CreateLiftDragPolarMap(WS,con,a,p,f,s,options)
%%% Description
% This function computes the total installed lift and drag coefficients for
% a range of mach numbers and thrust coefficient values, in order to obtain
% the L/D maps, based on which the optimum CL_iso for maximum L/D can be
% selected. The function is meant to be used after the power-loading
% diagrams, when WS has been determined, for a given flight condition.
%
% Input
%   - WS: design wing loading [N/m2] (NOT corrected to current aircraft
%       weight!)
%   - con: active flight condition (e.g. 'cr','cl'...)
%   - a,p,f,s: structures containing aicraft, propulsion system, anonymous
%       functions and program settings data (see WP_WS_diagram_DP.m)
%   - options: plot figure if equal to 1
%
% Output
%   - M_array: vector containing M values sampled
%   - Tc_array: vector containing Tc values sampled
%   - CL_iso_opt: 2D array (first dimension = Tc values, second dimension =
%       M values) containing the airframe lift coefficient corresponding
%       to maximum aircraft L/D (incl. aero-propulsive effects) for each
%       combination of M and Tc.
%   - LD_opt: 2D array (first dimension = Tc values, second dimension =
%       M values) containing the maximum lift-to-drag ratio (incl. 
%       aero-propulsive effects) for each combination of M and Tc
%   - fig: figure handle
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 21-03-18
%%% Last modified: 22-03-18


%% Pre-calculations

% Disk loading [m2/N]
D2W = f.D2W(WS,p.b_dp,p.N,p.dy,a.AR);

% Compute prop radius/wing chord ratio
Rc = 0.5*(D2W*WS*a.AR)^0.5;

% Initialize arrays
Tc_array = linspace(s.Polar.TcInterval(1),s.Polar.TcInterval(2),s.Polar.N);
M_array = linspace(s.Polar.MInterval(1),s.Polar.MInterval(2),s.Polar.N);
CLiso_array = linspace(s.Polar.CLisoInterval(1),...
                        s.Polar.CLisoInterval(2),s.Polar.N_CLiso);
CL_array = NaN(s.Polar.N,s.Polar.N,s.Polar.N_CLiso);
CD_array = NaN(s.Polar.N,s.Polar.N,s.Polar.N_CLiso);
LD_array = NaN(s.Polar.N,s.Polar.N,s.Polar.N_CLiso);
LDiso_array = NaN(s.Polar.N,s.Polar.N,s.Polar.N_CLiso);
LD_opt = NaN(s.Polar.N,s.Polar.N);
CLiso_opt = NaN(s.Polar.N,s.Polar.N);
LDiso_opt = NaN(s.Polar.N,s.Polar.N);


%%  Generate polar maps

% Loop over: 1st dimension = Tc, 2nd dimension = M, 3rd dimension = CLiso
for i = 1:s.Polar.N
    for j = 1:s.Polar.N
        for k = 1:s.Polar.N_CLiso

            % Compute deltas
            [dCl,dCd0,dCdi,~] = WingPropDeltas(Tc_array(i),Rc,...
                p.xp,a.AR,CLiso_array(k),...
                M_array(j),a.Lambda,p.(con).Gamma,p.(con).etap,0);
            dCL = dCl*p.b_dp;
            dCD0 = dCd0*p.b_dp;
            dCDi = dCdi*p.b_dp;
            
            % Lift and drag coefficients
            CL_array(i,j,k) = CLiso_array(k) + dCL;
            CD_array(i,j,k) = a.(con).CD0 + CLiso_array(k)^2/pi...
                                /a.AR/a.(con).e + dCD0 + dCDi;
            LD_array(i,j,k) = CL_array(i,j,k)/CD_array(i,j,k);
            LDiso_array(i,j,k) = CLiso_array(k)/(a.(con).CD0 + ...
                                 CLiso_array(k)^2/pi/a.AR/a.(con).e);
        end
        
        % Identify and save optimum
        [LD_opt(i,j),idx] = max(LD_array(i,j,:));
        CLiso_opt(i,j) = CLiso_array(idx);
        LDiso_opt(i,j) = LDiso_array(i,j,idx);
    end
end


%% Plotting
if options == 1
    
    % Generate figure
    fig = figure(s.figStart + size(s.figs,2));
    fig.Name = ['Aerodynamic polar maps (' con ')'];
    
    % CLiso-CL versus Tc (Take middle of interval as sample M)
    idxM = round(s.Polar.N/2);
    subplot(2,4,1)
    hold on; grid on; box on;
    surf(CLiso_array,Tc_array,squeeze(CL_array(:,idxM,:)),'edgecolor','none')
    xlabel('C_{L,iso}'); ylabel('T_{c}'); zlabel('C_{L}')
    title(['M = ' num2str(M_array(idxM),'%5.3f')])
    
    % CLiso-CD versus Tc
    subplot(2,4,2)
    hold on; grid on; box on;
    surf(CLiso_array,Tc_array,squeeze(CD_array(:,idxM,:)),'edgecolor','none')
    xlabel('C_{L,iso}'); ylabel('T_{c}'); zlabel('C_{D}')
    title(['M = ' num2str(M_array(idxM),'%5.3f')])
    
    % CLiso-LD versus Tc
    subplot(2,4,3)
    hold on; grid on; box on;
    surf(CLiso_array,Tc_array,squeeze(LD_array(:,idxM,:)),...
        'edgecolor','none');
    surf(CLiso_array,Tc_array,squeeze(LDiso_array(:,idxM,:)),...
        'facecolor','none','edgecolor',[0.5 0.5 0.5]);
    xlabel('C_{L,iso}'); ylabel('T_{c}'); zlabel('L/D')
    title(['M = ' num2str(M_array(idxM),'%5.3f')])
    
    % CL_iso corresponding to maximum L/D versus Tc and M
    subplot(2,4,4)
    hold on; grid on; box on;
    surf(Tc_array,M_array,CLiso_opt)
    xlabel('T_{c}'); ylabel('M'); zlabel('C_{L,iso,opt}')
    
    % CLiso-CL versus Mach number
    idxTc = round(s.Polar.N/2);
    subplot(2,4,5)
    hold on; grid on; box on;
    surf(CLiso_array,M_array,squeeze(CL_array(idxTc,:,:)),'edgecolor','none')
    xlabel('C_{L,iso}'); ylabel('M'); zlabel('C_{L}')
    title(['T_c = ' num2str(Tc_array(idxTc),'%5.3f')])
    
    % CLiso-CD versus Mach number
    subplot(2,4,6)
    hold on; grid on; box on;
    surf(CLiso_array,M_array,squeeze(CD_array(idxTc,:,:)),'edgecolor','none')
    xlabel('C_{L,iso}'); ylabel('M'); zlabel('C_{D}')
    title(['T_c = ' num2str(Tc_array(idxTc),'%5.3f')])
    
    % CLiso-LD versus Mach number
    subplot(2,4,7)
    hold on; grid on; box on;
    surf(CLiso_array,M_array,squeeze(LD_array(idxTc,:,:)),'edgecolor','none')
    surf(CLiso_array,M_array,squeeze(LDiso_array(idxTc,:,:)),...
        'facecolor','none','edgecolor',[0.5 0.5 0.5]);
    xlabel('C_{L,iso}'); ylabel('M'); zlabel('L/D')
    title(['T_c = ' num2str(Tc_array(idxTc),'%5.3f')])
    
    % maximum L/D versus Tc and M
    subplot(2,4,8)
    hold on; grid on; box on;
    surf(Tc_array,M_array,LD_opt)
    surf(Tc_array,M_array,LDiso_opt,...
        'facecolor','none','edgecolor',[0.5 0.5 0.5]);
    xlabel('T_c'); ylabel('M'); zlabel('L/D_{max}')
    
end





