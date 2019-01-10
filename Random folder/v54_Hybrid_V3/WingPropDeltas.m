function [dCl,dCd0,dCdi,detap] = WingPropDeltas(Tc,Rc,xc,AR,Cl,M,Lambda,Gamma,etap,options)
% This function computes the lift, zero-lift drag, and thrust-induced drag
% increase of a wing featuring wing-mounted distributed propulsions
% systems, as well as providing the thrust coefficient in order to verify
% whether the values are realistic or not.
%
% Input:
%   - Tc: thrust coefficient of a single propulsor,
%         defined as Tc = T/(rho*v^2*D*2) [-]
%   - Rc: Ratio between propulsor radius and wing chord [-]
%   - xc: axial propeller position as a function of chord:
%           xc < 0: tractor
%           0 < xc < 1: OTW
%           xc > 1: pusher
%           xc = Inf: No effect of propeller on wing
%   - AR: wing aspect ratio 
%   - Cl: estimation of isolated wing sectional lift coefficient [-]
%   - M: Mach number [-]
%   - Lambda: half-chord wing sweep angle [deg]
%   - Gamma: propeller angle of attack w.r.t. freestream [deg]
%   - etap: isolated propeller propulsive efficiency [-]
%   - options: shows plots if set to 1
%
% Output:
%   - dCl: sectional lift coefficient increased due to DP w.r.t. clean
%       wing. 
%   - dCd0: zero-lift/zero-thrust drag penalty due to e.g. pylons and
%       nacelles of DP system.
%   - dCdi: thrust-induced drag, i.e. sectional drag coefficient increase
%       due to the DP system in power-on conditions.
%   - detap: change in propulsive efficiency, etap_installed - etap_iso.
%
% All coefficients are two-dimensional and supposed constant throughout the
% DP system span interval (so not along the full wing). Correlations are
% based on empirical/semi-analytical methods. The propulsors are assumed
% far away from the wing tips.
%
% This version of the code assumes zero dCd0 for all configurations; see
% end of code.
%
%
%%% Reynard de Vries
%%% TU Delft
%%% Date created: 21-06-17
%%% Last modified: 06-02-18


%% OTW configurations
if xc > 0 && xc < 1
  
    % Values obtained from ObtainExperimentalResultsOTW_dCl_dCd.m function:
    % inserted polynomial coefficients here to improve speed
    pCl =    [6.143649865484650   2.457028957705863  -0.000000000000000
              3.384339343604002   1.371296849483363   0.000000000000000
              1.875306713643190   0.897026995137536  -0.000000000000000
              1.441958159774097   0.762397071062925   0.000000000000000
              1.365352092825532   0.984741283870551   0.000000000000000];
    
    pCd2 =    10000*[1.043872382670962   0.002902934229647   0.000000000000000
                     0.030317753144338   0.001197146779598  -0.000000000000000
                     0.030902619622957   0.001129879069199   0.000000000000000
                     0.172431418321773   0.002460827952950  -0.000000000000000
                     0.000000691465193   0.000000142894078   0.000000000000000];
    x = [0 0.2524 0.5019 0.758 1];
    
    % Limit the parabolas to a maximum Tc, afterwards suppose constant, to avoid divergence
    Tmax = 1;
    if Tc > Tmax 
        TC = Tmax;
    else
        TC = Tc;
    end
    for i = 1:length(x)
        dCl4(i,:) = (-pCl(i,2)+(pCl(i,2)^2-4*pCl(i,1)*(pCl(i,3)-TC)).^0.5)/(2*pCl(i,1));
        if i ==length(x)
            dCd4(i,:) = pCd2(i,1)*TC.^2 + pCd2(i,2)*TC + pCd2(i,3);
        else
            dCd4(i,:) = -(-pCd2(i,2)+(pCd2(i,2)^2-4*pCd2(i,1)*(pCd2(i,3)-TC)).^0.5)/(2*pCd2(i,1));
        end
    end
    dCl = interp1(x,dCl4,xc,'linear');
    dCdi = interp1(x,dCd4,xc,'linear');
    
    % 3. Sectional zero-lift drag coefficient increase (generic value)
    dCd0 = 0.000;    
    
    % 4. Propulsive efficiency: assume GENERIC 5% reduction in propulsive
    % efficiency, to demonstrate sensitivity of method. To be replaced.
    detap = -0.05*etap;
    
    % Plot
    if options == 1
        
        % Generate map for plotting
        Tc_map = linspace(0,Tmax,10);
        xc_map = x;
        dCl_map =zeros(length(xc_map),length(Tc_map));
        dCd_map =zeros(length(xc_map),length(Tc_map));
        for j = 1:length(Tc_map)
            TC = Tc_map(j);
            for i = 1:length(xc_map)
                dCl_map(i,j) = (-pCl(i,2)+(pCl(i,2)^2-4*pCl(i,1)*(pCl(i,3)-TC)).^0.5)/(2*pCl(i,1));
                if i ==length(xc_map)
                    dCd_map(i,j) = pCd2(i,1)*TC.^2 + pCd2(i,2)*TC + pCd2(i,3);
                else
                    dCd_map(i,j) = -(-pCd2(i,2)+(pCd2(i,2)^2-4*pCd2(i,1)*(pCd2(i,3)-TC)).^0.5)/(2*pCd2(i,1));
                end
            end
        end
        
        % Include original experimental data points
        xc_exp = [0.0006    0.2524    0.5019    0.7580    1.0038]';
        Tc_exp = [0    0.1700    0.5300];
        dCl_exp =  [0    0.0601    0.1554
                    0    0.0995    0.2420
                    0    0.1453    0.3438
                    0    0.1690    0.3970
                    0    0.1439    0.3593];
        dCd_exp =  [0   -0.0029   -0.0059
                    0   -0.0111   -0.0265
                    0   -0.0115   -0.0270
                    0   -0.0051   -0.0118
                    0    0.0004    0.0027];
        
        % Generate figure (surfaces)       
        figure('name','Lift/Drag coefficient increase, OTW configuration')
        subplot(1,2,1); grid on; hold on; box on;
        surf(xc_map,Tc_map,dCl_map')
        surf(xc_exp,Tc_exp,dCl_exp','facecolor','none','edgecolor','none','marker','o','markeredgecolor','k','markerfacecolor','k')
        scatter3(xc,Tc,dCl,36,'markerfacecolor','r')
        xlabel('Axial propeller position')
        ylabel('Thrust coefficient')
        zlabel('{\Delta}C_l')
        subplot(1,2,2); grid on; hold on; box on;
        surf(xc_exp,Tc_exp,dCd_exp','facecolor','none','edgecolor','none','marker','o','markeredgecolor','k','markerfacecolor','k')
        surf(xc_map,Tc_map,dCd_map')
        scatter3(xc,Tc,dCdi,36,'markerfacecolor','r')
        xlabel('Axial propeller position')
        ylabel('Thrust coefficient')
        zlabel('{\Delta}C_d')

%         % Generate figure (contours, for saving)
%         xc_exp = [0    0.2524    0.5019    0.7580    1]';
%         c = colormap(gray);
%         cc2 = c(end-30:end,:);
%         cc = flip(cc2,1);
%         fig = figure('name','Lift/Drag coefficient increase, OTW configuration','color',[1 1 1]);
%         subplot(1,2,1); hold on; box on;
%         contourf(xc_map+1e-6,Tc_map+1e-6,dCl_map','LevelList',-0.05:0.005:0.01)
%         [ccc,hhh] = contourf(xc_map,Tc_map,dCl_map','ShowText','on','LevelList',0:0.05:0.6,...
%             'TextList',0:0.1:0.6,'LabelSpacing',120);
%         clabel(ccc,hhh,'FontName','Times')
%         ax = gca;ax.Layer = 'top';
%         ax.Position = [0.09 0.27 0.38 0.72];
%         ax.XTick = 0:0.2:1;
%         ax.YTick = 0:0.2:1;
%         axis([0 1 0 1])
%         ax.XAxis.TickLabelFormat = '%.1f';
%         ax.YAxis.TickLabelFormat = '%.1f';
%         [xx,yy] = meshgrid(xc_exp,Tc_exp);
%         scatter(xx(:),yy(:),25,'markeredgecolor','k','markerfacecolor','k');
%         xlabel('Axial propeller position \it{x_{\rm{p}\it}} \rm[-]')
%         ylabel('Thrust coefficient \it{T_{\rm{c}\it}} \rm[-]')
%         tt = title('a) Lift coefficient increase \it{\Delta}{C_l} \rm[-]');
%         tt.Position = [0.5 -0.37 0];
%         ax.FontName = 'Times';
%         colormap(ax,cc)
%         subplot(1,2,2); hold on; box on;
%         contourf(xc_map+1e-6,Tc_map+1e-6,dCd_map','LevelList',-0.05:0.005:0.01)
%         [ccc,hhh] = contourf(xc_map,Tc_map,dCd_map','ShowText','on','LevelList',-0.05:0.005:0.01,...
%             'TextList',-0.05:0.01:0.01,'LabelSpacing',140);
%         clabel(ccc,hhh,'FontName','Times')
%         ax2 = gca; ax2.Layer = 'top';
%         ax2.Position = [0.6 0.27 0.38 0.72];
%         ax2.XTick = 0:0.2:1;
%         ax2.YTick = 0:0.2:1;
%         axis([0 1 0 1])
%         ax2.XAxis.TickLabelFormat = '%.1f';
%         ax2.YAxis.TickLabelFormat = '%.1f';
%         scatter(xx(:),yy(:),25,'markeredgecolor','k','markerfacecolor','k')
%         xlabel('Axial propeller position \it{x_{\rm{p}\it}} \rm[-]')
%         tt = title('b) Induced drag coeff. increase \it{\Delta}{C_{d,\rm{i}}} \rm[-]');
%         tt.Position = [0.5 -0.37 0];
%         ax2.FontName = 'Times';
%         fig.PaperUnits = 'inches';
%         fig.PaperSize = [6.1 3.1];
%         fig.PaperPosition = [0 0 6 3];
%         colormap(ax2,cc2)
%         % print('OTWDeltasMap.pdf','-dpdf')
    end

       
%% Tractor configurations
elseif xc < 0
    
    % 1. Sectional lift coefficient increase
    %
    % This method is based on an AD + flat plate assumption, described in
    % Patterson & German (2015), Patterson et al. (2016) and Patterson's 
    % PhD disseration (2016). An empirical correction/surrogate model
    % for finite slipstream height based on a surrogate model is applied.
    %
    % Assumptions/Notes
    %   - AD with AOA = 0 assumed for velocity increase at propeller disk and
    %       slipstream evolution. No nacelle included for slipstream evolution
    %   - Wing lift curve slope approximated by 2*pi: flat plate. The resulting
    %       Delta Cl value is not excessively sensitive to wing lift curve
    %       slope.
    %   - Effect of propeller on wing loading limited to slipstream (no effect
    %       on rest of wing span.
    %   - Wing is supposed fully immersed in slipstream (no effect of propeller
    %       height considered)
    %   - AD: effect of swirl is neglected
    %
    % Input values for comparison with Leo 2005 (p. 105) PROWIM experiment
    % Cl = 0.644;
    % v_inf = 49.5;
    % Tc = 0.2;
    % Gamma = 0;
    % rho = 1.225;
    % D = 0.236;
    % AR = 2*0.64/0.24;
    % xLE = 0.84;
    % N = 2;
    % c = 0.24;
    % v_p = 0.5*v_inf*((1+8/pi*Tc)^0.5-1);
    % hc = pi/4*kLE*D/2/c;
    
   
    % Velocity INCREASE at propeller disk [m/s] (v_effective = v_inf + v_p)
    % (Veldhuis 2005, Appendix A)
    a_p = 0.5*(sqrt(1+8/pi*Tc)-1);
    
    % From Patteson PhD 2016 (p. 93/94). 
    % Surrogate model is valid for SS velocity ratios of 1.25-2.25, i.e.
    % for v_p/v_inf values of 0.125-0.625, R/c values of 0.125-3, and xc
    % positions ranging from -0.25 to -3. Note: if the thrust coefficient
    % is too low, a < 0.125 may be obtained, which could lead to
    % unrealistically high beta values. In this case beta is limited to 1.
   
    % Axial position of propeller as a fraction of propeller radius, POSITIVE
    % if propeller is ahead of QUARTER CHORD. 
    xr = -(xc-1/4)/Rc;
    
    % Fraction between slipstream diameter and propeller diameter at wing 
    % quarter chord (Veldhuis 2005, Appendix C)
    kSS = @(a,x) ((1+a)./(1+a.*(1+x./(x.^2+1).^0.5))).^0.5;
    k_c4 = kSS(a_p,xr);
    
    % Conservation of mass: velocity at quarter chord [m/s]
    a_c4 = (1+a_p)/k_c4^2-1;

    % Column vector X. Uses the velocity ratio of the jet far downstream,
    % and according to AD theory a_farDownstream = 2*a_p. The model also
    % includes the position of the
    % propeller, supposedly to account for contraction of the slipstream
    % radius.
    X = [1 -xc (-xc)^2 -xc*(2*a_p+1) (2*a_p+1) (2*a_p+1)^2]';
    
    % Constants from surrogate model (Patterson, 2016)
    K0 = [0.378269 0.748135 -0.179986 -0.056464 -0.146746 -0.015255];
    K1 = [3.071020 -1.769885 0.436595 0.148643 -0.989332 0.197940];
    K2 = [-2.827730 2.054064 -0.467410 -0.277325 0.698981 -0.008226];
    K3 = [0.997936 -0.916118 0.199829 0.157810 -0.143368 -0.057385];
    K4 = [-0.127645 0.135543 -0.028919 -0.026546 0.010470 0.012221];

    % Beta factor
    Beta = K0*X + K1*X*Rc + K2*X*Rc^2 + K3*X*Rc^3 + K4*X*Rc^4;
    
    % Due to extrapolation/surrogate model, beta might be slightly larger 
    % than 1. Limit to 1.
    Beta = min([Beta 1]);
    
    % Angle between propeller axis and wing chord [deg]
    % Method for finite wings from Roskam: Methods dor estimation drag
    % polars of subsonic airplanes, S. 3.3.1. Assuming 2D lift curve slope
    % of airfoil sections = 2*pi.
    Mcorr = (1-M^2)^0.5;
    CLalpha = 2*pi*AR/(2+(AR^2*Mcorr^2*(1+(tand(Lambda))^2/Mcorr^2)+4)^0.5);
    alfa = rad2deg(Cl/CLalpha);
    ip = Gamma-alfa;
    
    % Sectional lift coefficient increase (Patterson 2015)
    dCl = 2*pi*((sind(alfa)-Beta*a_c4*sind(ip))*...
        (Beta^2*a_c4^2+2*Beta*a_c4*cosd(alfa+ip)+1)^0.5 - sind(alfa));

    % Plot
    if options == 1
        figure('name','Lift coefficient increase, tractor configuration')
        
        % Streamtube contraction
        subplot(1,2,1)
        xp = -5:0.05:5;
        k_c4_plot = kSS(a_p,xp);
        ac4_plot = a_p./k_c4_plot;
        grid on; hold on; box on
        plot(xp,k_c4_plot,'b')
        plot(xp,ac4_plot,'r')
        for aa = [0 0.2 0.4 0.6 0.8 1]
            k_c4_plot = kSS(aa,xp);
            ac4_plot = (1+aa)./k_c4_plot.^2-1;
            plot(xp,k_c4_plot,'--')
            plot(xp,ac4_plot,'-.')
        end
        xlabel('x/R')
        ylabel('Velocity & Area ratio')
        legend('Slipstream radius','V/V_\infty-1',...
            'Slipstream radius, a = 0.0','V/V_\infty-1, a = 0.0',...
            'Slipstream radius, a = 0.2','V/V_\infty-1, a = 0.2',...
            'Slipstream radius, a = 0.4','V/V_\infty-1, a = 0.4',...
            'Slipstream radius, a = 0.6','V/V_\infty-1, a = 0.6',...
            'Slipstream radius, a = 0.8','V/V_\infty-1, a = 0.8',...
            'Slipstream radius, a = 1.0','V/V_\infty-1, a = 1.0')
        plot(xr,k_c4,'ob')
        plot(xr,a_c4,'or')
        axis([min(xp) max(xp) 0 3])
        title(['x/R = ' num2str(xr,'%.3f') ', a = ' num2str(a_p,'%.3f') ...
            ', a_{p} = ' num2str(a_p,'%.3f') ', r_{c/4}/R = ' ...
            num2str(k_c4,'%.3f') ', R/c = ' num2str(Rc,'%.3f') ...
            ', T_c = ' num2str(Tc,'%.3f')])
        
        % Finite slipstream height correction
        subplot(1,2,2)
        rc_array = 0.125:0.125:3;
        a_array = 0.125:0.02:0.625;
        beta_array = zeros(length(rc_array),length(a_array));
        for i = 1:length(rc_array)
            for j = 1:length(a_array)
                X = [1 -xc (-xc)^2 -xc*(2*a_array(j)+1) (2*a_array(j)+1) (2*a_array(j)+1)^2]';
                beta_array(i,j) = K0*X + K1*X*rc_array(i) + K2*X*rc_array(i)^2 + K3*X*rc_array(i)^3 + K4*X*rc_array(i)^4;
                beta_array(i,j) = min([beta_array(i,j) 1]);
            end
        end
        hold on; box on; grid on
        h1 = surf(rc_array,a_array,beta_array');
        rc_array2 = 0:0.125:3.5;
        a_array2 = 0.005:0.02:0.705;
        beta_array2 = zeros(length(rc_array2),length(a_array2));
        for i = 1:length(rc_array2)
            for j = 1:length(a_array2)
                X = [1 -xc (-xc)^2 -xc*(2*a_array2(j)+1) (2*a_array2(j)+1) (2*a_array2(j)+1)^2]';
                beta_array2(i,j) = K0*X + K1*X*rc_array2(i) + K2*X*rc_array2(i)^2 + K3*X*rc_array2(i)^3 + K4*X*rc_array2(i)^4;
                beta_array2(i,j) = min([beta_array2(i,j) 1]);
            end
        end
        hold on; box on; grid on
        surf(rc_array2,a_array2,beta_array2','facecolor','none')
        h2 = scatter3(Rc,a_p,Beta,'markerfacecolor','r','markeredgecolor','r');
        xlabel('R/c')
        ylabel('a = v_p/v_\infty')
        zlabel('\beta = {\Delta}c_l/{\Delta}c_{l,h=\infty}')
        title(['{\Delta}C_l = ' num2str(dCl,'%.3f') ' (' num2str(100*dCl/Cl,'%.3f') '%)'])
        legend([h1,h2],'Surrogate model limits','Design point')
        
    end

    % 2. Sectional drag coefficient increase, following Biber (2011). 
    % The increase in drag due to propeller thrust has two components: a
    % change in friction drag due to increased dynamic pressure over the
    % wing surface covered by the slipstream, and and increase in
    % lift-induced drag as a consequence of the change in lift due to the
    % propeller slipstream
    Cf = 0.009; % Local skin friction coefficient (see Biber 2011) 
    dCdi_0 = Cf*a_c4^2;
    dCdi_l = (2*Cl*dCl)/pi/AR;
    dCdi = dCdi_0 + dCdi_l;
    
    % 3. Sectional zero-lift drag coefficient increase (generic value)
    dCd0 = 0.000;
    
    % 4. Propulsive efficiency not affected
    detap = 0;
    
    
%% Pusher configurations
elseif xc > 1 && xc < Inf
     
    % Velocity INCREASE at propeller disk [m/s] (v_effective = v_inf + v_p)
    % (Veldhuis 2005, Appendix A)
    a_p = 0.5*(sqrt(1+8/pi*Tc)-1);
    
    % Second method: from Patteson PhD 2016 (p. 93/94). 
    % Surrogate model is valid for v_p/v_inf values of 0.125-0.625, 
    % R/c values of 0.125-3, and xc
    % positions ranging from -0.25 to -3. Note: if the thrust coefficient
    % is too low, a < 0.125 may be obtained, which could lead to
    % unrealistically high beta values. In this case beta is limited to 1.
    % Since the surrogate model is not valid for pusher propellers, it is
    % slightle tweaked: beta is evaluated using a hypothetical propeller
    % located at xc = -0.25, of v_p = v_pusher_-0.125 and R =
    % R_ss_pusher_-0.125.

    % Prop radius over wing chord ratio
    
    % Axial position of propeller as a fraction of propeller radius, POSITIVE
    % if propeller is ahead of QUARTER CHORD. 
    xr = -(xc-1/4)/Rc;
    
    % Fraction between slipstream diameter and propeller diameter at wing 
    % quarter chord (Veldhuis 2005, Appendix C)
    kSS = @(a,x) ((1+a)./(1+a.*(1+x./(x.^2+1).^0.5))).^0.5;
    k_c4 = kSS(a_p,xr);
    
    % Conservation of mass: velocity at quarter chord [m/s]
    a_c4 = (1+a_p)/k_c4^2-1;

    % Induction and radius of hypothetical propeller at xc = -0.125
    xc_125 = -0.125;
    xr_125 = -(xc_125-1/4)/Rc;
    k_125 = kSS(a_p,xr_125);
    a_125 = a_p/k_125^2;
    rc_125 = Rc*k_125;
    
    % Column vector X. Uses the velocity ratio of the jet (i.e. the
    % velocity far downstream) but also includes the position of the
    % propeller, supposedly to account for contraction of the slipstream
    % radius. Here, an equivalent propeller (of radius = radius of the SS
    % at the LE) is considered to be positioned at the x/c = -0.125.
    X = [1 -xc_125 (-xc_125)^2 -xc_125*(a_125+1) (a_125+1) (a_125+1)^2]';
    
    % Constants from surrogate model (Patterson, 2016)
    K0 = [0.378269 0.748135 -0.179986 -0.056464 -0.146746 -0.015255];
    K1 = [3.071020 -1.769885 0.436595 0.148643 -0.989332 0.197940];
    K2 = [-2.827730 2.054064 -0.467410 -0.277325 0.698981 -0.008226];
    K3 = [0.997936 -0.916118 0.199829 0.157810 -0.143368 -0.057385];
    K4 = [-0.127645 0.135543 -0.028919 -0.026546 0.010470 0.012221];

    % Beta factor
    Beta = K0*X + K1*X*rc_125 + K2*X*rc_125^2 + K3*X*rc_125^3 + K4*X*rc_125^4;
    Beta = min([Beta 1]);
    
    % Angle between propeller axis and wing chord [deg]
    % Method for finite wings from Roskam: Methods dor estimation drag
    % polars of subsonic airplanes, S. 3.3.1. Assuming 2D lift curve slope
    % of airfoil sections = 2*pi.
    Mcorr = (1-M^2)^0.5;
    CLalpha = 2*pi*AR/(2+(AR^2*Mcorr^2*(1+(tand(Lambda))^2/Mcorr^2)+4)^0.5);
    alfa = rad2deg(Cl/CLalpha);
    ip = Gamma-alfa;
    
    % Sectional lift coefficient increase (Patterson 2015)
    dCl = 2*pi*((sind(alfa)-Beta*a_c4*sind(ip))*...
        (Beta^2*a_c4^2+2*Beta*a_c4*cosd(alfa+ip)+1)^0.5 - sind(alfa));

    % Plot
    if options == 1
        figure('name','Lift coefficient increase, tractor configuration')
        
        % Streamtube contraction
        subplot(1,2,1)
        xp = -5:0.05:5;
        k_c4_plot = kSS(a_p);
        ac4_plot = a_p./k_c_plot;
        grid on; hold on; box on
        plot(xp,k_c4_plot,'b')
        plot(xp,ac4_plot,'r')
        for aa = [0 0.2 0.4 0.6 0.8 1]
            k_c4_plot = kSS(aa,xp);
            ac4_plot = (1+aa)./k_c4_plot.^2-1;
            plot(xp,k_c4_plot,'--')
            plot(xp,ac4_plot,'-.')
        end
        xlabel('x/R')
        ylabel('Velocity & Area ratio')
        legend('Slipstream radius','V/V_\infty-1',...
            'Slipstream radius, a = 0.0','V/V_\infty-1, a = 0.0',...
            'Slipstream radius, a = 0.2','V/V_\infty-1, a = 0.2',...
            'Slipstream radius, a = 0.4','V/V_\infty-1, a = 0.4',...
            'Slipstream radius, a = 0.6','V/V_\infty-1, a = 0.6',...
            'Slipstream radius, a = 0.8','V/V_\infty-1, a = 0.8',...
            'Slipstream radius, a = 1.0','V/V_\infty-1, a = 1.0')
        plot(xr,k_c4,'ob')
        plot(xr,a_c4,'or')
        axis([min(xp) max(xp) 0 3])
        title(['x/R = ' num2str(xr,'%.3f') ', a = ' num2str(a_p,'%.3f') ...
            ', a_{p} = ' num2str(a_p,'%.3f') ', r_{c/4}/R = ' ...
            num2str(k_c4,'%.3f') ', R/c = ' num2str(Rc,'%.3f') ...
            ', T_c = ' num2str(Tc,'%.3f')])
        
        % Finite slipstream height correction
        subplot(1,2,2)
        rc_array = 0.125:0.125:3;
        a_array = 0.125:0.02:0.625;
        beta_array = zeros(length(rc_array),length(a_array));
        for i = 1:length(rc_array)
            for j = 1:length(a_array)
                X = [1 -xc_125 (-xc_125)^2 -xc_125*(a_array(j)+1) (a_array(j)+1) (a_array(j)+1)^2]';
                beta_array(i,j) = K0*X + K1*X*rc_array(i) + K2*X*rc_array(i)^2 + K3*X*rc_array(i)^3 + K4*X*rc_array(i)^4;
                beta_array(i,j) = min([beta_array(i,j) 1]);
            end
        end
        hold on; box on; grid on
        h1 = surf(rc_array,a_array,beta_array');
        rc_array2 = 0:0.125:3.5;
        a_array2 = 0.005:0.02:0.705;
        beta_array2 = zeros(length(rc_array2),length(a_array2));
        for i = 1:length(rc_array2)
            for j = 1:length(a_array2)
                X = [1 -xc_125 (-xc_125)^2 -xc_125*(a_array2(j)+1) (a_array2(j)+1) (a_array2(j)+1)^2]';
                beta_array2(i,j) = K0*X + K1*X*rc_array2(i) + K2*X*rc_array2(i)^2 + K3*X*rc_array2(i)^3 + K4*X*rc_array2(i)^4;
                beta_array2(i,j) = min([beta_array2(i,j) 1]);                
            end
        end
        hold on; box on; grid on
        surf(rc_array2,a_array2,beta_array2','facecolor','none')
        h2 = scatter3(rc_125,a_125,Beta,'markerfacecolor','r','markeredgecolor','r');
        xlabel('R/c')
        ylabel('a = v_p/v_\infty')
        zlabel('\beta = {\Delta}c_l/{\Delta}c_{l,h=\infty}')
        title(['{\Delta}C_l = ' num2str(dCl,'%.3f') ' (' num2str(100*dCl/Cl,'%.3f') '%)'])
        legend([h1,h2],'Surrogate model limits','Design point')
        
    end

    % 2. Sectional drag coefficient increase, following Biber (2011). 
    % The increase in drag due to propeller thrust has two components: a
    % change in friction drag due to increased dynamic pressure over the
    % wing surface covered by the slipstream, and and increase in
    % lift-induced drag as a consequence of the change in lift due to the
    % propeller slipstream
    Cf = 0.009; % Local skin friction coefficient (see Biber 2011) 
    dCdi_0 = Cf*a_c4^2;
    dCdi_l = (2*Cl*dCl)/pi/AR;
    dCdi = dCdi_0 + dCdi_l;
    
    % 3. Sectional zero-lift drag coefficient increase (generic value)
    dCd0 = 0.000;

    % 4. Propulsive efficiency not affected
    detap = 0;
   
    
%% No effect of propeller on wing 
elseif xc == Inf
    
    % Can be used for comparison (e.g. what would be the difference if
    % there were no effect of the DP system on wing performance), or for
    % e.g. aft-fuselage pylon-mounted propellers, where there is indeed
    % (almost) no effect on the wing. A larger penalty should be given for
    % dCd0, since larger pylons would be required.
    dCl = 0;
    dCdi = 0;
    dCd0 = 0;
    detap = 0;
   
% Wrongly specified propeller location  
else
    error(['Axial propeller position incorrectly specified (x/c = ' ...
            num2str(xc,'%.3f') ')'])
end

%% Check etap is not greater than 1 or smaller than 0
if (etap+detap) > 1 
    error(['Propulsive efficiency is greater than 1 (eta_p = ' ...
            numstr(etap) ', deta_p = ' num2str(detap) ')'])
elseif (etap+detap) <0
    error(['Propulsive efficiency is smaller than 0 (eta_p = ' ...
            numstr(etap) ', deta_p = ' num2str(detap) ')'])
end




