function [s] = PlotMissionAnalysis(p,m,s,c,aircraft,M,MA)
%% Settings

% Select only mission segments of MA array
names = {'cl','cr','de','Dcl','Dcr','Dde'};
Ncon = size(names,2);

% Colors for weight pie chart
piecolors = [1 1 0
    0 1 0
    0 0.7 0
    0 1 1
    1 0 0
    0 0.2 1
    0 0 0.9
    0 0 0.6];

% % Compute total battery energy available, in case batteries are size by
% % power and therefore some battery energy is left at the end of the mission
% EbatInstalledP = M.bat*p.SE.bat;
% EbatStart = MA.cl.Ebat(1);
% DeltaE = EbatInstalledP - EbatStart;
% if DeltaE < 0
%     DeltaE = 0;
% end


%% I: Mission parameters

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'MA: mission';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% Altitude (if any component limit has been exceeded, show mission segment 
% in dotted lines) 
subplot(2,2,1)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).h,'-b')
    if any(MA.(names{i}).limits==1)==1
        plot(MA.(names{i}).t(MA.(names{i}).limits==1)/60,...
            MA.(names{i}).h(MA.(names{i}).limits==1),'-','color',[1 1 1])
        H = plot(MA.(names{i}).t(MA.(names{i}).limits==1)/60,...
            MA.(names{i}).h(MA.(names{i}).limits==1),'--r');
        legend(H,'Warning: P_{required}>P_{installed}')
    end
end
xlabel('Time [min]')
ylabel('Altitude [m]')

% Mach number
subplot(2,2,2)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).M,'-b')
    if any(MA.(names{i}).limits==1)==1
        plot(MA.(names{i}).t(MA.(names{i}).limits==1)/60,...
            MA.(names{i}).M(MA.(names{i}).limits==1),'-','color',[1 1 1])
        H = plot(MA.(names{i}).t(MA.(names{i}).limits==1)/60,...
            MA.(names{i}).M(MA.(names{i}).limits==1),'--r');
        legend(H,'Warning: P_{required}>P_{installed}')
    end
end
plot([0 max(MA.cl.t)]/60,[m.cr.M m.cr.M],':k')
xlabel('Time [min]')
ylabel('Mach [-]')

% Aircraft mass
subplot(2,2,3)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).W/c.g/1000,'-b')
end
xlabel('Time [min]')
ylabel('Aircraft mass [t]')

% Energies remaining
subplot(2,2,4)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).Ef/1e9,'-r');
    h2 = plot(MA.(names{i}).t/60,(MA.(names{i}).Ebat)/1e9,'-b');
end
xlabel('Time [min]')
ylabel('Energy remaining [GJ]')
legend([h1 h2],'Fuel','Battery','Location','northeast')

% Add division lines of the different segments
for j = 1:4
    for i= 1:Ncon
        subplot(2,2,j)
        ax=gca;
        plot([max(MA.(names{i}).t/60) max(MA.(names{i}).t/60)],...
        ax.YLim,':k')
    end
end


%% II: Performance parameters

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'MA: performance';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% Climb rate
subplot(2,2,1)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).dhdt,'-b')
end
xlabel('Time [min]')
ylabel('Climb rate [m/s]')

% Acceleration
subplot(2,2,2)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).dVdt,'-b')
end
xlabel('Time [min]')
ylabel('Acceleration [m/s^2]')

% Lift coefficients
subplot(2,2,3)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).aero.CL,'-b');
    h2 = plot(MA.(names{i}).t/60,MA.(names{i}).aero.CLiso,'--b');
end
xlabel('Time [min]')
ylabel('Lift coefficient [-]')
legend([h1,h2],'Total','Airframe only','Location','northeast')

% Lift-to-drag ratios
subplot(2,2,4)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).aero.LD,'-b');
    h2 = plot(MA.(names{i}).t/60,MA.(names{i}).aero.CLiso./...
        MA.(names{i}).aero.CDiso,'--b');
end
xlabel('Time [min]')
ylabel('Lift-to-drag ratio [-]')
legend([h1 h2],'Total','Airframe only','Location','northeast')

% Add division lines of the different segments
for j = 1:4
    for i= 1:Ncon
        subplot(2,2,j)
        ax=gca;
        plot([max(MA.(names{i}).t/60) max(MA.(names{i}).t/60)],...
        ax.YLim,':k')
    end
end


%% III: Powertrain parameters

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'MA: powertrain';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% Power-control parameters
subplot(2,2,1)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).P.xi,'-b');
    h2 = plot(MA.(names{i}).t/60,MA.(names{i}).P.phi,'-r');
    h3 = plot(MA.(names{i}).t/60,MA.(names{i}).P.Phi,'-g');
end
xlabel('Time [min]')
ylabel('Power-control parameters [-]')
legend([h1 h2 h3],'Throttle','Supplied power ratio','Shaft power ratio',...
        'Location','northeast')

% Share between climb and accelerate SEP
subplot(2,2,2)
grid on; hold on; box on
for i = 1:Ncon
    plot(MA.(names{i}).t/60,MA.(names{i}).P.x,'-b')
end
xlabel('Time [min]')
ylabel('SEP share, dhdt/(SEP) [-]')

% Fuel and battery power
subplot(2,2,3)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).P.f/1e6,'-r');
    h2 = plot(MA.(names{i}).t/60,MA.(names{i}).P.bat/1e6,'-b');
end

xlabel('Time [min]')
ylabel('Source power [MW]')
legend([h1 h2],'Fuel','Battery','Location','northeast')

% Propulsive power breakdown
subplot(2,2,4)
grid on; hold on; box on
for i = 1:Ncon
    h1 = plot(MA.(names{i}).t/60,MA.(names{i}).P.p/1e6,'-r','linewidth',2);
    h2 = plot(MA.(names{i}).t/60,MA.(names{i}).P.drag/1e6,'-b');
    h3 = plot(MA.(names{i}).t/60,MA.(names{i}).P.climb/1e6,'-m');
    h4 = plot(MA.(names{i}).t/60,MA.(names{i}).P.acceleration/1e6,'-g');
end
xlabel('Time [min]')
ylabel('Flight power [MW]')
legend([h1 h2 h3 h4],'Total propulsive','Drag','Climb','Accelerate',...
    'Location','northeast')

% Add division lines of the different segments
for j = 1:4
    for i= 1:Ncon
        subplot(2,2,j)
        ax=gca;
        plot([max(MA.(names{i}).t/60) max(MA.(names{i}).t/60)],...
        ax.YLim,':k')
    end
end


%% IV: Convergence history
% Note: the DOH plotted in the convergence history is the one used in the
% MA, which assumes zero remaining energy. If the batteries are
% posteriorly sized by power, then this value will not coincide with the
% final DOH result, since there will be energy left in the batteries after
% the mission.

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'MA: convergence';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% Fuel fraction convergence
subplot(2,2,1)
hold on; grid on; box on
plot(0:length(MA.conv.FF)-1,MA.conv.FF,'.k')
xlabel('Iteration number [-]')
if strcmp(p.config,'e-1') || strcmp(p.config,'e-2') ||...
        strcmp(p.config,'dual-e')
    ylabel('Total battery mass fraction [-]')
else
    ylabel('Total fuel mass fraction [-]')
end

% DOH convergence
subplot(2,2,2)
hold on; grid on; box on
plot(0:length(MA.conv.DOH)-1,MA.conv.DOH,'.k')
xlabel('Iteration number [-]')
ylabel('Total degree of hybridization [-]')

% TOM convergence
subplot(2,2,3)
hold on; grid on; box on
plot(0:length(MA.conv.TOM)-1,MA.conv.TOM/1000,'.k')
xlabel('Iteration number [-]')
ylabel('Take-off mass [t]')

% Convergence error
subplot(2,2,4)
hold on; grid on; box on
plot(0:length(MA.conv.err)-1,MA.conv.err,'.k')
xlabel('Iteration number [-]')
ylabel('Error [-]')


%% V: Weight breakdown

% Generate figure
fig = figure(s.figStart + size(s.figs,2));
fig.Name = 'Weight breakdown';
s.figs(size(s.figs,2)+1) = fig;
fig.Color = [1 1 1];

% Prepare divisions and labels of pie chart
X = [M.PL M.OEM M.w M.bat M.f M.GT M.EM1 M.EM2];
Y = 100*X/M.TOM;
Z = {['Payload (' num2str(X(1),'%.0f') ' kg, ' num2str(Y(1),'%.1f') '%)'],...
    ['OEM excl. wing (' num2str(X(2),'%.0f') ' kg, ' num2str(Y(2),'%.1f') '%)'],...
    ['Wing (' num2str(X(3),'%.0f') ' kg, ' num2str(Y(3),'%.1f') '%)'],...
    ['Batteries (' num2str(X(4),'%.0f') ' kg, ' num2str(Y(4),'%.1f') '%)'],...
    ['Fuel (' num2str(X(5),'%.0f') ' kg, ' num2str(Y(5),'%.1f') '%)'],...
    ['Gas turbines (' num2str(X(6),'%.0f') ' kg, ' num2str(Y(6),'%.1f') '%)'],...
    ['Primary EM''s (' num2str(X(7),'%.0f') ' kg, ' num2str(Y(7),'%.1f') '%)'],...
    ['Secondary EM''s (' num2str(X(8),'%.0f') ' kg, ' num2str(Y(8),'%.1f') '%)']};
explode = [0 0 0 0 0 1 1 1];

% Remove divisions which are zero/not used
if sum(X==0)>0
    Z(X==0) = []; explode(X==0) = []; piecolors(X==0,:)=[]; X(X==0) = [];
end

% Generate chart
H = pie(X,explode,Z);
if M.bat == 0
    hT = title(['MTOM = ' num2str(M.TOM,'%.0f') ' kg, S_w = ' num2str(aircraft.Sw,'%.0f') ' m^2']);
else
    if M.bat == M.bat_E
        hT = title({['MTOM = ' num2str(M.TOM,'%.0f') ' kg, S_w = ' num2str(aircraft.Sw,'%.0f') ' m^2'],...
            '(Battery sized by energy requirements)'});
    elseif M.bat == M.bat_P
        hT = title({['MTOM = ' num2str(M.TOM,'%.0f') ' kg, S_w = ' num2str(aircraft.Sw,'%.0f') ' m^2'],...
            '(Battery sized by power requirements)'});
    end;
end
hT.Position = [0 -1.5 0];
if length(H) == 2*length(X)
    for i = 1:length(X)
        hAux = H(2*i-1);
        hAux.FaceColor = piecolors(i,:);
    end
end






