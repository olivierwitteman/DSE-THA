function createfigure(X1, yvector1, X2, Y1, X3, X4, YMatrix1, X5, YMatrix2, X6, Y2, X7, Y3)
%CREATEFIGURE(X1, yvector1, X2, Y1, X3, X4, YMatrix1, X5, YMatrix2, X6, Y2, X7, Y3)
%  X1:  area x
%  YVECTOR1:  area yvector
%  X2:  vector of x data
%  Y1:  vector of y data
%  X3:  vector of x data
%  X4:  vector of x data
%  YMATRIX1:  matrix of y data
%  X5:  vector of x data
%  YMATRIX2:  matrix of y data
%  X6:  vector of x data
%  Y2:  vector of y data
%  X7:  vector of x data
%  Y3:  vector of y data

%  Auto-generated by MATLAB on 21-Jan-2019 13:45:37

% Create figure
figure('Name',...
    'Corrected gas turbine comp. wing-loading/power-loading diagram',...
    'Color',[1 1 1]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create area
area(X1,yvector1,'DisplayName','Feasible design space',...
    'FaceColor',[0.9 0.9 0.9],...
    'EdgeColor','none');

% Create plot
% plot(X2,Y1,'DisplayName','Clean wing V_s','Color',[1 0 0]);

% Create plot
plot(X3,Y1,'DisplayName','Powered wing V_s','Color',[0 0 1]);

% Create multiple lines using matrix input to plot
plot1 = plot(X4,YMatrix1);
set(plot1(1),'DisplayName','Cruise speed','Color',[0 1 0]);
set(plot1(2),'DisplayName','Take off (TOP)','Color',[1 0 1]);
set(plot1(3),'DisplayName','Balked landing, OEI 1','Color',[0 1 1]);

% Create multiple lines using matrix input to plot
plot2 = plot(X5,YMatrix2,'Color','none');
set(plot2(1),'DisplayName','Design point for minimum wing size',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'Marker','o');
% set(plot2(2),...
%     'DisplayName','Design point for min. Corrected gas turbine comp. power',...
%     'MarkerEdgeColor',[0 1 0],...
%     'Marker','diamond');
% set(plot2(3),'DisplayName','Design point for min. Fuel power',...
%     'MarkerEdgeColor',[0.8 0 0.3],...
%     'Marker','+');

% Create plot
% plot(X6,Y2,...
%     'DisplayName','Design point for min. Corrected primary EM comp. power',...
%     'MarkerEdgeColor',[0 1 1],...
%     'Marker','v',...
%     'Color','none');

% Create plot
% plot(X7,Y3,...
%     'DisplayName','Design point for min. Corrected secondary EM comp. power',...
%     'MarkerEdgeColor',[0.3 0 0.8],...
%     'Marker','<',...
%     'Color','none');

% Create ylabel
ylabel('Corrected gas turbine comp. power loading \it{W_{\rm{TO}}\it/P_{\rm{GTM}}\rm}\rm [N/W]');

% Create xlabel
xlabel('Wing loading \it{W_{\rm{TO}}\it/S}\rm [N/m^2]');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 10000]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 0.3]);
box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'Layer','top','XTick',...
    [0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000],'YTick',...
    [0 0.05 0.1 0.15 0.2 0.25 0.3]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.727083370484887 0.37262161493274 0.197222222222222 0.293868921775898]);
