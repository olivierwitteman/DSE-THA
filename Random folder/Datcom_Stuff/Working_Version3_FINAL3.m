clc; clear all; 
close all

%Sources used:
%Lectures Flight Dynamic// AE3212-I
%Airplane Design 1985, J.Roskam
%Flight Dynamics 2013, J.A. Mulder / W.H.J.J. van Staveren / J.C. van der Vaart / E. de Weerdt /
%C.C. de Visser / A.C. in �t Veld / E. Mooij

%%Inputs
%CRUISE CONDITIONS FOR ALL MANEUVERS 

% Ix = 55229.17/1.3; %Moment of inertia around x-axis [kg*m^2] ????
% Iy = Ix/2.; %Moment of inertia around y-axis [kg*m^2] ???? 4909
% Iz = 58139.42/1.3; %Moment of inertia around z-axis [kg*m^2] ????
% Ixz = 627.60; %Moment of inertia around xz-plane [kg*m^2] ????


Ix = 3111.17; %Moment of inertia around x-axis [kg*m^2] ????
Iy = 4169; %Moment of inertia around y-axis [kg*m^2] ???? 4909
Iz = 6452; %Moment of inertia around z-axis [kg*m^2] ????
Ixz = 400.60; %Moment of inertia around xz-plane [kg*m^2] ????
C_la= 5.58603; %Lift curve slope, taken from data sheet provided by Martin dCl/dalpha[1/rad]  ????
S= 9.52; %Wing surface area [m^2] ???? STAN
c= 1.036; %Mean aerodynamic chord [m] ???? STAN
b= 9.76; %Wing span [m]         % STAN ????
A = 10; %Aspect ratio [-]
g= 9.81; %Gravitational acceleration [m/s^2]
m = 1435. ;%mass [kg}
W = m*g; %Weight [N]
V = 92.6; %Cruise velocity  [m/s] 
C_L = 0.35;% Lift coefficient at 0 deg AoA and -2 deg incidence [-]
rho = 0.966632; %Density at 5000m altitude [kg/m^3]
C_D = 0.0294; %Drag coefficient of the aircraft [-]
T_p = 180000.; %Thrust of the propuslive system [N] , as taken from Vlado
T_c = T_p/(0.5*rho*V^2*S); %Thrust coefficient of the propulsive system [-]
% T_c = -C_D
mu_c=m./(rho*S*c); %Dimesional parameter for symmetric motion 
mu_b=m./(rho*S*b); %Dimesional parameter for asymmetric motion 
Kx2 = (Ix/(rho*S*b^3))/mu_b; % Non-dimesional moment of inertia around x axis [kg*m^2]
Ky2 = (Iy/(rho*S*c^3))/mu_c; % Non-dimenstional moment of inertia around y axis [kg*m^2]
Kz2 = (Iz/(rho*S*b^3))/mu_b; %Non-dimensional moment of inertia around z axis [kg*m^2]
Kxz = (Ixz/(rho*S*b^3))/mu_b;%Non-dimensional moment of inertia around in the xz plane [kg*m^2]
e = 0.82; %Assumed oswald efficiency factor
H = 2400; %Cruise altitude [m]
de_da = 0.378; %Total downwash gradient [-]
Sh = 2.3347; %Surface area of horizontal tail [m^2]
lh = 4.3; %Distance between the aerodynamic center of the wing and the aerodynamic center of the horizontal tail [m]
Vh_V = sqrt(0.85); % Velocity ratio of horizontal tail to wing [-]
C_nha = 5.69; %Lift curve slope of horizontal tail, extrapolated from data by Martin [-]

% % % Derivatives for symmetric flight OUR
alpha_0 = 1.03*pi/180;
i_p = -0.57*pi/180;
T_c = C_D
C_x0 = -C_D + T_c; %Coefficient of X-force in X-direction in equilibrium condition
C_z0 = -C_L - C_D*(alpha_0 + i_p); %Coefficient of Z-force in Z-direction in equilibrium condition
C_xu = -0.001; %-0.001 .  -0.0279; % Derivative of coefficient of X-force in X-direction with respect to velocity
C_zu =  -2.*C_L + C_D*(-(alpha_0 + i_p)); %-0.3762;% Derivative of coefficient of Z-force in Z-direction with respect to velocity
C_mu = 0.0;%0.0699; % Derivative of the pitching moment coefficient with respect to velocity
C_xa =  C_L*(1-(2.*C_la)/(pi*A*e)); %-0.4797; % Derivative of coefficient of X-force in X-direction with respect to angle of attack
C_za = -C_la - C_D; %-5.7434;% Derivative of coefficient of Z-force in Z-direction with respect to angle of attack
C_ma =  -0.37; %-0.66-0.5626; % Derivative of the pitching moment coefficient with respect to angle of attack, from data by Martin
C_zq = -2.*C_nha*Vh_V^2*(Sh*lh)/(S*c)  ; %-5.6629; % Derivative of coefficient of Z-force in Z-direction with respect to pitching velocity
C_mq =   -1.15*C_nha*Vh_V^2*(Sh*lh^2)/(S*c^2); %-8.7941; % Derivative of the pitching moment coefficient with respect to pitching velocity
C_zq = C_mq*1.03/(1.1*lh)
C_zar = -0.2*de_da*Vh_V^2*(Sh*lh)/(S*c) ; %-0.0035; % Derivative of coefficient of Z-force in Z-direction with respect to the change in the angle of attack
C_mar = -0.2*de_da*Vh_V^2*(Sh*lh^2)/(S*c^2) ; %0.1780; % Derivative of the pitching moment coefficient with respect to the change in the angle of attack



% % 
% % % % Derivatives for symmetric flight SLAVENKA
% % C_x0 = -C_D + T_c; %Coefficient of X-force in X-direction in equilibrium condition
% % C_z0 = -C_L; %Coefficient of Z-force in Z-direction in equilibrium condition
% % C_xu = -2.*C_D; %-0.0279; % Derivative of coefficient of X-force in X-direction with respect to velocity
% % C_zu =  -2.*C_L ; %-0.3762;% Derivative of coefficient of Z-force in Z-direction with respect to velocity
% % C_mu = 0.0;%0.0699; % Derivative of the pitching moment coefficient with respect to velocity
% % C_xa =  C_L*(1-(2.*C_la)/(pi*A*e)); %-0.4797; % Derivative of coefficient of X-force in X-direction with respect to angle of attack
% % C_za = -C_la - C_D; %-5.7434;% Derivative of coefficient of Z-force in Z-direction with respect to angle of attack
% % C_ma =  -0.37 ; %-0.5626; % Derivative of the pitching moment coefficient with respect to angle of attack, from data by Martin
% % C_zq = -2.*C_nha*Vh_V^2*(Sh*lh)/(S*c)  ; %-5.6629; % Derivative of coefficient of Z-force in Z-direction with respect to pitching velocity
% % C_mq =   -1.15*C_nha*Vh_V^2*(Sh*lh^2)/(S*c^2); %-8.7941; % Derivative of the pitching moment coefficient with respect to pitching velocity
% % C_zar = -0.2*de_da*Vh_V^2*(Sh*lh)/(S*c) ; %-0.0035; % Derivative of coefficient of Z-force in Z-direction with respect to the change in the angle of attack
% % C_mar = -0.2*de_da*Vh_V^2*(Sh*lh^2)/(S*c^2) ; %0.1780; % Derivative of the pitching moment coefficient with respect to the change in the angle of attack





%Derivatives for asymmetric flight SLAVENKA
%Taken from cessna citation II as a reference aicraft with similar
%geometry,since 3D simulations required for correct derivatives
C_yb = -0.75; % Derivative of coefficient of the Y-force in Y-direction with resprect to the sideslip angle
C_lb = -0.1026;% Derivative of the rolling moment coefficient with respect to the sideslip angle
C_nb =  0.1348;% Derivative of the yawing moment coefficient with respect to the sideslip angle
C_nbr =  0.0; % Derivative of the yawing moment coefficient with respect to the change in the sideslip angle
C_yp = -0.0304;% Derivative of coefficient of the Y-force in Y-direction with resprect to the rolling velocity
C_lp = -0.7108;% Derivative of the rolling moment coefficient with respect to rolling velocity
C_np = -0.0602; % Derivative of the yawing moment coefficient with respect to the rolling velocity
C_yr = 0.8495; % Derivative of coefficient of the Y-force in Y-direction with resprect to the yawing velocity
C_lr = 0.2376; % Derivative of the rolling moment coefficient with respect to yawing velocity
C_nr = -0.2061;% Derivative of the yawing moment coefficient with respect to the yawing velocity



%Derivatives for asymmetric flight OUR VALUES CORRECT for METERS PLUS PLUS
%Taken from cessna citation II as a reference aicraft with similar
%geometry,since 3D simulations required for correct derivatives
C_yb = -0.756; % Derivative of coefficient of the Y-force in Y-direction with resprect to the sideslip angle
C_lb = -0.15378;% Derivative of the rolling moment coefficient with respect to the sideslip angle
C_nb =  +0.163;% Derivative of the yawing moment coefficient with respect to the sideslip angle
C_nbr =  0.0; % Derivative of the yawing moment coefficient with respect to the change in the sideslip angle ??????
C_yp = -0.015;% Derivative of coefficient of the Y-force in Y-direction with resprect to the rolling velocity
C_lp = -0.51;% Derivative of the rolling moment coefficient with respect to rolling velocity
C_np = -0.012; % Derivative of the yawing moment coefficient with respect to the rolling velocity
% C_yr = 0.8495; % Derivative of coefficient of the Y-force in Y-direction with resprect to the yawing velocity
C_lr = 0.1; % Derivative of the rolling moment coefficient with respect to yawing velocity
C_nr = -0.428;% Derivative of the yawing moment coefficient with respect to the yawing velocity


% % Derivatives for asymmetric flight OUR VALUES CORRECT for FEET MINUS
% MINUS
% % Taken from cessna citation II as a reference aicraft with similar
% % geometry,since 3D simulations required for correct derivatives
C_yb = -0.756; % Derivative of coefficient of the Y-force in Y-direction with resprect to the sideslip angle
C_lb = -0.15378;% Derivative of the rolling moment coefficient with respect to the sideslip angle
C_nb =  +0.154;% Derivative of the yawing moment coefficient with respect to the sideslip angle
C_nbr =  0.0; % Derivative of the yawing moment coefficient with respect to the change in the sideslip angle ??????
C_yp = -0.012;% Derivative of coefficient of the Y-force in Y-direction with resprect to the rolling velocity
C_lp = -0.51;% Derivative of the rolling moment coefficient with respect to rolling velocity
C_np = -0.012; % Derivative of the yawing moment coefficient with respect to the rolling velocity
C_yr = 0.8495; % Derivative of coefficient of the Y-force in Y-direction with resprect to the yawing velocity
C_lr = 0.1; % Derivative of the rolling moment coefficient with respect to yawing velocity
C_nr = -0.428;% Derivative of the yawing moment coefficient with respect to the yawing velocity


%%
disp("1: Short Period")
disp("2: Phugoid")
disp("3: Aperiodic Roll")
disp("4: Dutch Roll")
disp("5: Spiral")
prompt=input('Manuevre: ');

%% Short Oscillation
if prompt==1
    %slide 50
    disp('Short Period: ')
    short_period=[2*mu_c*Ky2*(2*mu_c-C_zar), -2*mu_c*Ky2*C_za-(2*mu_c+C_zq)*C_mar-(2*mu_c-C_zar)*C_mq,...
        C_za*C_mq-(2*mu_c+C_zq)*C_ma];
    root_short=roots(short_period);
    root_short=root_short*V/c;
    %Frequency, Period and Oscillations
    w0=sqrt(real(root_short(1)).^(2)+imag(root_short(1)).^(2))*V/c;
    big_zeta=-real(root_short(1))./sqrt(real(root_short(1)).^(2)+imag(root_short(1)).^(2));
    wn=w0*sqrt(1-big_zeta.^(2));
    P=2*pi./wn;
    T12=-0.693/real(root_short(1))*c/V;
    C12=-0.11*(imag(root_short(1))./real(root_short(1)));
    
    X1 = [real(root_short(1));real(root_short(2))];
    YMatrix1 = [imag(root_short(1)),imag(root_short(2))];
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot(X1,YMatrix1,'DisplayName','Symmetric Short Period','MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[1 1 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
    xlabel('Real Axis');
    title('Eigenvalues of a Short Period Motion');
    ylabel('Imaginary Axis');
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
    legend1 = legend(axes1,'show');
    set(legend1,'Location','northwest');
    
    
end
%% Phugoid
if prompt==2
    disp('Phugoid: ')
    phugoid=[-4.*mu_c.^(2), 2.*mu_c.*C_xu, -C_zu.*C_z0];
    root_phugoid=roots(phugoid);
    root_phugoid=root_phugoid*V/c;
    %Frequency, Period and Oscillations
    w0=sqrt(real(root_phugoid(1)).^(2)+imag(root_phugoid(1)).^(2))*V/c;
    big_zeta=-real(root_phugoid(1))./sqrt(real(root_phugoid(1)).^(2)+imag(root_phugoid(1)).^(2));
    wn=w0*sqrt(1-big_zeta.^(2));
    P=2*pi./wn;
    T12=-0.693/real(root_phugoid(1))*c./V;
    C12=-0.11*(imag(root_phugoid(1))./real(root_phugoid(1)));
    
    X2 = [real(root_phugoid(1));real(root_phugoid(2))];
    YMatrix2 = [imag(root_phugoid(1)),imag(root_phugoid(2))];
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot(X2,YMatrix2,'DisplayName','Symmetric Phugoid','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]); 
    xlabel('Real Axis');
    title('Eigenvalues of a Phugoid Motion');
    ylabel('Imaginary Axis');
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
    legend1 = legend(axes1,'show');
    set(legend1,'Location','northwest');
end
%% Aperiodic roll
if prompt==3
    disp('Aperiodic roll: ')
    aperiodic_roll=C_lp./(4*mu_b*Kx2);
    aperiodic_roll = aperiodic_roll*V./b;
    % Time for damp
    T12=-0.693./aperiodic_roll*b./V;
    
    X3 = [real(aperiodic_roll(1))];
    YMatrix3 = [0];
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot(X3,YMatrix3,'DisplayName','Asymmetric Aperiodic Roll','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
    xlabel('Real Axis');
    title('Eigenvalues of a Aperiodic Roll Motion');
    ylabel('Imaginary Axis');
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
    legend1 = legend(axes1,'show');
    set(legend1,'Location','northwest');
end
%% Dutch roll
if prompt==4
    disp('Dutch roll: ')
    dutch_roll=[8*mu_b.^(2)*Kz2, -2*mu_b*(C_nr+2*Kz2*C_yb), 4*mu_b*C_nb+C_yb*C_nr];
    root_dutch=roots(dutch_roll);
    root_dutch = root_dutch.*V./b;
    % Period, Frequency Time damping
    w0=sqrt(real(root_dutch(1)).^(2)+imag(root_dutch(1)).^(2))*V/c;
    big_zeta=-real(root_dutch(1))./sqrt(real(root_dutch(1)).^(2)+imag(root_dutch(1)).^(2));
    wn=w0*sqrt(1-big_zeta.^(2));
    P=2*pi./wn;
    T12=-0.693/real(root_dutch(1))*b/V;
    C12=-0.11*(imag(root_dutch(1))./real(root_dutch(1))); 
    
    X4 = [real(root_dutch(1));real(root_dutch(2))];
    YMatrix4 = [imag(root_dutch(1)),imag(root_dutch(2))];
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot(X4,YMatrix4,'DisplayName','Asymmetric Dutch Roll','MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
    xlabel('Real Axis');
    title('Eigenvalues of a Dutch Roll Motion');
    ylabel('Imaginary Axis');
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
    legend1 = legend(axes1,'show');
    set(legend1,'Location','northwest');
end
%% Spiral 
if prompt==5
    disp('Spiral eigenvalue: ')
    root_spiral=(2*C_L*(C_lb*C_nr-C_nb*C_lr))./(C_lp*(C_yb*C_nr+4*mu_b*C_nb)-C_np*(C_yb*C_lr+4*mu_b*C_lb));
    root_spiral = root_spiral*V./b;
    disp('Spiral Numerator >0 is convergence');  
    (2*C_L*(C_lb*C_nr-C_nb*C_lr));
    disp('Numerator negative,spiral is unstable');
    % Time damping period
    T12=abs(-0.693/real(root_spiral)*b/V);
    
    X5 = [real(root_spiral(1))];
    YMatrix5 = [0];
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    plot(X5,YMatrix5,'DisplayName','Asymmetric Spiral','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
    xlabel('Real Axis');
    title('Eigenvalues of a Spiral Motion');
    ylabel('Imaginary Axis');
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
    legend1 = legend(axes1,'show');
    set(legend1,'Location','northwest');
end

%PLOT ALL MOTIONS ON A SINGLE
%GRAPH-----------------------------------------------------------------------------------------------
short_period=[2*mu_c*Ky2*(2*mu_c-C_zar), -2*mu_c*Ky2*C_za-(2*mu_c+C_zq)*C_mar-(2*mu_c-C_zar)*C_mq,...
    C_za*C_mq-(2*mu_c+C_zq)*C_ma];
root_short=roots(short_period);
root_short=root_short*V/c;

phugoid=[-4.*mu_c.^(2), 2.*mu_c.*C_xu, -C_zu.*C_z0];
root_phugoid=roots(phugoid);
root_phugoid=root_phugoid*V/c;

aperiodic_roll=C_lp./(4*mu_b*Kx2);
aperiodic_roll = aperiodic_roll*V./b;

dutch_roll=[8*mu_b.^(2)*Kz2, -2*mu_b*(C_nr+2*Kz2*C_yb), 4*mu_b*C_nb+C_yb*C_nr];
root_dutch=roots(dutch_roll);
root_dutch = root_dutch.*V./b;

root_spiral=(2*C_L*(C_lb*C_nr-C_nb*C_lr))./(C_lp*(C_yb*C_nr+4*mu_b*C_nb)-C_np*(C_yb*C_lr+4*mu_b*C_lb));
root_spiral = root_spiral*V./b;

X1 = [real(root_short(1));real(root_short(2))];
YMatrix1 = [imag(root_short(1)),imag(root_short(2))];
X2 = [real(root_phugoid(1));real(root_phugoid(2))];
YMatrix2 = [imag(root_phugoid(1)),imag(root_phugoid(2))];
X3 = [real(aperiodic_roll(1))];
YMatrix3 = [0];
X4 = [real(root_dutch(1));real(root_dutch(2))];
YMatrix4 = [imag(root_dutch(1)),imag(root_dutch(2))];
X5 = [real(root_spiral(1))];
YMatrix5 = [0];

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(X1,YMatrix1,'DisplayName','Symmetric Short Period','MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[1 1 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]); hold on;
plot(X2,YMatrix2,'DisplayName','Symmetric Phugoid','MarkerFaceColor',[1 0 0],'MarkerSize',8,'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]); hold on;
plot(X3,YMatrix3,'DisplayName','Asymmetric Aperiodic Roll','MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]); hold on;
plot(X4,YMatrix4,'DisplayName','Asymmetric Dutch Roll','MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]); hold on;
plot(X5,YMatrix5,'DisplayName','Asymmetric Spiral','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);

xlabel('Real Axis');
title('Eigenvalues');
ylabel('Imaginary Axis');
box(axes1,'on');
grid(axes1,'on');
set(axes1,'XAxisLocation','origin','YAxisLocation','origin')
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest');
%----------------------------------------------------------------------------------------------------------------





