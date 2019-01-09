%Rate of climb steady and unsteady conditions.
%inputs
R=287;
g=8.80665;
lambda=-0.0065;
rho=1.225;
S=22;
h=0;          %INPUT   [m]
T=288.15+lambda*h;
W=10000;      %INPUT   [N]
T=5000;       %INPUT   [N]
D=3000;       %INPUT   [N]
alhph_T =0;   %INPUT   [radian] Assumed to be zero.
h_alt=2400;   %INPUT   [m] 
V_climb=10;   %INPUT   [m/s]
a=sqrt(1.4*R*T);
M=V_climb/a;
%% Steady Climb
L=W;
ROC_steady=(L-D)*V_climb/W ; %Rate of climb [m/s]
flight_angle_steady=asin(ROC_steady/V_climb) ;  %flight angle [radians]

%% Unsteady Climb
% assumption: quasi rectilinear change in flight angle is small, but change
% in velocity is not negated. This is only for a certain moment in time.
%for constant EAS
L=W;
ROC_unsteady_EAS=ROC_steady*(1+1.4*M^2/2*(1+R/g*lambda));
flight_angle_unsteady_EAS=asin(sin(flight_angle_steady)*(1+1.4*M^2/2*(1+R/g*lambda)));
% for constant mach number
ROC_unsteady_mach=ROC_steady*(1+0.5*1.4*M^2*R/g*lambda);
flight_angle_unsteady_mach=asin( sin(flight_angle_steady)*(1+0.5*1.4*M^2*R/g*lambda));

%% Gliding
CLopt=1.5;
CD=0.2;
CL=0.5;
glide_angle=asin(CD/CL);
V_glide=sqrt(W/S*2/rho*1/CLopt);

%% Energy height gedoe for minimum time to climb.
%0 to 1400 m with 
fuel_consuption=2;
Ve1=10;
Ve2=20;
H1=0;
H2=2400;
He1=H1+Ve1^2/(2*g);
He2=H2+Ve2^2/(2*g);
t_climb=(He2-He1)/ROC_steady;
S_horizontal=(He2-He1)/tan(flight_angle_steady);

%% Path performance of unsteady climb. 
h=0

