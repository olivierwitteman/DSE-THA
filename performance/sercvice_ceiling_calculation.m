%% Maximum rate of climb

R_isa=287.;              %??
g_isa=9.80665;          %gravitational acceleration
lambda_isa=-0.0065;     %temperature gradient
rho0_isa=1.225;         %sea level density
P0_isa=101325;          %sea level pressure
h=1000;                      %INPUT   [m]
T_0_isa=288.15;         %sea level temperature
T_isa=288.15 -h*lambda_isa ;  %temperature

P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa/(R_isa*T_isa);
S=9;
V=linspace(1,200,200);
CD0=0.0338;
CD0_takeoff=0.0338+0.015;
A=10.;
e=0.8;
W=1526*9.80665;
CL=1./rho_isa*W./S*2./V.^2
D=CD0+CL.^2./(pi*A*e);
PR=0.5*rho_isa*V.^3*S*CD0+W^2./(0.5*rho_isa*V*S)*(1./(pi*A*e))
Ta=10000;
plot(V,PR)
min(PR)