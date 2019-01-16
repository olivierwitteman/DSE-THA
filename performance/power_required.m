%% Maximum rate of climb
close all
clear all
R_isa=287.;              %??
g_isa=9.80665;          %gravitational acceleration
lambda_isa=-0.0065;     %temperature gradient
rho0_isa=1.225;         %sea level density
P0_isa=101325;          %sea level pressure
h=2400     ;               %INPUT   [m]
T_0_isa=288.15;         %sea level temperature
T_isa=288.15 -h*lambda_isa ;  %temperature

P_isa = P0_isa * ((1 + (lambda_isa*h) ./ T_0_isa) .^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa./(R_isa*T_isa);
S=8.8;
V=linspace(0,200);
CD0=0.0338;
CD0_takeoff=0.0338+0.015;
A=10.;
e=0.77;
W=1300*9.80665;
%CL=1./rho_isa*W./S*2./V.^2;
%D=CD0+CL.^2./(pi*A*e);
x=linspace(1,100);
mat=[];
rho_isa=1.225;
mat1=[]
mat3=[]
for i =V
    Drag=0.5*rho_isa.*i^2*S*CD0_takeoff;
    induced= W^2/(0.5*rho_isa*i^2*S)*(1/(pi*A*e));
    %cdi=((W/(0.5*S*i^2*rho_isa)))^2/(pi*A*e)
    cl=(W/(0.5*S*i^2*rho_isa))
    Pr=(Drag+induced)*i;
    %Pr5= 0.5*rho_isa.*i.^3*S*CD0_takeoff+W^2./(0.5*rho_isa.*i.*S)*(1./(pi*A*e))
    mat1=[mat1, cdi];
    mat3=[mat3, cl];
    mat=[mat, Pr];
    
end
plot(cdi, cl)
%plot(V,mat)

