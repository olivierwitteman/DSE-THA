%% Speed and climb things
%ISA inputs
R_isa=287.;              %??
g_isa=9.80665;          %gravitational acceleration
lambda_isa=-0.0065;     %temperature gradient
rho0_isa=1.225;         %sea level density
P0_isa=101325;          %sea level pressure
%h=[0., 2400, 1524];                               %INPUT   [m]
h=linspace(0,11000,100)
T_0_isa=288.15;%sea level temperature
Parray=[]
rhoarray=[]

for i = h
    T=288.15+lambda_isa*i;
    P_isa = P0_isa * ((1 + (lambda_isa*i) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
    rho_isa=P_isa/R_isa/T;
    Parray=[Parray, P_isa];
    rhoarray=[rhoarray, rho_isa]
end
%% Stall speed and maximum speed
CLmax=1.5;              %INPUT
W=1557*9.80665;         %INPUT The weight
S=20;                   %INPUT
Pbr=100;               %INPUT Bbrake power of the propellor.
np=1;                 %INPUT Efficiency of the propellor
A=10;
e=0.95;                 %INPUT
T=5000;                 %INPUT
CD0=0.0355;             %INPUT
CLopt_Vmax=sqrt(3/2*CD0*1/(pi*A*e));%INPUT
CD=0.0355+CLopt_Vmax^2/(pi*A*e);
V_stall=sqrt (W./S*2./rhoarray*1./CLmax);
Vmax=sqrt( T./(rhoarray*CD*S)*(1+sqrt(4*CD./(pi*A*e)*(W./T)^2)))
%f1 = figure
plot(V_stall,h)
 
%%Cp= Pbr/rho
%% Maximum climb angle
CD0=0.0355+0.015;
A=1;
e=0.95;
CL=sqrt(CD0*pi*A*e);
CD=CD0+CL^2/(pi*A*e)
clim_angle_max=asin(T_W-1/(CL/CD))*180/pi
k=1/(pi*A*e)
P_W=19;
W_S=1263;
speed=4*W_S*k/(1.225*0.8*P_W)
speed=linspace(-100,100)
iets=speed.^4+np*P_W*W_S./(1.225*CD0)*speed-4*W_S.^2*k./(1.225^2*CD0);
f2 = figure
plot(speed,iets);
disp(4*W_S.^2*k./(1.225^2*CD0))
y=linspace(-100,100);
hold on
%plot(0,y)
line([-500,500], [0, 0], "LineStyle", "--")
%ylim([-20,500])
hold off