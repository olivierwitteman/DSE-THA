%% max ROC
clear all
close all
clc 
%ISA
lambda_isa=-0.0065;
P0_isa=101325;
T_0_isa=288.15;
g_isa=9.80665;
R_isa=287
%other
A=10;
e=0.8;
S=8.8;
W=1700*9.80665;
PA=250000;
rho=1.225;
CD0=0.0388;
CD0_takeoff=0.0377+0.015;
CL=sqrt(3*pi*A*e*CD0_takeoff);
CD=CD0_takeoff+CL^2/(pi*A*e);
max=CL^(3/2)/CD;
h=linspace(0,2400)
ROC_mat=[]
for i =h
T=288.15+lambda_isa*i;
P_isa = P0_isa * ((1 + (lambda_isa*i) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa/R_isa/T;
ROC=PA/W-sqrt(W/S)*sqrt(2)/(max*sqrt(rho_isa));
ROC_mat=[ROC_mat, ROC];
end
plot(h,ROC_mat)
%% Minimum rate of Descent
ROD=-sqrt(W/S)*sqrt(2)/(max*sqrt(rho))


%% Time to climb
V_stall=32.5;
V_takeoff=1.2*V_stall;
V_cruise=92.6;
He1=15.24+(V_takeoff^2/(2*g_isa));
He2=2400+(V_cruise^2/(2*g_isa));

ROC_avg=(ROC_mat(1)+ROC_mat(100))/2;
t_climb=(He2-He1)/ROC_avg