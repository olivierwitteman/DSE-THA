%% Turning performance.
%INPUTS
close all
clear clc
lambda_isa=-0.0065;
P0_isa=101325;
T_0_isa=288.15;
g_isa=9.80665;
R_isa=287;
h=2400;
%other
h2=2400;
h1=15.24;
A=8;
e=0.8;
S=8.8;
W=1300*9.80665;
%important
PA=290000;
CL_max=1.4;

T=288.15+lambda_isa*h;
P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa/R_isa/T;
V_stall=sqrt(W/S*2/rho_isa*1/CL_max);
CD=(0.0388+CL_max^2/(pi*A*e));
V_PA_PR=(PA/S*2/rho_isa*1/CD)^(1/3)

Vrange1=linspace(V_stall+1,V_PA_PR,500);
L=0.5*rho_isa*S*Vrange1.^2*CL_max;
n=L./W;
bankangle1=acos(1./n);
R=Vrange1.^2./(g_isa*tan(bankangle1));
T_pi=pi.*R./Vrange1;
CDlimit1=0.0388+CL_max^2/(pi*A*e)
CDlimit2=PA/(0.5*rho_isa*120^3*S)
Vrange2=linspace(V_PA_PR, 120,500);

%CD2=PA./(0.5*rho_isa*Vrange2.^3*S);
CD2=linspace(CD, CDlimit2+0.005,500);
CL=sqrt((CD2-0.0388).*pi.*A.*e)
%CL=linspace(1.4,0.25,500)
L2=0.5*rho_isa*S.*Vrange2.^2.*CL;
n2=L2/W;
bankangle2=acos(1./n2);
R2=Vrange2.^2./(g_isa*tan(bankangle2));
T_pi2=pi.*R2./Vrange2;
V=[Vrange1, Vrange2];
bankangletotal=[bankangle1,bankangle2];
R_final=[R, R2];
T_pi_final=[T_pi, T_pi2];
nfinal=[n,n2];
%V=Vrange2
plot(V,R_final)
figure; plot(Vrange2, R2)
hold on
plot(Vrange1, R)
hold off
figure; plot(V,T_pi_final)
figure; plot(V,nfinal)
figure; plot(V,bankangletotal*180/pi)

%plot(V,R_final)
% R_mat=[];
% T_mat=[];
% n_mat=[];
% V=linspace(V_stall+10, 90, 100);
% CL=1.4;
% CD0=0.0488
% CD=CD0+CL^2/(pi*A*e)
%  for i = V
%      Pr=CD*0.5*rho_isa*i^3*S;
%      if Pr>=PA
%          CD=PA/(0.5*rho_isa*i^3*S);
%          CL=sqrt((CD-CD0)*pi*A*e);
%      end
%      L=0.5*CL*rho_isa*i^2*S;
%      n=L/W;
%      bankangle=acos(1/n);
%      R=i^2/(g_isa*tan(bankangle));
%      T_pi=pi*R/i;
%      R_mat=[R_mat,R];
%      T_mat=[T_mat,T_pi];
%      n_mat=[n_mat, n];
%  end
%  
% plot(V,R_mat)
% figure; plot(V,T_mat)
% figure; plot(V,n_mat)