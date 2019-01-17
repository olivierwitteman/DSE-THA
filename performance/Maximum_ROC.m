%% max ROC
clear all
close all
clc 
%ISA
lambda_isa=-0.0065;
P0_isa=101325;
T_0_isa=288.15;
g_isa=9.80665;
R_isa=287;
%other
h2=2400;
h1=15.24;
A=10;
e=0.8;
S=8.8;
W=1700*9.80665;
%important
PA=240000;
%
rho=1.225;
CD0=0.0388;
CD0_takeoff=CD0+0.015;
CD0_land=CD0+0.055;
CL=sqrt(3*pi*A*e*CD0_takeoff);
CD=CD0_takeoff+CL^2/(pi*A*e);
max5=CL^(3/2)/CD;
h=linspace(0,2400)
ROC_mat=[]
for i =h
T=288.15+lambda_isa*i;
P_isa = P0_isa * ((1 + (lambda_isa*i) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa/R_isa/T;
ROC=PA/W-sqrt(W/S)*sqrt(2)/(max5*sqrt(rho_isa));
ROC_mat=[ROC_mat, ROC];
end
disp(max(ROC_mat))
%plot(h,ROC_mat)
%% Minimum rate of Descent
ROD=-sqrt(W/S)*sqrt(2)/(max5*sqrt(rho))


%% Time to climb
V_stall=32.5;
V_takeoff=1.2*V_stall;
V_cruise=92.6;
He1=h1+V_takeoff^2/(2*g_isa);
He2=h2++(V_cruise^2/(2*g_isa));

ROC_avg=(ROC_mat(1)+ROC_mat(100))/2;
t_climb=(He2-He1)/ROC_avg;
t_climb1=(2400-15.24)/ROC_avg;



%% maximum climb angle.

V_theata= 4*W/S*1/(pi*A*e)/(1.225*0.8*(PA/W));
PRmin=W*sqrt(W/S*2/rho_isa*CD^2/CL^3);
Vcl=sqrt(W/S*2/rho_isa*1/CL);
max_angle=(PA-PRmin)/(W*Vcl)*180/pi


%% horizontal distance during climb and descent
normal_angle= (max_angle-5)/180*pi;
save('normal_angle')
s_horizontal=(h2-h1)/tan(normal_angle)
s_horizontal2=(He2-He1)/tan(normal_angle)

%% minimum glide angle
CL_glide= sqrt(CD0_land*A*pi*e);
CD_glide= CD0_land+CL_glide^2/(pi*A*e);
glide_angle=asin(-1/(CL/CD))*180/pi

%% MAximum speed
V=linspace(0,200,1000);
prmat=[];
 
for i=V
    CL=W./S*2./rho_isa*1./(i^2);
    CD=CD0+CL.^2/(pi*A*e);
    PR=W*sqrt(W./S*2/rho_isa*CD.^2./CL.^3);
    prmat=[prmat, PR];
    if PR<1.001*PA && PR>0.999*PA % || PR<0.95*PA
        V_max=i
    end
end

plot(V,prmat)
line([0 200],[PA PA])

