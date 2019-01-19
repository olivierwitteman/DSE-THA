%% max ROC
clear all
close all
clc 
%ISA
lambda_isa=-0.0065;
P0_isa=101325;
T_0_isa=288.15;
rho_0_isa=1.225;
h2=2400;
g_isa=9.80665;
R_isa=287;
%other
h1=15.24;
A=10;
e_clean=0.82;
e_landing=0.73;
e_takeoff=0.77;
W_to=1434.6*9.80665;
S=10.23;
W=1395*9.80665;
%important
PA_total=180000;
PA_secondary= 30000;
PA_primary=150000
%
CL_takeoff=1.7;
CL_land=2.;
rho=1.225;
CD0=0.032;
CD0_takeoff=0.0470;
CD0_land=0.0870;
CL=sqrt(3*pi*A*e_takeoff*CD0_takeoff);
CD=CD0_takeoff+CL^2/(pi*A*e_takeoff);
max5=CL^(3/2)/CD;
h=linspace(0,3880,5000);
ROC_mat=[];
rc=(110-190)*0.75*1000/((15000-2500)*0.3048);
for i =h
T=288.15+lambda_isa*i;
P_isa = P0_isa * ((1 + (lambda_isa*i) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
rho_isa=P_isa/R_isa/T;
PA1=rho_isa/rho_0_isa*PA_total+rc*i;
ROC=PA1/W_to-sqrt(W_to/S)*sqrt(2)/(max5*sqrt(rho_isa));
ROC_mat=[ROC_mat, ROC];
end
%disp(max(ROC_mat))
plot(h,ROC_mat)
%% Minimum rate of Descent (T=0)
W_land=1336*9.80665;
ROD=-sqrt(W_land/S)*sqrt(2)/(max5*sqrt(rho))


%% Time to climb
V_stall=sqrt(W_to/S*2/rho_isa*1/CL_takeoff);
V_takeoff=1.2*V_stall;
V_cruise=92.6;
He1=h1+V_takeoff^2/(2*g_isa);
He2=h2++(V_cruise^2/(2*g_isa));

ROC_avg=(ROC_mat(1)+ROC_mat(100))/2;
t_climb=(He2-He1)/ROC_avg;
t_climb1=(2400-15.24)/ROC_avg;



%% maximum climb angle.

V_theata= 4*W_to/S*1/(pi*A*e_takeoff)/(1.225*(PA_total/W_to));
PRmin=W_to*sqrt(W_to/S*2/rho*CD^2/CL^3);
Vcl=sqrt(W_to/S*2/rho*1/CL_takeoff);
max_angle=asin((PA_total-PRmin)/(W_to*Vcl))*180/pi;


%% horizontal distance during climb and descent
normal_angle= (max_angle-6)/180*pi;
save('normal_angle')
s_horizontal=(h2-h1)/tan(normal_angle);
s_horizontal2=(He2-He1)/tan(normal_angle);

%% minimum glide angle
CL_glide= sqrt(CD0_land*A*pi*e_landing);
CD_glide= CD0_land+CL_glide^2/(pi*A*e_landing);
glide_angle=asin(-1/(CL_glide/CD_glide))*180/pi

%% MAximum speed
V=linspace(0,200,10000);
prmat=[];
 
for i=V
    CL=W./S*2./rho_isa*1./(i^2);
    CD=CD0+CL.^2/(pi*A*e_clean);
    PR=W*sqrt(W./S*2/rho_isa*CD.^2./CL.^3);
    prmat=[prmat, PR];
    if PR<1.001*PA_total && PR>0.999*PA_total % || PR<0.95*PA
        V_max=i;
    end
end

%plot(V,prmat)
%line([0 200],[PA PA])

