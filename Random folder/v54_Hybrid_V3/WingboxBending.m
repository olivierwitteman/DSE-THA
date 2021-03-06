clc 
clear

tic
%function [BOOMAREA] = WingboxBending(c_r, c_t, b, W_en)
c_r = 1.3945;
c_t = c_r*0.4;
b= 9.76;
W_engine = 115;
%Torsion = 10000;
rho = 0.9664;
MCF = -0.06 ;
V_cruise = 92.56;

%function [BOOMAREA] = WingboxBending(c_r, c_t, b)

Torsion = 10000;

 

%A = 0.3; %<---input Boom area
 
x1 = 0.15;  %Positions boom
x2 = 0.295;
x3 = 0.44;
x4 = 0.585 ;
x5 = 0.73 ;
 
y1 = 0.065;
y2 = 0.0804;
y3 = 0.0755;
y4 = 0.0654;
y5 = 0.0473;
y6 = -0.0197;
y7 = -0.0288;
y8 = -0.036;
y9 = -0.0414;
y10 = -0.0412;

dc = (c_r-c_t)/(b/2);
dz = 0.1;
z = [0:dz:b/2];
chl = c_r - dc*z;

mat_x = [x1, x2, x3, x4, x5, x5, x4, x3, x2, x1];
mat_y = [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10];
 
centroid_x = sum(mat_x) /(length(mat_x));
centroid_y = sum(mat_y) /(length(mat_y));
%MOI 
IXX = [];
IYY = [];
IXY = [];

i = 1;
syms A
while i < length(z)+1
 
    I_yy = sum( (((mat_x - centroid_x)*chl(i)).^2))*A;
    I_xx = sum( (((mat_y - centroid_y)*chl(i)).^2))*A;
    I_xy = sum( (((mat_x - centroid_x).*(mat_y - centroid_y)*chl(i))))*A;
    IXX = [IXX, I_xx];
    IYY = [IYY, I_yy];
    IXY = [IXY, I_xy];
   
    i = i+1;
end

%OBTAINING MY DISTRIBUTION
M_y = 4113*(0.5*b-z); %4113 is the max thrust of one engine.

%OBTAINING MX DISTRIBUTION
filename = 'forcesandmoments.xlsx';
sheet = 1;
COLS = 'E2:E485';
MXLIST = xlsread(filename, sheet, COLS);
MXLIST = flip(MXLIST);
MX = 4.4.*(MXLIST(1:10:end))'-4.4*W_engine*9.81*(0.5*b-z);

%OBTAINING TORSION DISTRIBUTION
TMAT = ones(length(z));
TRIMAT = triu(TMAT);
TDIS = 0.5*rho*V_cruise^2*chl*MCF*dz;
TZ = TRIMAT*TDIS';
TORSION = TZ;

%OBTAINING SY DISTRIBUTION
COLSSY = 'B2:B485';
SYLIST = xlsread(filename, sheet, COLSSY);
SYLIST = flip(SYLIST);
SYFRAC = SYLIST(1:10:end)*dz;
SY = TRIMAT*SYFRAC-W_engine*9.81*ones(1,length(z));

%OBTAINING SX DISTRIBUTION
SX = 4113*ones(1,length(z));

%Bending stress 
% g = 9.80665;
% W_wing = 200 * g;%input from layout
% W_engine = 100 * g; %input from layout
% W_hl = 10 * g; %weight of high lift device (all of the weights have to be ajdusted)
% z_L = 3;%input from layout
% z_D = 3;%input from layout
% z_1 = 2;%input from layout
% z_2 = 4;%input from layout
% z_3 = 6;%input from layout
% z_4 = 8;%input from layout
% z_wing = 3.5; %<-- input from CATIA
% T_hl = 200;
% T_E = 5000;

y_max = max(mat_y); % highest distance from  centroid (FROM dat points in excel)_
x_max = max(mat_x); % furtherst distance from centroid

y = (y_max - centroid_y);
x = (x_max - centroid_x);
yb = (mat_y - centroid_y);
xb = (mat_x - centroid_x);
 
%M_x = (L * z_L - W_hl * z_1 - W_hl * z_2 - W_hl * z_3 - W_hl * z_4 - W_wing * z_wing - W_engine * b);
%M_y = (T_hl*(z_1+z_2+z_3+z_4)+T_E*b-D*z_D);
%M_y = 4113*(0.5*b-z); %4113 is the max thrust of one engine.
    
j = 1;
MAXSHEAR = [];
TAU = [];

while j < length(z)+1

    [shear_flow, shear_stress, max_shear, A ]= function_shear_calc(chl(j), SX(j), SY(j), IXX(j), IYY(j), IXY(j), 0.02, A, TORSION(j));

    MAXSHEAR = [MAXSHEAR, max_shear];
    TAU = [TAU, shear_stress];
    j = j+1;
end
% M_z = L*(x_cg-x_ac) + T_E *(y_cg-y_T) + D*(y_cg-y_D) + T_hl*(y_cg-y_Thl);
%sigma_x_max = ((M_z*I_yy-M_y*I_yz)*y+(M_y*I_zz-M_z*I_yz)*z)/(I_yy*I_zz-I_yz^2);
%sigma_y_max = ((M_x*I_zz-M_z*I_xz)*z+(M_z*I_xx-M_x*I_xz)*x)/(I_zz*I_xx-I_xz^2);

k = 1;
BOOMAREA = [];
sigmayield= 78*10^6;
j = 1;
j=linspace(1,10,10);
BAM = [];
while k < length(z)+1
    BOOMAREA = [];
    for m = j
        sigzZ = ((MX(k).*IYY(k)-M_y(k).*IXY(k)).*yb(m)+(M_y(k).*IXX(k)-MX(k).*IXY(k)).*xb(m))/((IXX(k).*IYY(k)-IXY(k).^2));
        eqn1 = sigmayield == sqrt(sigzZ.^2+3.*TAU(m,k).^2);
        BOOMAREA = [BOOMAREA, double(solve(eqn1,A))];  
    end
    BAM = [BAM, max(max(BOOMAREA))]
    k = k+1
end

BAM = BAM(BAM>=0);
ASIZE = max(BAM)*10^6 %MM2

figure
plot(z,MX,'LineWidth', 3)
xlabel('z [m]')
ylabel('Mx [Nm]')
figure
plot(z,M_y,'LineWidth', 3)
xlabel('z [m]')
ylabel('My [Nm]')
figure
plot(z,SX,'LineWidth', 3)
xlabel('z [m]')
ylabel('Sx [Nm]')
figure
plot(z,SY,'LineWidth', 3)
xlabel('z [m]')
ylabel('Sy [Nm]')

BOOMAREA = BOOMAREA(BOOMAREA>=0)
N = 333
toc