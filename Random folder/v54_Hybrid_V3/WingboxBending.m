clc 
clear

tic
%function [BOOMAREA] = WingboxBending(c_r, c_t, b)
c_r = 1.6;
c_t = 1.6*0.4;
b= 11;
<<<<<<< HEAD
Torsion=10000
=======
Torsion = 10000;
 
>>>>>>> 10ad22b343e8790e8a4176a19861e02fda8a15fc
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
z = [0:0.1:b/2];
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

%Bending stress 
g = 9.80665;
L = 10000; %input from xlrf5 (robel)
b = 10.51; %input from layout
W_wing = 200 * g;%input from layout
W_engine = 100 * g; %input from layout
W_hl = 10 * g; %weight of high lift device (all of the weights have to be ajdusted)
z_L = 3;%input from layout
z_D = 3;%input from layout
z_1 = 2;%input from layout
z_2 = 4;%input from layout
z_3 = 6;%input from layout
z_4 = 8;%input from layout
z_wing = 3.5; %<-- input from CATIA
T_hl = 200;
T_E = 5000;
D = 1000; %xflr5 (robel)
y_max = max(mat_y); % highest distance from  centroid (FROM dat points in excel)_
x_max = max(mat_x); % furtherst distance from centroid
x_cg = 0; %Centre of gravity of wing on airfoil (CATIA Stan?)
y_cg = 0; %idem
x_ac = 0; %aerodynamic centre from xflr5 (robel)
y_ac = 0; %idem
y_T = 0; %vertical placing of engines (CATIA STAN)
y_Thl = 0; %vertical placing of high lift engines (CATIA Stan) (probably the same as y_T)
y_D = 0; %vertical position of drag (xflr5 robel)
y = (y_max - centroid_y);
x = (x_max - centroid_x);
A_y = 4*W_hl+W_wing+W_engine-L;
 
M_x = (L * z_L - W_hl * z_1 - W_hl * z_2 - W_hl * z_3 - W_hl * z_4 - W_wing * z_wing - W_engine * b);
M_y = (T_hl*(z_1+z_2+z_3+z_4)+T_E*b-D*z_D);

%sigma_z_max = 511240;
%    eqn1 = sigma_z_max == ((M_x*I_yy-M_y*I_xy)*y+(M_y*I_xx-M_x*I_xy)*x)/(I_xx*I_yy-I_xy^2);
%    BOOMAREA = double(solve(eqn1,A))
    
j = 1;
MAXSHEAR = [];
while j < length(z)+1
<<<<<<< HEAD
    [shear_flow, shear_stress, max_shear, A ]= function_shear_calc(chl(j), 1000, 1000, IXX(j), IYY(j), IXY(j), 0.25, 0.3, 0.25, 0.02,  0.02, A, Torsion);
=======
    [shear_flow, shear_stress, max_shear, A ]= function_shear_calc(chl(j), 1000, 1000, IXX(j), IYY(j), IXY(j), 0.25, 0.3, 0.25, 0.02, 0.02, A, Torsion);
>>>>>>> 10ad22b343e8790e8a4176a19861e02fda8a15fc
    MAXSHEAR = [MAXSHEAR, max_shear];
    j = j+1;
end
% M_z = L*(x_cg-x_ac) + T_E *(y_cg-y_T) + D*(y_cg-y_D) + T_hl*(y_cg-y_Thl);
%sigma_x_max = ((M_z*I_yy-M_y*I_yz)*y+(M_y*I_zz-M_z*I_yz)*z)/(I_yy*I_zz-I_yz^2);
%sigma_y_max = ((M_x*I_zz-M_z*I_xz)*z+(M_z*I_xx-M_x*I_xz)*x)/(I_zz*I_xx-I_xz^2);

k = 1;
BOOMAREA = [];
sigmayield= 511240;

while k < length(z)+1
    sigma=(((M_x.*IYY(k)-M_y.*IXY(k))*y+(M_y.*IXX(k)-M_x.*IXY(k)).*x)/(IXX(k).*IYY(k)-IXY(k).^2)).^2;
    eqn1 = sigmayield == sqrt((((M_x.*IYY(k)-M_y.*IXY(k))*y+(M_y.*IXX(k)-M_x.*IXY(k)).*x)/(IXX(k).*IYY(k)-IXY(k).^2)).^2+3.*MAXSHEAR(k).^2);
    BOOMAREA = [BOOMAREA, double(solve(eqn1,A))];
    k = k+1
end
BOOMAREA = BOOMAREA%(BOOMAREA>=0)
N = 333
toc