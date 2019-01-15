clc; clear all
% MOI
%Centroid  
%assuming constant boom area 
h=1.650
w = 1.416
fct=h/w
 
c_r = 1.416
c_t = 0.800
b   = 1.200
 
A = 0.003; %<---input Boom area
 
x1 = -0.50;  %Positions boom (as fct of the w) 
x2 = -0.30;
x3 = -0.10;
x4 = 0.10;
x5 = 0.30;
x6 = 0.50; 
 
y1 = x1 * fct;
y2 = x2 * fct;
y3 = x3 * fct;
y4 = x4 * fct;
y5 = x5 * fct;
y6 = x6 * fct;
 
 
dc = (c_r-c_t)/(b);%taper ratio
z1 = [0:0.01:b];
 
chl = c_r - dc*z1;
 
mat_x = [x4, x5, x6, x6, x6, x6, x6, x6, x5, x4, x3, x2, x1, x1, x1, x1, x1, x1, x2, x3];
mat_y = [y1, y1, y1, y2, y3, y4, y5, y6, y6, y6, y6, y6, y6, y5, y4, y3, y2, y1, y1, y1];
 
 
centroid_x = sum(mat_x * A) /(A*length(mat_x));
centroid_y = sum(mat_y * A) /(A*length(mat_y));
 
tot_Q = [];
%MOI 
IXX = [];
IYY = [];
IXY = [];
Q   = [];
S_y = -8 ; %enter shear force distribution
i = 1;
for i = 1:1: length(z1)
 
    I_yy = sum( (((mat_x - centroid_x)*chl(i)).^2))*A;
    I_xx = sum( (((mat_y - centroid_y)*chl(i)).^2))*A;
    I_xy = sum( (((mat_x - centroid_x).*(mat_y - centroid_y)*chl(i))))*A;
    IXX = [IXX, I_xx];
    IYY = [IYY, I_yy];
    IXY = [IXY, I_xy];
    Q   = [ Q , S_y/I_xx]; % is this per section
    
    q_small = [0]
    counter = 1
    for j = 1:1:length(mat_y)
        
        q   =  -Q(i) * A * (mat_y(j) )*chl(i)
        q_small = [q_small, q + q_small(counter)]
        j = j+1;
        counter = counter + 1
        if counter == 21
            continue
        end
    end
    
    i = i+1;
    tot_Q = [tot_Q, q_small];
end
 
 
 
 
%Bending stress 
 
g = 9.80665;
Load = 10000; %input from xlrf5 (robel)
L_F = 10.51; %input from layout
W_wing = 200 * g;%input from layout
%= 100 * g; %input from layout
W_hl = 10 * g; %weight of high lift device (all of the weights have to be ajdusted)
z_L = 3; %input from layout
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





M_x = L * z_L - W_hl * z_1 - W_hl * z_2 - W_hl * z_3 - W_hl * z_4 - W_wing * z_wing - W_engine * b;
M_y = T_hl*(z_1+z_2+z_3+z_4)+T_E*b-D*z_D;
% 
% M_z = L*(x_cg-x_ac) + T_E *(y_cg-y_T) + D*(y_cg-y_D) + T_hl*(y_cg-y_Thl);
 
sigma_z_max = ((M_x*I_yy-M_y*I_xy)*y+(M_y*I_xx-M_x*I_xy)*x)/(I_xx*I_yy-I_xy^2);
%sigma_x_max = ((M_z*I_yy-M_y*I_yz)*y+(M_y*I_zz-M_z*I_yz)*z)/(I_yy*I_zz-I_yz^2);
%sigma_y_max = ((M_x*I_zz-M_z*I_xz)*z+(M_z*I_xx-M_x*I_xz)*x)/(I_zz*I_xx-I_xz^2);
 
 
 


