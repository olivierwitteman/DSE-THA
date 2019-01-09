%% Maximum climb angle and endurance stuff maybe???
np=0.8;
V=50;
cp=0.05;
g=9.80665;
L_D=10;
W4_3=1.5;

Pr=298000
W=1555*g

E=np/(V*g*cp)*L_D*log(W4_3)
%% Maximum climb angle
CD0=0.0355+0.015;
A=10;
e=0.95;
CL=sqrt(CD0*pi*A*e);