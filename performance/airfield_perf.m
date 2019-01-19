%%% AIRFIELD TO PERFORMANCE
clear all
clear clc

%% Assumptions

%% TAKE OFF
%% INPUTS
W_to=1434.6*9.81          ;                                       %%INPUT 
S=9.529 ; 
PA=180000;                                       %%INPUT 
rho=1.225               ;                                       %%INPUT    
%Cl_max=2.               ;                                       %%INPUT 
Cd_to=0.0470            ;                                       %%INPUT
g=9.80665                 ;                                       %%INPUT    


gamma_to= load('normal_angle.mat')     ;    %%INPUT  
gamma_to1=gamma_to.normal_angle;
H_scr=15.24         ;          	%M]  %input from guide                %%INPUT
mu = 0.02          ;

Cl=1.7;
A=10;
e_takeoff=0.77;

CD_to=Cd_to+Cl^2/(pi*A*e_takeoff)
Vlimit=sqrt(W_to/S*2/rho*1/Cl);
V=linspace(0,Vlimit);

Dg=mu*(W_to-Cl*0.5*rho*V.^2*S);
DG=(Dg(1)+Dg(100));
%D=0.5*rho*V.^2*S*CD_to;
%D_av_to=mean(D)
%% Calculations
%Ground distance
V_min = sqrt(2*W_to/(S*rho*Cl));

V_lof = 1.05*V_min;
T_av_to=PA/(V_lof*0.8);
D_av_to = 0.5 * rho * S * (Cd_to+Cl^2/(pi*A*e_takeoff)) * (V_lof/sqrt(2))^2              
DG=(W_to-0.5*rho*S*Cl*(V_lof/sqrt(2))^2)*mu
a_av = (g/W_to)*(T_av_to-D_av_to-DG);

X_gd_to = V_lof^2/(2*a_av);

% Transition phase
X_tp_to = V_lof^2*sin(gamma_to1)/(0.15*g);

% Climb distance
X_cd_to = (H_scr-(1-cos(gamma_to1))*V_lof^2/(0.15*g))/tan(gamma_to1)
X_cd_to=(H_scr/tan(gamma_to1))

% TO total distance
X_TO = X_gd_to + X_tp_to + X_cd_to



%% LANDING
%% INPUTS

W_la = 1384.6*g              ;                                       %%INPUT

Cl_land=2.0  ; %%INPUT
V_ap = 1.3*V_min;%%INPUT
T_av_rev=PA/(V_ap*0.4)
V_land=linspace(V_ap, 0,500);
e_landing=0.73;
CD=0.0388+0.055+Cl_land^2/(pi*A*e_landing);
D_land=0.5*rho*S*V_land.^2*CD;
D_av_la=mean(D_land);


gamma_la = 3*pi/180            ;                                       %%INPUT
Trans_time = 3          ; %[sec] time to transition             
mu = 0.05                  ;                                    
delta_n = 0.15 *g    ; %[] 0.15*g for                                       %%INPUT

Cl_land=2.0             ;


%% Calculations
DG_land=mu*(W_la-0.5*rho*S*V_land.^2*Cl_land);
DG_la=DG_land(DG_land>=0);
DG_la=mean(DG_la);
%DG_la=mu*W_la
R = 1.3^2*sqrt(W_la*2/(S*rho*Cl_land))/(delta_n*g);

X_ab_la = R * sin(gamma_la) + (H_scr-(1-cos(gamma_la))*R)/tan(gamma_la);

% Transition phase
X_tp_la = Trans_time * V_ap;

% Braking phase

X_br_la = W_la^2*2*1.3^2/(2*g*S*rho*Cl_land*(T_av_rev+D_av_la + DG_la));
 
X_LA = X_br_la + X_tp_la + X_ab_la






