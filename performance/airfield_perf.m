%%% AIRFIELD TO PERFORMANCE
clear all
clear clc

%% Assumptions

%% TAKE OFF
%% INPUTS
W_to=1300*9.81          ;                                       %%INPUT 
S=8.8                  ;                                       %%INPUT 
rho=1.225               ;                                       %%INPUT    
%Cl_max=2.               ;                                       %%INPUT 
Cd_to=0.0388+0.015            ;                                       %%INPUT
g=9.80665                 ;                                       %%INPUT    
T_av_to=10000               ;   %  341kW                             %%INPUT                                           %%INPUT 
Dg_av_to=2000             ;                                       %%INPUT 
gamma_to= load('normal_angle.mat')     ;    %%INPUT  
gamma_to=gamma_to.max_angle*pi/180;
H_scr=15.24         ;          	%M]  %input from guide                %%INPUT
mu = 0.02           ;
W=1300*g   ;
Cl=1.7;
A=10;
e_takeoff=0.77;
CD_land=0.0388+0.015+Cl^2/(pi*A*e_takeoff);
Vlimit=sqrt(W/S*2/rho*1/Cl);
V=linspace(0,Vlimit);
Dg=mu*(W-Cl*0.5*rho*V.^2*S);
DG=mean(Dg);
D=0.5*rho*V.^2*S*CD_land;
D_av_to=mean(D);
%% Calculations
%Ground distance
V_min = sqrt(2*W_to/(S*rho*Cl));

V_lof = 1.05*V_min;

D_av_to = 0.5 * rho * S * Cd_to * (V_lof/sqrt(2))^2  ;            

a_av = (g/W_to)*(T_av_to-D_av_to-DG);

X_gd_to = V_lof^2/(2*a_av);

% Transition phase
X_tp_to = V_lof^2*sin(gamma_to)/(0.15*g);

% Climb distance
X_cd_to = (H_scr-(1-cos(gamma_to))*V_lof^2/(0.15*g))/tan(gamma_to);


% TO total distance
X_TO = X_gd_to + X_tp_to + X_cd_to;



%% LANDING
%% INPUTS

W_la = 1000*g               ;                                       %%INPUT
T_av_la = 2000             ;
Cl_land=2.0  ; %%INPUT
V_ap = 1.3*V_min;%%INPUT
V_land=linspace(V_ap, 0,500);
e_landing=0.72;
CD=0.0388+0.055+Cl_land^2/(pi*A*e_landing);
D_land=0.5*rho*S*V_land.^2*CD;
D_av_la=mean(D_land);


gamma_la = 1            ;                                       %%INPUT
Trans_time = 2          ; %[sec] time to transition             %%INPUT
mu = 0.02                  ;                                       %%INPUT
delta_n = 0.15*g     ; %[] 0.15*g for prop                   %%INPUT
T_av_rev = 1            ;                                       %%INPUT
L_av_la = 1             ; 
Cl_land=2.0             ;
V_ap = 1.3*V_min;%%INPUT

%% Calculations
DG_land=mu*(W_la-0.5*rho*S*V_land.^2*Cl_land);
DG_la=DG_land(DG_land>=0);
DG_la=mean(DG_la);
R = 1.3^2*W_la*2/(S*rho*Cl_land*delta_n*g);

X_ab_la = R * sin(gamma_la) + (H_scr-(1-cos(gamma_la))*R)/tan(gamma_la);

% Transition phase
X_tp_la = Trans_time * V_ap;

% Braking phase

X_br_la = W_la^2*2*1.3^2/(2*g*S*rho*Cl_land*(T_av_rev+D_av_la + DG_la));
 
X_LA = X_br_la + X_tp_la + X_ab_la






