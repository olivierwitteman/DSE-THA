%% Inputs
%ISA inputs
R_isa=287.;              %??
g_isa=9.80665;          %gravitational acceleration
lambda_isa=-0.0065;     %temperature gradient
rho0_isa=1.225;         %sea level density
P0_isa=101325;          %sea level pressure
h=1000.;                               %INPUT   [m]
T_0_isa=288.15;         %sea level temperature
T_isa=288.15  ;  %temperature
%% Other inputs
W=20000.;               %INPUT   [N]
Thrust=5000.;           %INPUT   [N]
D=300.;                %INPUT   [N]
alhph_T =0;             %INPUT   [radian] Assumed to be zero.
h_alt=2400;             %INPUT   [m] 
V_climb=10.;             %INPUT   [m/s]
t=0.;
Ve=linspace(0,200);
Pa=10000.;  %constant
fuel_flow=1.
tarray=[];
harray=[];
ROCarray=[];
P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));

delta_h=100;
for j=Ve
    P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
    M=j*sqrt(rho0_isa/(1.4*P_isa))
    %L=W;
    ROC_steady=(Pa-D*j)/W; %Rate of climb [m/s]
    %flight_angle_steady=asin(ROC_steady/V_climb);   %flight angle [radians]
    ROC_unsteady_steady=1/(1+1.4/2*M^2*(1+R_isa/g_isa*lambda_isa)) ;%unsteady climb 
        %ROC_unsteady_steady=1/(1+0.56*M^2);
    
        
  
    %ROC_unsteady_EAS=ROC_steady*ROC_unsteady_steady;
    %flight_angle_unsteady_EAS=asin(sin(flight_angle_steady)*(1+1.4*M^2/2*(1+R_isa/g_isa*lambda_isa)));
    
    h=h+delta_h    ;          %height change
    %delta_t=delta_h/ROC_unsteady_EAS;
    %t=t+delta_t;%Time Change
    %W=W-delta_t*fuel_flow; %Weight change
    %tarray=[tarray, t];
    harray=[harray,h];
    ROCarray=[ROCarray, ROC_unsteady_steady]
    
end
disp(ROCarray)


%plot(tarray,harray);
plot(Ve, ROCarray);