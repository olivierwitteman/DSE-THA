%% Inputs
%ISA inputs
R_isa=287.;              %??
g_isa=9.80665;          %gravitational acceleration
lambda_isa=-0.0065;     %temperature gradient
rho0_isa=1.225;         %sea level density
P0_isa=101325;          %sea level pressure
                              %INPUT   [m]
T_0_isa=288.15;         %sea level temperature
T_isa=288.15  ;  %temperature
%% Other inputs
%W=100.;               %INPUT   [N]
Thrust=5000.;           %INPUT   [N]
D=300.;                %INPUT   [N]
alhph_T =0;             %INPUT   [radian] Assumed to be zero.
h_alt=2400;             %INPUT   [m] 
V_climb=10.;             %INPUT   [m/s]
t=0.;
Ve=linspace(0,100);
Pa=10000.;  %constant
fuel_flow=1.
tarray=[];
harray=[];
ROCarray=[];
h=0.;
P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)))
 
delta_h=1000;
for j=Ve
    
    M=j*sqrt(rho0_isa/(1.4*P_isa));
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
    ROCarray=[ROCarray, ROC_unsteady_steady];
    
end
h=1000;
ROC2array=[];
for j=Ve
    P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
    M=j*sqrt(rho0_isa/(1.4*P_isa));
    %ROC_steady=(Pa-D*j)/W; %Rate of climb [m/s]
    ROC_unsteady_steady=1/(1+1.4/2*M^2*(1+R_isa/g_isa*lambda_isa)) ;%unsteady climb   
    harray=[harray,h];
    ROC2array=[ROC2array, ROC_unsteady_steady];
    
end
h=2000;
ROC3array=[];
for j=Ve
    P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
    M=j*sqrt(rho0_isa/(1.4*P_isa));
    %ROC_steady=(Pa-D*j)/W; %Rate of climb [m/s]
    ROC_unsteady_steady=1/(1+1.4/2*M^2*(1+R_isa/g_isa*lambda_isa)) ;%unsteady climb   
    harray=[harray,h];
    ROC3array=[ROC3array, ROC_unsteady_steady];
    
end
h=2400;
ROC4array=[];
for j=Ve
    P_isa = P0_isa * ((1 + (lambda_isa*h) / T_0_isa) ^ (-g_isa / (R_isa*lambda_isa)));
    M=j*sqrt(rho0_isa/(1.4*P_isa));
    %ROC_steady=(Pa-D*j)/W; %Rate of climb [m/s]
    ROC_unsteady_steady=1/(1+1.4/2*M^2*(1+R_isa/g_isa*lambda_isa)) ;%unsteady climb   
    harray=[harray,h];
    ROC4array=[ROC4array, ROC_unsteady_steady];
    
end
%plot(tarray,harray);
plot(Ve, ROCarray, 'LineWidth',2);
hold on
plot(Ve, ROC2array,'LineWidth',2)
plot(Ve, ROC3array,'LineWidth',2)
plot(Ve, ROC4array,'LineWidth',2)
legend('h=0', 'h=1000', 'h=2000', 'h=2400');
ylabel('RC/RC_s') 
xlabel('Equivalent airspeed [m/s]');
ax = gca;
ax.FontSize = 20;
hold off