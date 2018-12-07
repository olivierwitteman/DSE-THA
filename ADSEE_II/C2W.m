classdef C2W
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end

    methods (Static)
        function W_breakdown = calculation(W_dg,N_z,N_gear,S_w,A,t_over_c,lambda,LAMBDA,S_f,L_D,W_fw,V_cruise,rho,S_ht,LAMBDA_ht,A_ht,lambda_h,H_t_over_H_v,S_vt,LAMBDA_vt,A_vt,lambda_vt,L_t,W_press,W_l,L_m,L_n,W_en,N_en,V_t,V_i,N_t,L,B_w,W_uav,N_p,M)
            q = 0.02088547 * 0.5 * rho * V_cruise^2;
            N_l = 1.5*N_gear;
            W_wing = 0.036*S_w^(0.758) * W_fw^(0.0035) * (A/(cos(LAMBDA))^2.)^(0.6) * q^(0.006) * lambda^(0.04) * (100*t_over_c/cos(LAMBDA))^(-0.3) * (N_z*W_dg)^0.49;
            W_horizontaltail = 0.016*(N_z*W_dg)^(0.414) * q^(0.168) * (10.76 * S_ht)^(0.896) * (100 * t_over_c / cos(LAMBDA))^(-0.12) * (A_ht/(cos(LAMBDA_ht))^2)^(0.043) * lambda_h^(-0.02);
            W_verticaltail = 0.073*(1+0.2*H_t_over_H_v)*(N_z*W_dg)^0.376*q^0.122*(10.76 * S_vt)^0.873*(100*t_over_c/cos(LAMBDA_vt))^-0.49*(A_vt/(cos(LAMBDA_vt))^2)^0.357*lambda_vt^0.039;
            W_fuselage =0.052*S_f^1.086*(N_z*W_dg)^0.177*L_t^-0.051*L_D^-0.072*q^0.241+W_press;
            W_mainlandinggear = 0.095*(N_l*W_l)^0.768*(L_m/12)^0.409;
            W_noselandinggear = 0.125*(N_l*W_l)^0.566*(L_n/12)^0.845;
            W_installedengines = 2.575 * W_en^(0.922) * N_en;
            W_fuelsystem = 2.49 * V_t^(0.726) * (1/(1 + V_i/V_t))^(0.363) * N_t^(0.242) * N_en^(0.157);
            W_flightcontrols = 0.053 * L^(1.536) * B_w^(0.371) * (N_z * W_dg * 10E-4)^(0.8);
            W_hydraulics = 0.001*W_dg;
            W_avionics = 2.117*W_uav^0.933;
            W_electrical = 12.57 * (W_fuelsystem + W_avionics)^0.51;
            W_airco_and_anti_ice = 0.265*W_dg^0.52*N_p^0.68*W_avionics^0.17*M^0.08;
            W_furnishings = 0.0582 * W_dg - 65;
            W_breakdown = [W_wing, W_horizontaltail, W_verticaltail, W_fuselage, W_mainlandinggear, W_noselandinggear, W_installedengines, W_fuelsystem, W_flightcontrols, W_hydraulics, W_avionics, W_electrical, W_airco_and_anti_ice, W_furnishings]/2.2;
        end    
    end
end

