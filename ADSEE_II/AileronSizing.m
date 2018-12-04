% % classdef ADSEE_II_Aileron
% %     properties
% %     end
% %     methods (Static)
% %         function
% %             
% %         end
% %     end
% % end
% 
% clc; clear all
% 
% 
% %clc; clear all
% %%% ADSEE II - Lecture 4 - Drag coefficient estimation
% P_req = degtorad(60)/1.3;%requirement of roll rate
% 
% %Input here your wing  parameters
% c_r = 1.67; % root chord
% c_t = 0.67; %tip chord
% lambda = 0.; % sweep at leading edge in degrees (positive number)
% theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
% c_l_alpha = 0.32; % Airfoil lift curve slope
% S_ref = 12.3; % Wing surface in square meters
% c_d0 = 0.02; % 2D zero lift drag coefficient
% V = 190.; %speed in m/s
% b = 10.51; %wingspan in meters
% %%%
% %b1 = [0:0.5:(b/2-aileron_length)]; %   the length in meters where the aileron starts measured from the wing root
% %b2 = b1+aileron_length ; % end aileron '
% %%%Aileron geometry input (DO NOT CHANGE)!%%%
% aileron_length = [0:0.05:b/2]; % aileron length in meters
% %b1 = b/2-aileron_length % the length in meters where the aileron starts measured from the wing root
% %b2 = b1+aileron_length % end aileron ''
% tau = 0.6 ; % Function of ratio of the aileron chord over the wing chord (aileron effectiveness) (See slide 10 of ADSEE-II lecture 4 of 2016 for the graph, or look in aircraft design by Mohammed Sadraey)
%             % The aileron should be placed after the rear spar, this
%             % determines the maximum chord ratio
% chordratio_ail_total = [0.075, 0.19, 0.41, 0.7];            
% %chordratio_ail_total = [0.075, 0.19, 0.41, 0.7];
% %tau = [0.2, 0.4, 0.6, 0.8];
% da_max = 30. ; %maximum aileron deflection angle in degrees (reference Mohammed Sadraey)
% %%%
% 
% %function tau, P, min_h_aileron_l,   
% %'For tau =' j 'The turn rate equals:' P,'The minimum (horizontal) aileron length should be (in meters):', aileron_length(counter)),disp('The aileron length itself (inside the wing, induced by sweep) should be at least (in meters):'), disp(aileron_length(counter)/cos(degtorad(theta))) ,disp('The aileron inner chord equals (meters):'), disp((c_r - b1(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)), disp('The aileron outer chord equals (meters):'), disp((c_r - b2(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)) 
%  
% % 
% % for j = [1:1:1]
% %     counter = 1;
% %     for i = [1:1:length(aileron_length)]
% %         syms y;
% %         cy = (c_r - y*(tan(degtorad(lambda))+tan(degtorad(theta))));
% %         ail_sur = int(y*cy,[b1(i) b2(i)]);
% %         C_l_dda = 2.*c_l_alpha*tau(j)/(S_ref*b)*ail_sur; %Derative of the rolling moment coefficient w.r.t. aileron deflection
% %         ail_vol = int(cy*y^2,[0 b/2]);
% %         C_l_p = -4.*(c_l_alpha+c_d0)/(S_ref*b)*ail_vol;
% %         P = -C_l_dda/C_l_p*degtorad(da_max)*(2*V/b);
% %         %disp(double(P));
% %         if P>=P_req
% %             disp('For tau ='),disp(tau(j)),disp('The turn rate equals:'),disp(double(P)), disp('The aileron length itself (inside the wing, induced by sweep) should be at least (in meters):'), disp(aileron_length(counter)/cos(degtorad(theta))) ,disp('The aileron inner chord equals (meters):'), disp((c_r - b1(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)), disp('The aileron outer chord equals (meters):'), disp((c_r - b2(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)) 
% %             break
% %         end
% %         counter= counter + 1; 
% %         
% %     end
% % 
% % end
% %P.S. this code doesn't include differentiable ailerons, nor does it take
% %into account the twist of the wing (which on its place reduces the lift
% %and thus rolling moment) caused by the deflection of the aileron (aileron reversal).
% % It may be better to place the aileron inboard, rather than outboard
% % because we only fly at low speeds. This is done, not to suffer from
% % aileron reversal.
% 
% %After your aileron is designed, check it with reference aircraft, the
% %table in ADSEE-II lecture 2 slide 22 has one.
% 
% b2 = b/2;
% b1 = b2/2;
% P = Intergral(lambda, theta, b1, b2, c_l_alpha, tau, S_ref, b, c_d0, c_r, da_max, V);

classdef AileronSizing
    properties
    end
    methods (Static)
    function [b1, Inner_Ail_Chord, Outer_Ail_Chord] = Iteration(b1_0,b2_0,LAMBDA, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v, P, P_req, chordratio_ail_total) 
        b1 = 0.5*(b1_0 + b2_0)
        while P_req-abs(double(P)) > 0.005
            if double(P)>P_req
                 b1 = 0.5*(b1_0+b1);
            end
            if double(P)<P_req
                b1=0.5*(b1+b2_0);
            end
            b1
            double(P)
            P = AileronSizing.Intergral(LAMBDA, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v);
            if b1 < 0.005
                disp('Aileron cant be implemented')
                break
            end
        end

        Inner_Ail_Chord = (c_r - b1*(tan(degtorad(LAMBDA))+tan(degtorad(theta))))*chordratio_ail_total; %disp('The aileron inner chord equals (meters):'), 
        Outer_Ail_Chord = (c_r - b2*(tan(degtorad(LAMBDA))+tan(degtorad(theta))))*chordratio_ail_total; %disp('The aileron outer chord equals (meters):') 
    end
    
    function P = Intergral(LAMBDA, theta, b1, b2, c_l_alpha, tau, S_ref, b, total_cD0, c_r, da_max, v) 
        syms y;
        cy = (c_r - y*(tan(degtorad(LAMBDA))+tan(degtorad(theta))));
        ail_sur = int(y*cy,[b1 b2]);
        C_l_dda = 2.*c_l_alpha*tau/(S_ref*b)*ail_sur; %Derative of the rolling moment coefficient w.r.t. aileron deflection
        ail_vol = int(cy*y^2,[0 b/2]);
        C_l_p = -4.*(c_l_alpha+total_cD0)/(S_ref*b)*ail_vol;
        P = -C_l_dda/C_l_p*degtorad(da_max)*(2*v/b);
    end
    end
end

%aileron_location_ = b2, b1
%disp('For tau ='),disp(tau(j)),disp('The turn rate equals:'),disp(double(P)), disp('The aileron length itself (inside the wing, induced by sweep) should be at least (in meters):'), disp(aileron_length(counter)/cos(degtorad(theta))) ,disp('The aileron inner chord equals (meters):'), disp((c_r - b1(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)), disp('The aileron outer chord equals (meters):'), disp((c_r - b2(i)*(tan(degtorad(lambda))+tan(degtorad(theta))))*chordratio_ail_total(j)) 
% if double(P)<P_req
%     b1=b1-b1/2
% end
% end




