function [aileron_length] = AileronNEW(c_r, c_t, sweep_LE, theta, c_l_alpha,...
    S_ref, c_d0, V_stall, b)
disp('AileronNEW')
P_req = degtorad(60)/1.3 ;%requirement of roll rate

% % %Input here your wing  parameters
% c_r = 1.67; % root chord
% c_t = 0.67; %tip chord
% sweep_LE = 0.; % sweep at leading edge in degrees (positive number)
% theta = 10.7773; %sweep at trailing edge in degrees (positive number) (If sweep at leading edge is zero, this equals "atan((c_r-c_t)/(b/2.))"
% c_l_alpha = 0.32; % Airfoil lift curve slope
% S_ref = 12.3; % Wing surface in square meters
% c_d0 = 0.02; % 2D zero lift drag coefficient
% V = 190.; %speed in m/s
% b = 10.51; %wingspan in meters
%%%
%b1 = [0:0.5:(b/2-aileron_length)]; %   the length in meters where the aileron starts measured from the wing root
%b2 = b1+aileron_length ; % end aileron '
%%%Aileron geometry input (DO NOT CHANGE)!%%%
aileron_length = [0:0.05:b/2]; % aileron length in meters
%b1 = b/2-aileron_length; %   the length in meters where the aileron starts measured from the wing root
%b2 = b1+aileron_length; % end aileron ''
tau = [0.5] ; % Function of ratio of the aileron chord over the wing chord (aileron effectiveness) (See slide 10 of ADSEE-II lecture 4 of 2016 for the graph, or look in aircraft design by Mohammed Sadraey)
            % The aileron should be placed after the rear spar, this
            % determines the maximum chord ratio
chordratio_ail_total = [0.286];
da_max = 25; %maximum aileron deflection angle in degrees (reference Mohammed Sadraey)
i=1
V_stall = 190
while i < length(aileron_length)
b1 = b/2-aileron_length; %   the length in meters where the aileron starts measured from the wing root
b2 = b1+aileron_length;  
    
syms y;
cy = (c_r - y*(tan(sweep_LE)+tan(degtorad(theta))));
ail_sur = int(y*cy,[b1(i) b2(i)]);
C_l_dda = 2.*c_l_alpha*tau(1)/(S_ref*b)*ail_sur; %Derative of the rolling moment coefficient w.r.t. aileron deflection
ail_vol = int(cy*y^2,[0 b/2]);
C_l_p = -4.*(c_l_alpha+c_d0)/(S_ref*b)*ail_vol;
P = -C_l_dda/C_l_p*degtorad(da_max)*(2*V_stall/b);

if P>=P_req
    disp('GELUKT')
    disp('For tau ='),disp(tau(1)),disp('The turn rate equals:'),disp(double(P)), disp('The minimum (horizontal) aileron length should be (in meters):'), disp(aileron_length(counter)),disp('The aileron length itself (inside the wing, induced by sweep) should be at least (in meters):'), disp(aileron_length(counter)/cos(degtorad(theta))) ,disp('The aileron inner chord equals (meters):'), disp((c_r - b1(i)*(tan(sweep_LE)+tan(degtorad(theta))))*chordratio_ail_total(1)), disp('The aileron outer chord equals (meters):'), disp((c_r - b2(i)*(tan(sweep_LE)+tan(degtorad(theta))))*chordratio_ail_total(1))
    break
else i = i+1        
end

end

%P.S. this code doesn't include differentiable ailerons, nor does it take
%into account the twist of the wing (which on its place reduces the lift
%and thus rolling moment) caused by the deflection of the aileron (aileron reversal).
% It may be better to place the aileron inboard, rather than outboard
% because we only fly at low speeds. This is done, not to suffer from
% aileron reversal.

%After your aileron is designed, check it with reference aircraft, the
%table in ADSEE-II lecture 2 slide 22 has one.
end


