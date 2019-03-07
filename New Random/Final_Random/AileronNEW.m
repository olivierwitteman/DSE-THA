function [aileron_length] = AileronNEW(c_r, c_t, sweep_LE, sweep_TE, c_l_alpha,...
    S_ref, c_d0, V_stall, b)

P_req = degtorad(60)/1.3 ;%requirement of roll rate
da_max = 25. ; %maximum aileron deflection angle in degrees (reference Mohammed Sadraey)
tau = 0.5; % Aileron effectiveness
tip_margin = 0; % How much it is away from the tip
b2 = b/2 - tip_margin;
K = 2*c_l_alpha*tau/(S_ref*b);
Q = -4*(c_l_alpha+c_d0)/(S_ref*b);

syms B
eqn_verification = P_req == -K*da_max*(2*1.2*V_stall/b)/Q*(b2^2*c_r/2-b2^3/3*(tan(abs(sweep_LE))+tan(abs(sweep_TE)))-B^2*c_r/2+B^3/3*(tan(abs(sweep_LE))...
    +tan(abs(sweep_TE)))/(b^3/24*c_r-b^4/64*(tan(abs(sweep_LE))+tan(abs(sweep_TE)))));
B=double(solve(eqn_verification, B));

B = B(B>0 & B<b/2);
aileron_length = b2 - B
disp("Aileron starts @"),disp(B),disp('m');
disp('aileron length ='), disp(aileron_length), disp('m');

end


