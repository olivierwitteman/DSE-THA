function [dCLMAX, TVOL] = planformlayout(b, c_ac_c, c_r, c_t, bf, aileron_length, W_payload, W_fuel, S)

litersperkg = 1/0.721;  %density of avgas inverted
aileronchord_over_chord = 0.27; %follows from aileron effectiesness of 0.5 (see aileron slides adsee II)

fuswidth = bf;
wwidth = b/2 - fuswidth/2;
dc = (c_r-c_t)/(b/2);
z = [0:0.001:b/2];

chl = c_r - dc*z;
chlavgflap = (c_r - dc*0.7 + c_r-dc*(b/2-aileron_length))/2;
fueltanksize = litersperkg*(W_fuel+0.5*W_payload)%
%plot(z,chl)
fls = (wwidth-aileron_length)*chlavgflap;
Swf_S = 2*fls/S;
%AWB = 0.05*chl %Enclosed area wingbox as a function of z assumption

dCLMAX = 0.9*1.3*c_ac_c*Swf_S
TVOL = sum(AWB*0.01) %wingbox volume in m3

end
