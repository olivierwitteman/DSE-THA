function [dCLMAX, TVOL, flwidth] = planformlayout(b, c_r, c_t, bf, aileron_length, W_payload, W_fuel, S)

litersperkg = 1/0.721;  %density of avgas inverted
aileronchord_over_chord = 0.27; %follows from aileron effectiesness of 0.5 (see aileron slides adsee II)
c_ac_c = 1.1;

fuswidth = bf;
wwidth = b/2 - fuswidth/2;
dc = (c_r-c_t)/(b/2);
z = [0:0.001:b/2];

chl = c_r - dc*z;
chlavgflap = (c_r - dc*0.7 + c_r-dc*(b/2-aileron_length))/2;
fueltanksize = litersperkg*(W_fuel+0.5*W_payload)%
%plot(z,chl)
flwidth = wwidth - aileron_length
fls = (wwidth-aileron_length)*chlavgflap;
Swf_S = 2*fls/S;
AWB = 0.05*chl %Enclosed area wingbox as a function of z assumption

<<<<<<< HEAD

dCLMAX = 0.9*0.9*Swf_S
TVOL = sum(AWB*0.01) %wingbox volume in m3

end



=======
dCLMAX = 0.9*0.9*Swf_S
TVOL = sum(AWB*0.01) %wingbox volume in m3


end
>>>>>>> 11cfa611f46c24842b1bbf1f1912c4939d27fba1
