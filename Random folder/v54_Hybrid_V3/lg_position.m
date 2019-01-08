% function 
perc_mac = 0.4; % 40% of mac assumption
 
l_plane = 7.5;       % <---- INPUT from adsee 1
tailcone_length = 3.2; % <---- INPUT ""
most_aft_cg = 3.55;  % <---- INPUT ""
x_lemac = 3.2;      % <---- INPUT ""
MAC = 1.4;          % <---- INPUT ""
 
l_gear_n = x_lemac + perc_mac*MAC; %assuming sweep is so small it can be neglected
 
lgearrange = l_gear_n - most_aft_cg;
 
length_strut = lgearrange/tan(15*pi/180)
 
if length_strut/(l_plane - l_gear_n) < tan(15*pi/180)
    disp("Good to go")
    disp(length_strut)
    disp(l_gear_n)
else
    disp("Iterate again")
end
    
Dfuse = 1.54; % import fmor ADSEE 1
Tangle = 55*pi/180;
Pm = 0.93*1750*9.81/2; %MTOW from ADSEE assumption
Pn = 0.07*1750*9.81; %MTOW from adsee assumption
 
ln = 2*Pm*lgearrange/(Pn);
ln1 = most_aft_cg - ln
 
% for low wing, the cg is considered to be located at 0.27% of the fuselage
% diameter
 
prompt_winghigh = 'for highwing 2 for lowwing 1: ';
winghigh = double(input(prompt_winghigh))
 
if  winghigh == 1
    z = 0.27*Dfuse+length_strut;
    Ymlg = (ln + lgearrange)/(((ln^2 + (tan(Tangle))^2)/z^2)-1);
    disp(Ymlg)
    
end
 
%for high wing, the cg is considered to be located at 0.626% of the fuselage diameter  
 
if winghigh == 2
    z = 0.626*Dfuse+length_strut;
    Ymlg = (ln + lgearrange)/(((ln^2 + (tan(Tangle))^2)/z^2)-1);
    disp(Ymlg)
end

length_strut

Ymlg