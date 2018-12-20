clc;

perc_mac = 0.4; % 40% of mac

l_plane = 7.5;       % <---- INPUT from adsee 1
tailcone_length = 3; % <---- INPUT ""
most_aft_cg = 3.33;  % <---- INPUT ""
x_lemac = 2.95;      % <---- INPUT ""
MAC = 1.41;          % <---- INPUT ""

l_gear_n = x_lemac + perc_mac*MAC; %assuming sweep is so small it can be neglected

lgearrange = l_gear_n - most_aft_cg;

length_strut = lgearrange/tan(15*pi/180)

if length_strut/(l_plane - l_gear_n) < tan(15*pi/180)
    disp("Good to go")
else
    disp("Iterate again")
end
    
Dfuse = 1.4; % import fmor ADSEE 1
Tangle = 55*pi/180;
Pn = 0.92*Wto/2; %MTOW from ADSEE 
Pm = 0.08*Wto; %MTOW from adsee

ln = 2*Pn*lgearrange/(Pm)

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

    


