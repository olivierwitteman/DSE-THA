clc;

l_plane = 7.5;
tailcone_length = 3;


perc_mac = 0.4; % 40% of mac
most_aft_cg = 3.33;

x_lemac = 2.95;
MAC = 1.41;

l_gear_n = x_lemac + perc_mac*MAC;

something = l_gear_n - most_aft_cg;

length_strut = something/tan(15*pi/180)

if length_strut/(l_plane - l_gear_n) < tan(15*pi/180)
    disp("Good to go")
else
    disp("Iterate again")
end