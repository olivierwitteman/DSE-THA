% aircraft drawing example

l_fus = 8;
width_fus = 1.6;
height_fus = 1.85;

pilot_cabin_length = 1.4;

tailcone_length = 1.7;

seat_length = 0.4;
seat_heigth = 1.25;

seat_pilot = 1.4;
seat_row1 = 2.5;
seat_row2 = 3.5;

lg_extended = 0.4; 

wing_pos = "top"
wing_pos = lg_extended + height_fus + 0.1
x_lemac = 2.5; % <--- IMPORTANT



line([0, l_fus-tailcone_length], [lg_extended, lg_extended], "LineWidth", 7); % bottom of fuselage
line([pilot_cabin_length, l_fus], ...
    [lg_extended + height_fus, lg_extended + height_fus], "LineWidth", 7); % top of fuselage
line([0, pilot_cabin_length], [lg_extended, lg_extended + height_fus ], "LineWidth", 7); % front diag of fuselage
line([l_fus-tailcone_length, l_fus], [lg_extended,lg_extended + height_fus ], "LineWidth", 7); % back diag of fuselage

xlim([-1, l_fus + 0.5])
ylim([-0.5, lg_extended + height_fus + 0.5])



