function xyz = potato_ruben_123(l_fus, x_lemac, seat_pilot,  seat_row1, seat_row2, location_cargo, location_fuel)
%% POTATO PLOT

vars = load('../ADSEE_I/variables_ADSEE_I.mat');
MAC = double(vars.MAC)


counter = 1
x_lemac = 2.0   % <---- INPUT
l_fus = 6.6;    % <---- INPUT

lbs_to_kg = 0.45359237;
mass_pax=175;                    %lbs               <-----        INPUT (fixed)
mass_pax = mass_pax*lbs_to_kg;

mass_bags = 25;                    %lbs             <-----        INPUT (fixed)
mass_bags = mass_bags * lbs_to_kg;

mass_fuel=vars.W_fuel_total;                   %lbs               <-----        INPUT (fixed)
% mass_fuel = mass_fuel * lbs_to_kg;

W_OEW = double(vars.OEW);                      %lbs               <-----        INPUT (fixed)  
% W_OEW = W_OEW * lbs_to_kg;

% cg_OEW=3.7;                     %c.g. Position@OEW          <-----        INPUT (fixed)
% cg_OEW = 2.5177;    % 30% of mac
cg_OEW = 3.3;

% seat_pilot=2 ;                    %c.g. Position Pilot       <-----    INPUT
% seat_row1=4  ;                    %c.g. Position Row 1       <-----    INPUT
% seat_row2=6  ;                    %c.g. Position Row 2       <-----    INPUT
% location_cargo = 2.5;             %c.g. Position Baggage     <-----    INPUT
% location_fuel = 5;                %c.g. Position Fuel        <-----    INPUT
location_fuel = cg_OEW + 0.10*l_fus;

%cargo
W_OEW_cargo=W_OEW+4*mass_bags;
cg_OEW_cargo=((cg_OEW*W_OEW)+(location_cargo*4*mass_bags))/(W_OEW_cargo);



%back to front
W_OEW_1pax=1*mass_pax+W_OEW_cargo;
W_OEW_2pax=2*mass_pax+W_OEW_cargo;
W_OEW_3pax=3*mass_pax+W_OEW_cargo;
W_OEW_4pax=4*mass_pax+W_OEW_cargo;

cg_btf_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row2*mass_pax))/(W_OEW_1pax);
cg_btf_2=((cg_btf_1*W_OEW_1pax)+(seat_row2*mass_pax))/(W_OEW_2pax);
cg_btf_3=((cg_btf_2*W_OEW_2pax)+(seat_row1*mass_pax))/(W_OEW_3pax);
cg_btf_4=((cg_btf_3*W_OEW_3pax)+(seat_row1*mass_pax))/(W_OEW_4pax);


%front to back
cg_ftb_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row1*mass_pax))/(W_OEW_1pax);
cg_ftb_2=((cg_ftb_1*W_OEW_1pax)+(seat_row1*mass_pax))/(W_OEW_2pax);
cg_ftb_3=((cg_ftb_2*W_OEW_2pax)+(seat_row2*mass_pax))/(W_OEW_3pax);
cg_ftb_4=((cg_ftb_3*W_OEW_3pax)+(seat_row2*mass_pax))/(W_OEW_4pax);


%include fuel
cg_nofuel=cg_btf_4;
cg_fuel=((cg_nofuel*W_OEW_4pax)+(location_fuel*mass_fuel))/(mass_fuel+W_OEW_4pax);


%% graph
figure
line([([cg_OEW,cg_OEW_cargo] - x_lemac)/x_lemac],[[W_OEW,W_OEW_cargo]],'Color','green');
line([([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4]-x_lemac)/x_lemac],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');
line([([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4]-x_lemac)/x_lemac],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');
line([([cg_nofuel,cg_fuel]-x_lemac)/x_lemac],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');


cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])
cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])

xyz = "Ruben's potato finished";
xlabel('x_{cg}/MAC [%]')
ylabel('Mass [kg]')
legend("Cargo", "Back to front", "Front to back", "Fuel")
end