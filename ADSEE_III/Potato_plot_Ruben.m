%%POTATO PLOT
lbs_to_kg = 0.45359237;
mass_pax=175;                    %lbs                       INPUT (fixed)
mass_pax = mass_pax*lbs_to_kg;

mass_bags=25;                    %lbs                       INPUT (fixed)
mass_bags = mass_bags * lbs_to_kg;

mass_fuel=250;                   %lbs                       INPUT
mass_fuel = mass_fuel * lbs_to_kg;

W_OEW=2900;                      %lbs                       INPUT   
W_OEW = W_OEW * lbs_to_kg;
    
cg_OEW=3.7;                     %c.g. Position@OEW          INPUT 

seat_pilot=2 ;                  %c.g. Position Pilot        INPUT 
seat_row1=4  ;                  %c.g. Position Row 1        INPUT 
seat_row2=6  ;                  %c.g. Position Row 2        INPUT 
location_cargo=2.5;             %c.g. Position Baggage      INPUT
location_fuel=5;                %c.g. Position Fuel         INPUT

%cargo
W_OEW_cargo=W_OEW+4*mass_bags;
cg_OEW_cargo=((cg_OEW*W_OEW)+(location_cargo*4*mass_bags))/(W_OEW_cargo);

line([cg_OEW,cg_OEW_cargo],[W_OEW,W_OEW_cargo],'Color','green');

%back to front
W_OEW_1pax=1*mass_pax+W_OEW_cargo;
W_OEW_2pax=2*mass_pax+W_OEW_cargo;
W_OEW_3pax=3*mass_pax+W_OEW_cargo;
W_OEW_4pax=4*mass_pax+W_OEW_cargo;

cg_btf_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row2*mass_pax))/(W_OEW_1pax);
cg_btf_2=((cg_btf_1*W_OEW_1pax)+(seat_row2*mass_pax))/(W_OEW_2pax);
cg_btf_3=((cg_btf_2*W_OEW_2pax)+(seat_row1*mass_pax))/(W_OEW_3pax);
cg_btf_4=((cg_btf_3*W_OEW_3pax)+(seat_row1*mass_pax))/(W_OEW_4pax);

line([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');

%front to back
cg_ftb_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row1*mass_pax))/(W_OEW_1pax);
cg_ftb_2=((cg_ftb_1*W_OEW_1pax)+(seat_row1*mass_pax))/(W_OEW_2pax);
cg_ftb_3=((cg_ftb_2*W_OEW_2pax)+(seat_row2*mass_pax))/(W_OEW_3pax);
cg_ftb_4=((cg_ftb_3*W_OEW_3pax)+(seat_row2*mass_pax))/(W_OEW_4pax);

line([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');

%include fuel

cg_nofuel=cg_btf_4;
cg_fuel=((cg_nofuel*W_OEW_4pax)+(location_fuel*mass_fuel))/(mass_fuel+W_OEW_4pax);

line([cg_nofuel,cg_fuel],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');

title('Potato Plot')
ylabel('mass [kg]')
xlabel('c.g. position from nose [m]')
legend('Cargo','Back to Front','Front to Back','Fuel')

cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])
cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])

%%SCISSOR PLOT