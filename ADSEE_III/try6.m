i = [1: 0.01: 4]
vars = load('../ADSEE_I/variables_ADSEE_I.mat');
counter = 1
cg_mat = zeros(length(i),2);
l_fus = 7.5;

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


%     x_lemac = 2.9; % <------ BIGGEST INPUT
MAC = double(vars.MAC);


cg_OEW=((0.4*MAC + x_lemac) - x_lemac)/MAC;                     %c.g. Position@OEW          INPUT 

seat_pilot=(2.2 - x_lemac)/MAC;                  %c.g. Position Pilot        INPUT 
seat_row1=(2.7 - x_lemac)/MAC;                   %c.g. Position Row 1        INPUT 
seat_row2=(4.0 - x_lemac)/MAC;                   %c.g. Position Row 2        INPUT 
location_cargo= (4.6 - x_lemac)/MAC;             %c.g. Position Baggage      INPUT
location_fuel=((x_lemac+MAC*0.4) - x_lemac)/MAC;                %c.g. Position Fuel         INPUT

%cargo
W_OEW_cargo=((W_OEW+4*mass_bags));
cg_OEW_cargo=((cg_OEW*W_OEW)+(location_cargo*4*mass_bags))/(W_OEW_cargo);

% figure
% hold on
%line([cg_OEW,cg_OEW_cargo],[W_OEW,W_OEW_cargo],'Color','green');

%back to front
W_OEW_1pax=1*mass_pax+W_OEW_cargo;
W_OEW_2pax=2*mass_pax+W_OEW_cargo;
W_OEW_3pax=3*mass_pax+W_OEW_cargo;
W_OEW_4pax=4*mass_pax+W_OEW_cargo;

cg_btf_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row2*mass_pax))/(W_OEW_1pax);
cg_btf_2=((cg_btf_1*W_OEW_1pax)+(seat_row2*mass_pax))/(W_OEW_2pax);
cg_btf_3=((cg_btf_2*W_OEW_2pax)+(seat_row1*mass_pax))/(W_OEW_3pax);
cg_btf_4=((cg_btf_3*W_OEW_3pax)+(seat_row1*mass_pax))/(W_OEW_4pax);


%line([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');

%front to back
cg_ftb_1=((cg_OEW_cargo*W_OEW_cargo)+(seat_row1*mass_pax))/(W_OEW_1pax);
cg_ftb_2=((cg_ftb_1*W_OEW_1pax)+(seat_row1*mass_pax))/(W_OEW_2pax);
cg_ftb_3=((cg_ftb_2*W_OEW_2pax)+(seat_row2*mass_pax))/(W_OEW_3pax);
cg_ftb_4=((cg_ftb_3*W_OEW_3pax)+(seat_row2*mass_pax))/(W_OEW_4pax);

%line([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');

%include fuel

cg_nofuel=cg_btf_4;
cg_fuel=((cg_nofuel*W_OEW_4pax)+(location_fuel*mass_fuel))/(mass_fuel+W_OEW_4pax);

%line([cg_nofuel,cg_fuel],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');

title('Potato Plot')
ylabel('mass [kg]')
xlabel('c.g. position from nose [m]')
legend('Cargo','Back to Front','Front to Back','Fuel')

cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);
cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);



cg_mat(counter, 1) = cg_min;
cg_mat(counter, 2) = cg_max;
counter = counter + 1;
%%SCISSOR PLOT


x = linspace(-1,1,302)
figure
i = transpose(i)
plot(i, cg_mat(:,1))
hold on
plot(i, cg_mat(:,2))