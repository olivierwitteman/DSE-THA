%% POTATO PLOT
MAC = 1.6;
vars = load('../ADSEE_I/variables_ADSEE_I.mat');
OEW = double(vars.OEW)
x_lemac = [1: 0.01: 5];

cg_mat = zeros(length(x_lemac),2);
counter = 1
for i  = x_lemac
    l_fus = 6.6;                                      % <----- INPUT m 
    lbs_to_kg = 0.45359237;
    mass_pax=175;                    %lbs               <----- INPUT (fixed)
    mass_pax = mass_pax*lbs_to_kg;

    mass_bags = 25;                                %lbs <----- INPUT (fixed)
    mass_bags = mass_bags * lbs_to_kg;

    mass_fuel=(vars.W_fuel_total);                 %lbs <----- INPUT
%     mass_fuel = mass_fuel * lbs_to_kg;

    W_OEW = OEW;                                   %kg  <----- INPUT   


    cg_OEW=3.7;                     %c.g. Position@OEW  <----- INPUT 
    % cg_OEW = 2.5177;    % 30% of mac
%     cg_OEW = 1.6049+0.20;
    cg_OEW = 1.5;
    % cg_OEW = cg_OEW - x_lemac

    seat_pilot=2 ;                  %c.g. Position Pilot   <----- INPUT 
    seat_row1=4  ;                  %c.g. Position Row 1   <----- INPUT 
    seat_row2=6  ;                  %c.g. Position Row 2   <----- INPUT 
    location_cargo=2.5;             %c.g. Position Baggage <----- INPUT
    location_fuel=5;                %c.g. Position Fuel    <----- INPUT
    location_fue = cg_OEW + 0.10*l_fus;

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
    % figure
    % line([([cg_OEW,cg_OEW_cargo]- x_lemac)/MAC],[W_OEW,W_OEW_cargo],'Color','green');
    % line([([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4]-x_lemac)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');
    % line([([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4]-x_lemac)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');
    % line([([cg_nofuel,cg_fuel]-x_lemac)/MAC],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');

    
    
%     figure
%     line([([cg_OEW,cg_OEW_cargo]- x_lemac)/i],[[W_OEW,W_OEW_cargo]],'Color','green');
%     line([([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4]-x_lemac)/i],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');
%     line([([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4]-x_lemac)/i],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');
%     line([([cg_nofuel,cg_fuel]-x_lemac)/i],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');





%     title('Potato Plot')
%     ylabel('mass [kg]')
%     xlabel('c.g. position from nose [m]')
%     legend('Cargo','Back to Front','Front to Back','Fuel')

    cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])
    cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel])

    cg_max = (cg_max - i)/MAC;
    cg_min = (cg_min - i)/MAC;
    cg_mat(counter, 1) = cg_min;
    cg_mat(counter, 2) = cg_max;
    
    counter = counter + 1;
    %%SCISSOR PLOT
end

figure
plot([cg_mat(:,1), cg_mat(:,2)], x_lemac/l_fus)

hold on
% plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C)
xlabel("xc_{cg}/mac")
ylabel("x_{LEMAC}/L_{FUS}")

