clear all
vars = load('../ADSEE_I/variables_ADSEE_I.mat');
thick = 2;
% Weights based on engine selection and preliminary design (might change)
OEW = double(vars.OEW)*9.81;  %operational empty weight [N] <---- INPUT
MTOW = double(vars.MTOW)*9.81; %maximum take-off weight [N] <---- INPUT
Payload  = 363*9.81; %total payload [N]
Fuel = double(vars.W_fuel_total)*9.81;                    % <---- INPUT
%Dimensions and parameters (fixed) 
b = double(vars.b); %span [m]
S = double(vars.S); %reference surface area [m^2]
A = double(vars.S);   %aspect ratio [-]
M = double(vars.M_cruise); %mach number at cruise [-]
beta=sqrt(1-M^2);
bf = 1.7; %fuselage width [m]                   ?????? <----   INPUT
hf = 1.65; %fuselage height [m]                  ?????? <----   INPUT
bh = 4; %horizontal tail span [m]             ?????? <----   INPUT

lf = 7.5; % total length of fuselage [m] .      ?????? <-----  INPUT
Sh = double(vars.S_h); %horizontal tail area [m^2]
Sh = 2.4;


Sv = double(vars.S_v);  %vertical tail area [m^2]
Sv = 1.65;


Cr = double(vars.cr); %main wing root chord [m]
Ct = double(vars.ct); %main wing tip chord [m]
Cr_h = 0.4 * Cr; %horizontal tail root chord [m]   ??????????????????????


Ct_h = 0.2 * Cr_h; %horizontal tail tip chord [m]  ??????????????????????


lfn = 2.2; %distance from nose to leading edge of root chord [m] ?????????????????????????????????????????????????

cg = S/b; %average constant chord [m]

ln = 2.34; %distance from engine to quater chord mac [m] .   ?????????????????????????????????????????????????

bn = 0.0; %width of nacelles (engines) [m]                   
lambda_h = Ct_h/Cr_h; %taper ratio of horizontal tail [-]
lambda = Ct/Cr; %taper ratio of main wing [-]
% Ah = bh^2/Sh; %aspect ratio of horizontal tail [-]
Ah = double(vars.A_h); %aspect ratio of horizontal tail [-] 
Ah = 5.6;

Snet = S-bf*((Cr+1.35)/2); %net area (excluding the eclosed wing area in the fuselage) [m^2] ??????
eta = 0.95; %airfoil efficiency factor [-]
l_press = 0; % length of pressurized area (assumed) [m] 
%-------------------------------------------------------------------------
% Calculate wing sweep angle
sweep_LE = double(vars.sweep_LE)%*pi/180; %sweep at leading edge [rad] %needs to be changed
sweep_14 = atan(tan(sweep_LE) + (Cr/(2*b))*(lambda -1)); %sweep at quater chord [rad]
sweep_12 = atan(tan(sweep_LE) - (4/A)*(0.5*((1-lambda)/(1+lambda)))); %sweep at half chord [rad]
%-------------------------------------------------------------------------
% Calculate horizontal tail sweep angle
sweep_LE_h = double(vars.sweep_LE)%*pi/180; %sweep at leading edge [rad] %needs to be changed
sweep_14_h = atan(tan(sweep_LE_h) + (Cr_h/(2*bh))*(lambda_h - 1)); %sweep at quater chord [rad]
sweep_12_h = atan(tan(sweep_LE_h) - (4/Ah)*(0.5*((1-lambda_h)/(1+lambda_h)))); %sweep at half chord [rad]
%-------------------------------------------------------------------------
% Calculate wing mean aerodynamic chord (mac)
mac = (2/3)*Cr*( (1 + lambda + lambda^2)/(1+lambda)); % mean aerodynamic chord [m]
mac = double(vars.MAC); %%%%% from class I
y_mac = (b/6)*((1 + 2*lambda)/(1 + lambda)); %y location of mac [m]
x_mac = y_mac*tan(sweep_LE);          %x location of mac [m] 
x_datum = 2.7 ;                  % <----- INPUT    measured from planform for given geometry (from nose to wing) [m] ??????? ASK SUMANT
%--------------------------------------------------------------------------
%Calculate horizontal tail mean aerodynamic chord
mac_h = (2/3)*Cr_h*( (1 + lambda_h + lambda_h^2)/(1+lambda_h)); % mean aerodynamic chord of the horizontal tail [m] 
y_mac_h = (bh/6)*((1 + 2*lambda_h)/(1 + lambda_h)); % y location of mac of horizontal tail [m]
x_mac_h = y_mac_h*tan(sweep_LE_h); %x location of mac of horizontal tail [m] 
x_datum_h =7.6 ;                  % <------ INPUT     measured from planform for given geometry (from nose to horizontal tail) [m]
%--------------------------------------------------------------------------
lh = x_datum_h + x_mac_h + 0.25*mac_h - (x_datum + x_mac + 0.25*mac); % ??????? distance between aerodynamic center of main wing and horizontal tail [m] 
lh = lh + 1 % <----- INPUT ???????????????????
lh = 3.8; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x_lemac = x_datum + x_mac;          % x location of leading edge mean aerodynamic chord [m] ???????????
x_lemac = 3.25;

x_OEW = x_lemac + 0.375*mac;        % ?????? assumed CG of operational empty weight [m] >>>>>>>>> HAS TO BE SAME
x_OEW = (3.4 - x_lemac)/mac


x_Cargo = 5.5;                       % ?????  assumed CG of cargo in meters [m]

x_Fuel = x_lemac + 0.5*mac;         % ?????? CG of fuel [m]
% x_Fuel = 2.63;
%--------------------------------------------------------------------------

%% scissor plot
%%% Variables Stability (CRUISE CONFIGURATION *- deg Aoa, - deg Incidence, V = m/s)
SM = 0.05;
%stability margin for safety given as percentage/100
CLa_w = (2*pi*A)/(2 + sqrt(4+(A/eta)^2*(1 +(tan(sweep_12)^2/beta^2)) )) ; 
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for main wing
CLa_h = (2*pi*Ah)/(2+ sqrt(4 + (Ah/eta)^2 *(1+(tan(sweep_12_h)^2/beta^2)) ));
% CLa_h = 2.0; %(for verification)
% dCl/dalpha using DATCOM method [1/rad]
%compressibility ignored due to low speeds

%for horizontal tail
CLa_A_h = CLa_w*(1+(2.15*bf/b))*(Snet/S) + ((pi/2)*(bf^2/S)); 
%CLa_A_h = 8.0; %(for verification)
%dClDalpha for tail-less aicraft 

zh = 0.600 + 0.4881 ;%vertical distance between wing and tail root chord taken from current geometry [m]
m_tv = 2*zh/b; % distance factor between horizontal tail and vortex shed plane of main wing [-]
r = 2*lh/b; % distance factor quarter chord main wing and tail [-]

%ADDITION FOR PROPELLER 
rho= 0.736; % density at given altitude [kg/m^3]
Pbr= 132; %shaft horse power of one engine 132HP = 99000W (Rotax 915)
Cl= 0.645 ;%lift coefficient at given altitude for AoA = 0 with incidence -2 deg!
phi= asin(m_tv/r)*180/pi; %angle between r and m_tv

delta_s_de_da = 6.5*((rho*Pbr^2*S^3+Cl^3)/(lh^4*MTOW^3))^(0.25)*(sin(6*phi))^(2.5); %downwash propeller factor

ked = ((0.1121+0.1265*sweep_14+0.1766*sweep_14^2) / r^2 ) + 0.1024/r + 2;  %downwash corrective coefficient 
ked0 = 0.1124/r^2 + 0.1024/r + 2; %downwash corrective coefficient 
de_da = delta_s_de_da*0 + (ked/ked0)*( (r/(r^2 + m_tv^2))*(0.4876/sqrt(r^2 + 0.6319 + m_tv^2))+...
    (1+(r^2/(r^2 + 0.7915+5.0734*m_tv^2))^0.3113)*(1-sqrt(m_tv^2/(1+m_tv^2))))...
    *(CLa_w/(pi*A)); %total downwash with added delta_s controbution for propeller !!!!!! delta_s_de_da 000


kn = -4.0;% for an engine positioned in front of the lemac/nose propeller
x_ac_w = 0.25; %aerodynamic center of wing (assumed at 0.4mac)
x_ac_c = x_ac_w - ((1.8/CLa_A_h)*(bf*hf*lfn/(S*mac))) +...
    ((0.273/(1+lambda))*((bf*cg*(b-bf))/(mac^2*(b+2.15*bf))))*tan(sweep_14) +...
    2*kn*((bn^2*ln)/(S*mac*CLa_A_h)); %total aircraft aerodynamic center

Vh_V = sqrt(0.85); %flow velocity ratio between H-tail and main wing [-]

%% controllability

CL_A_h = 1.04 + 0.3;% lift coefficient of wing+fuselage (without tail, landing configuration)
% lecture 4 slide 37

dCl_max = 2.0; % change from zero to Clmax
CL_h = -0.35*Ah;% for fixed

CL_0 = 0.8563; % 
Cm0 = -0.216; %  (for main wing) [-]

         %Cmac_w = Cm0*((A*cos(sweep_14)^2)/(A+2*cos(sweep_14))); %pitching moment coefficient at aerodynamic center 
%for wing

         %Cmac_nac = 0.2; %assumed pitching moment coefficient at aerodynamic center (please estimate correctly when engine data available)
%for nacelle

         %Cmac_fus = -1.8*(1-(2.5*bf/lf))*((pi*bf*hf*lf)/(4*S*mac))*(CL_0/CLa_A_h); %pitching moment coefficient at aerodynamic center 
%for fuselage

         %Cm_ac = Cmac_w + Cmac_fus + Cmac_nac; %total moment coefficient at aerodynamic center
         
Cm_ac = -0.082; %(ac assumed to be within +-10% of the neutral point)
%------------------------------------------------------------------------

x_cg_c = (-1:0.01:1);

% Stability Curve
Sh_S = ((x_cg_c) + SM - (x_ac_c))/((CLa_h/CLa_A_h)*(1-de_da)*(lh/mac)*Vh_V^2);

% Neutral Stability Curve (including stability margin)
Sh_S_NS = ((x_cg_c) - (x_ac_c))/((CLa_h/CLa_A_h)*(1-de_da)*(lh/mac)*Vh_V^2 );

% Controllablity Curve
Sh_S_C = ((Cm_ac/CL_A_h) + (x_cg_c)-(x_ac_c)) / ((CL_h/CL_A_h)*(lh/mac)*Vh_V^2);

%-------------------------------------------------------------------------

figure
yyaxis right
% plots
% figure
plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C, "LineWidth", thick)
% title('Scissors-plot: Stability & Controllability Curve')
% xlabel('x_{cg}/MAC [%]')
ylabel('S_h/S [-]', "FontSize", 30)
axis([-1 1 -0.5 0.6])
% legend('Stability','Neutral Stability','Controllability')
hold on


%% POTATO PLOT
MAC = double(vars.MAC);
vars = load('../ADSEE_I/variables_ADSEE_I.mat');
OEW = double(vars.OEW);
x_lemac = [1: 0.01: 5];

cg_mat = zeros(length(x_lemac),2);
counter = 1
cg_OEW = 2.9 - 0.9
for i  = x_lemac
    lf = 7.5;                                      % <----- INPUT m 
    lbs_to_kg = 0.45359237;
    mass_pax=175;                    %lbs               <----- INPUT (fixed)
    mass_pax = mass_pax*lbs_to_kg;

    mass_bags = 25;                                %lbs <----- INPUT (fixed)
    mass_bags = mass_bags * lbs_to_kg;

    mass_fuel=(vars.W_fuel_total);                 %kg  <----- INPUT

    W_OEW = OEW;                                   %kg  <----- INPUT   


%     cg_OEW=2.9 - 2.9 +counter*0.01 + 0.7*0;                     %c.g. Position@OEW  <----- INPUT ????????????????????????????????
%     cg_OEW = 2.5177;    % 30% of mac
%     cg_OEW = 1.6049+0.20;
%     cg_OEW = 1.5;
    % cg_OEW = cg_OEW - x_lemac
%     cg_OEW = 1.0 - 0.01*counter*lf
%     cg_OEW = 1 + 0.01*counter*lf
%     cg_OEW = 1 + i*0.01*counter

%     cg_OEW = cg_OEW*1.001
    cg_OEW = cg_OEW + 0.001*lf
    cg_OEW = 2.9
   

    seat_pilot=2.2 ;                    %c.g. Position Pilot   <----- INPUT 
    seat_row1=2.7  ;                    %c.g. Position Row 1   <----- INPUT 
    seat_row2=4.0  ;                    %c.g. Position Row 2   <----- INPUT 
    location_cargo = 4.6;             %c.g. Position Baggage <----- INPUT
    location_fuel = 2.63;                %c.g. Position Fuel    <----- INPUT
%     location_fuel = cg_OEW + 0.10*l_fus;
%     location_feul = x_Fuel;
    location_fuel  = i + 0.4*MAC

    %cargo
    W_OEW_cargo = W_OEW+4*mass_bags;
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

    cg_max=max([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);
    cg_min=min([cg_OEW,cg_OEW_cargo,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4,cg_fuel]);

    cg_max = (cg_max - i)/MAC;
    cg_min = (cg_min - i)/MAC;
    cg_mat(counter, 1) = cg_min;
    cg_mat(counter, 2) = cg_max;
    
    counter = counter + 1;
end

% figure
yyaxis left
plot([cg_mat(:,1), cg_mat(:,2)], x_lemac/lf, "LineWidth", thick)

hold on
% plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C)
xlabel("xc_{cg}/mac", "FontSize", 30)
ylabel("x_{LEMAC}/L_{FUS}", "FontSize", 30)

legend("FORWARD CG", "AFT CG", "Stability", "Neutral Stability", "Controlability")
set(gca,'FontSize',30);

% yyaxis right
% % plots
% % figure
% plot(x_cg_c,Sh_S,x_cg_c,Sh_S_NS,x_cg_c,Sh_S_C)
% % title('Scissors-plot: Stability & Controllability Curve')
% % xlabel('x_{cg}/MAC  [%]')
% ylabel('S_h/S [-]')
% axis([-1 1 -0.5 0.6])
% % legend('Stability','Neutral Stability','Controllability')
% hold on


figure
x_lemac = 2.9;
line([([cg_OEW,cg_OEW_cargo]- x_lemac)/MAC],[[W_OEW,W_OEW_cargo]],'Color','green');
line([([cg_OEW_cargo,cg_btf_1,cg_btf_2,cg_btf_3,cg_btf_4]-x_lemac)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','blue');
line([([cg_OEW_cargo,cg_ftb_1,cg_ftb_2,cg_ftb_3,cg_ftb_4]-x_lemac)/MAC],[W_OEW_cargo,W_OEW_1pax,W_OEW_2pax,W_OEW_3pax,W_OEW_4pax],'Color','red');
line([([cg_nofuel,cg_fuel]-x_lemac)/MAC],[W_OEW_4pax, W_OEW_4pax+mass_fuel],'Color','black');

xlabel("x_{cg}/MAC")
ylabel("Mass [kg]")
legend("Cargo", "Front to back", "Back to fron", "Fuel")


prompt_for_cg = 'What is the forward cg position: ';
cg_foward = double(input(prompt_for_cg))


prompt_aft_cg = 'What is the aft cg position: ';
cg_aft = double(input(prompt_aft_cg))

prompt_lemac = 'What is the lemac position: ';
lemac_pos = double(input(prompt_lemac))


prompt_surface = 'What is your area ratio: ';
Sh_ratio = double(input(prompt_surface))

Sh_2 = Sh_ratio * double(vars.S)



