function [] = PR_func(MTOW, OEW , A_payload, L_D, percent_emptiness_payload , M_batt, eta_total, Payload, cb)
%% Input . group 3
g = 9.81;
MTOW=5298*9.80665
M_batt=2053
% MTOW = 1736; %kg <------ INPUT 
% OEW = 1130; %kg  <------ INPUT
% W_Fuel = 232; % kg <------ INPUT
% A_payload = 362 % kg <------ INPUT

% range calculations for propeller aircraft
np = 0.8;   % assumed from class 1 weight estimation excel
   % assumed from class 1 weight estimation excel

% w4w5A = 1./(0.966); % assumed from class 1 weight estimation excel  <------ INPUT
  % <------ INPUT
% L_D = 9; % assumed from slide 42  <------ INPUT

% percent_emptiness_payload = 0.5; % <------ Assumption 


%% Point A Second Try

M_fuel_cruiseA = MTOW ;
range_A = 250; % nmi from requirements
range_A = range_A * 1.852

%% Point B Second Try
L_D=15
range_B= cb*L_D/g*M_batt/(MTOW-percent_emptiness_payload*A_payload)*eta_total/1000 

%% Point C Ferry range

range_C=cb*L_D/g*M_batt/(MTOW - percent_emptiness_payload*A_payload - percent_emptiness_payload*A_payload)*eta_total/1000 

f3 = figure;
plot([0,range_A], [A_payload,A_payload], "Color","blue",'LineWidth',9);
hold on
plot([range_A,range_B], [A_payload, percent_emptiness_payload*A_payload], "Color","red",'LineWidth',9);
plot([range_B,range_C], [percent_emptiness_payload*A_payload, 0],'LineWidth',9);

plot([range_A, range_A], [A_payload, 0],'--', "Color","blue",'LineWidth',9);
plot([range_B, range_B], [percent_emptiness_payload*A_payload, 0],'--', "Color","red",'LineWidth',9);
plot([0, range_B], [percent_emptiness_payload*A_payload, percent_emptiness_payload*A_payload],...
    '--', "Color","red",'LineWidth',7);

ylim([0, 400])
xlabel("Range [km]","FontSize", 40)
ylabel("Payload [kg]","FontSize", 40)
% title("Payload Range Diagram","FontSize", 25)


text(20, A_payload+15, "MTOW with full payload", "FontSize", 40)
text(20, percent_emptiness_payload*A_payload+15, "50% payload", "FontSize", 40)
text(range_C-40, 20, "Ferry range", "FontSize", 40 )


set(gca,'FontSize',70)
end




