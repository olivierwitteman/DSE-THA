function [] = PR_func(MTOW, OEW, W_Fuel, A_payload, L_D, W4W5, percent_emptiness_payload)
%% Input . group 3
g = 9.81;
% MTOW = 1736; %kg <------ INPUT 
% OEW = 1130; %kg  <------ INPUT
% W_Fuel = 232; % kg <------ INPUT
% A_payload = 362 % kg <------ INPUT


% range calculations for propeller aircraft
np = 0.8;   % assumed from class 1 weight estimation excel
cp = 0.6;   % assumed from class 1 weight estimation excel

% w4w5A = 1./(0.966); % assumed from class 1 weight estimation excel  <------ INPUT
W4W5 = 1./(W4W5);   % <------ INPUT
% L_D = 9; % assumed from slide 42  <------ INPUT

% percent_emptiness_payload = 0.5; % <------ Assumption 


%% Point A Second Try
M_endcruiseA = MTOW/W4W5;
M_fuel_cruiseA = MTOW - M_endcruiseA;
range_A = 250; % nmi from requirements
range_A = range_A * 1.852

%% Point B Second Try
M_endcruiseB = MTOW/(MTOW - M_fuel_cruiseA - percent_emptiness_payload*A_payload);
range_B = 375*np./(cp).*L_D*log(M_endcruiseB);
range_B = range_B / 1.150779
range_B = range_B * 1.852

%% Point C Ferry range
M_endcruiseC = (MTOW - percent_emptiness_payload*A_payload)./...
    (MTOW - percent_emptiness_payload*A_payload - percent_emptiness_payload*A_payload - M_fuel_cruiseA)
range_C = 375*np./(cp).*L_D*log(M_endcruiseC);
range_C = range_C / 1.150779
range_C = range_C * 1.852

f3 = figure;
plot([0,range_A], [A_payload,A_payload], "Color","blue",'LineWidth',2);
hold on
plot([range_A,range_B], [A_payload, percent_emptiness_payload*A_payload], "Color","red",'LineWidth',2);
plot([range_B,range_C], [percent_emptiness_payload*A_payload, 0],'LineWidth',2);

plot([range_A, range_A], [A_payload, 0],'--', "Color","blue",'LineWidth',2);
plot([range_B, range_B], [percent_emptiness_payload*A_payload, 0],'--', "Color","red",'LineWidth',2);
plot([0, range_B], [percent_emptiness_payload*A_payload, percent_emptiness_payload*A_payload],...
    '--', "Color","red",'LineWidth',2);

ylim([0,450])
xlabel("Range [km]","FontSize", 10)
ylabel("Payload [kg]","FontSize", 1)
% title("Payload Range Diagram","FontSize", 25)


text(20, A_payload+40, "MTOW with full payload", "FontSize", 40)
text(20, percent_emptiness_payload*A_payload+15, "MTOW with full tanks and 50% payload", "FontSize", 40)
text(range_C-40, 20, "Ferry range", "FontSize", 40 )


set(gca,'FontSize',60)
saveas(f3,'PayloadRange.fig');
close all
end




