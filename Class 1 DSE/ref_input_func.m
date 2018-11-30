function [OEW_result, Slope, Intercept, average_MTOW] = ref_input_func(chosen_MTOW)
%% reference aircraft input
clc
format shortG

filename = 'Reference_Olivier.xlsx';
sheet = 1;

OEW_weights_set = 'B2:B13'; % where does matlab have to look in the excel
MTOW_weights_set = 'C2:C13'; % where does matlab have to look in the excel
OEW_weights = xlsread(filename,sheet,OEW_weights_set);
MTOW_weights = xlsread(filename,sheet,MTOW_weights_set);

average_MTOW = mean(MTOW_weights)

%% Scatter and regression line
f1 = figure; %%%%%%%%%%%%%
scatter(MTOW_weights, OEW_weights); % 
hold on
hl = lsline;    % plot regression line
set(hl,'color','k')
B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:); % get the slope
Slope = B(2);
Intercept = B(1);


xlabel("MTOW [kg]")
ylabel("OEW [kg]")
title("Reference aircraft MTOW/OEW in kg")

%% Calculating chosen OEW

% chosen_MTOW = 2000    % < -------- Input
OEW_result = chosen_MTOW.*Slope + Intercept


scatter(chosen_MTOW, OEW_result, 200, 'filled')
plot([0, chosen_MTOW], [OEW_result,OEW_result],"--", "Color", "red"); % horiz line
plot([chosen_MTOW, chosen_MTOW], [0, OEW_result],"--", "Color", "red"); % vert line


txt_eq = "Equation: "+string(Slope) + " * MTOW" + " + " + string(Intercept)
text(150, 2000, txt_eq, "FontSize", 20)

txt_point_MTOW = "MTOW: " + string(chosen_MTOW);
text(150, 1800,txt_point_MTOW,"FontSize", 20)

txt_point_OEW = "OE: " + string(round(OEW_result));
text(150, 1600,txt_point_OEW,"FontSize", 20)


end
