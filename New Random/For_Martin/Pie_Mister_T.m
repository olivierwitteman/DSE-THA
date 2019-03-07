%% Composite part


engineering = 84518.48835;
flight_test = 342.8701172;
manufacturing = 355751.8182;
tooling = 49883.19465;
quality_control = 70834.10455;
materials = 19067.16062;
% lg = -7500;
avionics = 70000;
powertrain = 150274;

total = 800672;


labels = {"Engineering " + num2str(round(engineering/total*100, 1)) + "%"...
    ,"Flight Test "+ num2str(round(flight_test/total*100, 1)) + "%"...
    ,"Manufacturing " + num2str(round(manufacturing/total*100, 1))+ "%"...
    ,"Tooling " + num2str(round(tooling/total*100, 1))+ "%"...
    ,"Quality Control " + num2str(round(quality_control/total*100, 1))+ "%",...
    "Materials " + num2str(round(materials/total*100, 1))+ "%"...
    ,"Avionics " + num2str(round(avionics/total*100, 1))+ "%"...
    ,"Powertrain " + num2str(round(powertrain/total*100, 1))+ "%"};
values = [engineering, flight_test, manufacturing, tooling, quality_control,...
    materials, avionics, powertrain];


s = 30
p1 = pie(values,labels)
title("Fully Composite", "FontSize", s)
% saveas("PieComposite.png")


% set(gca, "FontSize", 30)
%% Aluminium


engineering = 56421;
flight_test = 523.8701172;
manufacturing = 230559.8182;
tooling = 32973.19465;
quality_control = 30947.10455;
materials = 24525.16062;
% lg = -7500;
avionics = 70000;
powertrain = 150274;

total = 596226;


labels = {"Engineering " + num2str(round(engineering/total*100, 1)) + "%"...
    ,"Flight Test "+ num2str(round(flight_test/total*100, 1)) + "%"...
    ,"Manufacturing " + num2str(round(manufacturing/total*100, 1))+ "%"...
    ,"Tooling " + num2str(round(tooling/total*100, 1))+ "%"...
    ,"Quality Control " + num2str(round(quality_control/total*100, 1))+ "%",...
    "Materials " + num2str(round(materials/total*100, 1))+ "%"...
    ,"Avionics " + num2str(round(avionics/total*100, 1))+ "%"...
    ,"Powertrain " + num2str(round(powertrain/total*100, 1))+ "%"};
values = [engineering, flight_test, manufacturing, tooling, quality_control,...
    materials, avionics, powertrain];


s = 30
p2 = pie(values,labels)
title("Fully Aluminium", "FontSize", s)

% set(gca, "FontSize", 30)





%% Aluminium body and composite wing


engineering = 69517;
flight_test = 634.8701172;
manufacturing = 283298.8182;
tooling = 40445.19465;
quality_control = 39334.10455;
materials = 27486.16062;
% lg = -7500;
avionics = 70000;
powertrain = 150274;

total = 680992;


labels = {"Engineering " + num2str(round(engineering/total*100, 1)) + "%"...
    ,"Flight Test "+ num2str(round(flight_test/total*100, 1)) + "%"...
    ,"Manufacturing " + num2str(round(manufacturing/total*100, 1))+ "%"...
    ,"Tooling " + num2str(round(tooling/total*100, 1))+ "%"...
    ,"Quality Control " + num2str(round(quality_control/total*100, 1))+ "%",...
    "Materials " + num2str(round(materials/total*100, 1))+ "%"...
    ,"Avionics " + num2str(round(avionics/total*100, 1))+ "%"...
    ,"Powertrain " + num2str(round(powertrain/total*100, 1))+ "%"};
values = [engineering, flight_test, manufacturing, tooling, quality_control,...
    materials, avionics, powertrain];


s = 30
p3 = pie(values,labels)
title("Aluminium Body With Composite Wing", "FontSize", s)

% set(gca, "FontSize", 30)


%% Operational Cost Fully composite

depreciation = 29.89;
interest = 28.23;
pilot = 42;
maintenance = 105;
insurance = 5.6;
landing_fee = 7.33;
fuel_cost = 21.99;
engine_overhaul = 19.99;
electricity = 1.21;
battery_overhaul = 2.745;

total = 264.04;

values = [depreciation, interest, pilot, maintenance,...
    insurance, landing_fee, fuel_cost, engine_overhaul,...
    electricity, battery_overhaul]

labels2 = {"Depreciation " + num2str(round(depreciation/total*100,1)) + "%",...
    "Interest "+ num2str(round(interest/total*100,1)) + "%",...
    "Pilot " + num2str(round(pilot/total*100,1)) + "%",...
    "Maintenance " + num2str(round(maintenance/total*100,1)) + "%",...
    "Insurance "+ num2str(round(insurance/total*100,1)) + "%",...
    "Landing Fee "+ num2str(round(landing_fee/total*100,1)) + "%",...
    "Fuel Cost "+ num2str(round(fuel_cost/total*100,1)) + "%",...
    "Engine Overhaul "+ num2str(round(engine_overhaul/total*100,1)) + "%",...
    "Electricity "+ num2str(round(electricity/total*100,1)) + "%",...
    "Battery Overhaul "+ num2str(round(battery_overhaul/total*100,1)) + "%"}
pie(values, labels2)
title("Fully Composite","FontSize", s)

%% Operational Cost Fully Aluminum

depreciation = 22.25;
interest = 21.02;
pilot = 42;
maintenance = 120;
insurance = 3.87;
landing_fee = 11.4;
fuel_cost = 31.99;
engine_overhaul = 19.99;
electricity = 1.21;
battery_overhaul = 2.745;

total = 275.04;

values = [depreciation, interest, pilot, maintenance,...
    insurance, landing_fee, fuel_cost, engine_overhaul,...
    electricity, battery_overhaul]

labels2 = {"Depreciation " + num2str(round(depreciation/total*100,1)) + "%",...
    "Interest "+ num2str(round(interest/total*100,1)) + "%",...
    "Pilot " + num2str(round(pilot/total*100,1)) + "%",...
    "Maintenance " + num2str(round(maintenance/total*100,1)) + "%",...
    "Insurance "+ num2str(round(insurance/total*100,1)) + "%",...
    "Landing Fee "+ num2str(round(landing_fee/total*100,1)) + "%",...
    "Fuel Cost "+ num2str(round(fuel_cost/total*100,1)) + "%",...
    "Engine Overhaul "+ num2str(round(engine_overhaul/total*100,1)) + "%",...
    "Electricity "+ num2str(round(electricity/total*100,1)) + "%",...
    "Battery Overhaul "+ num2str(round(battery_overhaul/total*100,1)) + "%"}
pie(values, labels2)
title("Fully Aluminum", "FontSize", s)


%% Operational Cost Aluminum + Composite

depreciation = 25.42;
interest = 24.02;
pilot = 42;
maintenance = 115;
insurance = 4.62;
landing_fee = 12.76;
fuel_cost = 28.83;
engine_overhaul = 19.99;
electricity = 1.21;
battery_overhaul = 2.745;

total = 276.04;

values = [depreciation, interest, pilot, maintenance,...
    insurance, landing_fee, fuel_cost, engine_overhaul,...
    electricity, battery_overhaul]

labels2 = {"Depreciation " + num2str(round(depreciation/total*100,1)) + "%",...
    "Interest "+ num2str(round(interest/total*100,1)) + "%",...
    "Pilot " + num2str(round(pilot/total*100,1)) + "%",...
    "Maintenance " + num2str(round(maintenance/total*100,1)) + "%",...
    "Insurance "+ num2str(round(insurance/total*100,1)) + "%",...
    "Landing Fee "+ num2str(round(landing_fee/total*100,1)) + "%",...
    "Fuel Cost "+ num2str(round(fuel_cost/total*100,1)) + "%",...
    "Engine Overhaul "+ num2str(round(engine_overhaul/total*100,1)) + "%",...
    "Electricity "+ num2str(round(electricity/total*100,1)) + "%",...
    "Battery Overhaul "+ num2str(round(battery_overhaul/total*100,1)) + "%"}
pie(values, labels2)
title("Aluminum With Composite Wing", "FontSize", s)