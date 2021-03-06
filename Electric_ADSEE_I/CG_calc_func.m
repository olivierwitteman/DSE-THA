function [x_lemac, most_aft_cg, most_forward_cg] = CG_calc_func(MAC, payload, fus_length, m_fuel, MTOW, OEW, X_oew, X_payload, xc_oewcg, xc_wcg, wing_x, empen_x, fus_x, nacell_x)
% IT doesn not work good for wing mounted engines, but should be fine :D

% fractions for each subsystem based on roskam V page 131 SINGLE ENGINE
% PROPELLER
strucutre_grp = 0.3216;
power_plant_grp = 0.22425;
fixed_eq_grp = 0.11475;
empty_weight_grp = 0.6258;
wing_grp = 0.1074;
empen_grp = 0.025;
fus_grp = 0.1128;
nacelle_grp = 0.01525;


% mass fractions matrix
mass_fractions_subsyst = [wing_grp, empen_grp...
    , fus_grp, nacelle_grp...
    , power_plant_grp, fixed_eq_grp];


% dist form nose
wing_x = wing_x*fus_length;     % <---- INPUT from func input

empen_x = empen_x*fus_length;   % <---- INPUT from func input

fus_x = fus_x*fus_length;       % <---- INPUT from func input

nacell_x = nacell_x*fus_length; % <---- INPUT from func input

propul_x = nacell_x;            % <---- INPUT from func input

fixed_x = fus_x;                % <---- INPUT from func input

% dist matrix
dist_subsyst = [wing_x,
    empen_x,
    fus_x,
    nacell_x,
    propul_x,
    fixed_x];



% moment calculation
moment = mass_fractions_subsyst .* dist_subsyst.';


prompt_prop_pos = 'Where to put the enignes? wing or fuselage . 1/2: ';
prop_pos = double(input(prompt_prop_pos));

if prop_pos == 2    % fuselage
    M_fus_grp = [empen_grp, fus_grp, nacelle_grp, power_plant_grp, fixed_eq_grp];
    X_f = [empen_x, fus_x, nacell_x, propul_x, fixed_x];
    X_fcg = sum(M_fus_grp .* X_f)/(sum(M_fus_grp)); % moments formula
    
    M_wing_grp = wing_grp;
    X_wg = wing_x ;              
end

if prop_pos == 1    % wing
    M_fus_grp = [empen_grp, fus_grp, fixed_eq_grp]
    X_f = [empen_x, fus_x, fixed_eq_grp];           
    X_fcg = sum(M_fus_grp .* X_f)/(sum(M_fus_grp));
    
    
    X_wg = [wing_x, propul_x, nacell_x];         
    M_wing_grp = [wing_grp , power_plant_grp , nacelle_grp];
    X_wg = sum(M_wing_grp.*X_wg)/(sum(M_wing_grp)); % moments formula
end
M_wing_grp = sum(M_wing_grp); % make the group as values for X_lemac formula
M_fus_grp = sum(M_fus_grp);   % make the group as values for X_lemac formula

disp("XXXXXX_fcg")
disp(X_fcg)



% X LEMAC CALC
x_lemac = X_fcg + MAC*(xc_oewcg*M_wing_grp/M_fus_grp - xc_oewcg*(1 + M_wing_grp/M_fus_grp));


X_oew = fus_length * X_oew; % . in the middle for now % <---- INPUT from func input

X_payload = fus_length * X_payload;    % <---- INPUT from func input 

X_fuel = x_lemac + 0.05 * fus_length;  % < --- INPUT 5% of fuselage aft from the X_lemac

% change to xlemac + 0.15*mac



% final table
Mi_fin = [OEW, payload, m_fuel]/MTOW;
Xi_fin = [X_oew, X_payload, X_fuel];
MX_fin = Mi_fin .* Xi_fin;


M_woe_wp = Mi_fin(1) + Mi_fin(2);
M_woe_wp_wf = sum(Mi_fin);
M_woe_wf = Mi_fin(1) + Mi_fin(3);

X_woe_wp = (Xi_fin(1) + Xi_fin(2))/2;
X_woe_wp_wf = (Xi_fin(1) + Xi_fin(2) + Xi_fin(3))/3;
X_woe_wf = (Xi_fin(1) + Xi_fin(3))/2;

% CG EXCURSION
most_aft_cg = max([X_woe_wp,X_woe_wp_wf,X_woe_wf, X_oew]);
most_forward_cg = min([X_woe_wp,X_woe_wp_wf,X_woe_wf, X_oew]);

f5 = figure
plot([X_woe_wp] , [M_woe_wp],'-k*', 'MarkerSize',20)
hold on
plot([X_woe_wp_wf] , [M_woe_wp_wf],'-b*', 'MarkerSize',20)
plot([X_woe_wf] , [M_woe_wf],'-r*', 'MarkerSize',20)
plot([X_oew] , [OEW/MTOW],'-g*', 'MarkerSize',20)
plot([x_lemac, x_lemac + MAC], [0.6, 0.6], "LineWidth", 20); % airfoil
ylim([0.55, 1.06])


legend({'OEW + WP'...
    ,'OEW + WP + WF'...
    ,'OEW + WF'...
    ,'OEW'...
    , "MAC"
},'Location','northeast', 'FontSize', 16)


f6 = figure;
weights_for_pie = [payload, OEW, m_fuel];
perc_weights = weights_for_pie/MTOW;
labels = {"Payload " + string(perc_weights(1)),...
    "OEW " + string(perc_weights(2)),...
    "Fuel "+ string(perc_weights(3))}
pie(weights_for_pie, labels)
end