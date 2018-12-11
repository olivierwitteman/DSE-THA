clc; clear all
cb=680000; %specific density  in J/kg 188.89 Wh
V_TO=60;
rho=1.225;
MTOW=1750;
A=8.2;
b=10.5948;
S= b.^2/A;
e=0.75;
cD_0=0.0287; %cD_0 take-off from Roelof 
cL_TO=2.1;
cD=cD_0 + (cL_TO.^2)/(pi*A*e);

p_req_TO = 0.5*rho*S*cD*V_TO.^3

P=0.61*10.^6 %Watt


P_req_TO_Bat = p_req_TO*0.2 %20% of TO power in Watt
Specific_Power_lithium_ion = 500 %W/kg
Battery_mass = P/Specific_Power_lithium_ion






