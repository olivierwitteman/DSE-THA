clc; clear all

%Maneuver and Gust Diagrams 


%% Initial inputs
g = 9.81; %gravitational acceleration [m/s^2]
m = 1435.0; %total mass of aircraft [kg]
MTOW = m*g; %maximum take-off weight [N]
H = 2400.; %cruise altitude [m]
R = 287.; %Gas constant for air [m2/s2/K]
rho = 0.9664; %air density at sea level[kg/m^3]
L = -0.0065; %temperature lapse rate at troposphere [deg/m]
T = 288.15; %sea level temperature [K]
T_alt = T+L*H; %temperatute at chosen altitude [K]
rho_alt = rho*(T_alt/T)^(-(g/(L*R)+1)); %density at chosen altitude [kg/m^3]
Cl_max = 2.1; %max lift coefficient of aircraft [-]
Cl = 0.355; %lift coefficient at cruise [-], 0 deg angle of attack, -2 deg incidence
S = 9.52; %wing surface area [m^2]
n_max = 4.4; %maximum loading factor (CS23.333)
n_min = -0.4*n_max; % minimum loading factor (CS23.333)

%% V-n Diagram - axis and square curve

V_S = sqrt((2*MTOW)/(Cl_max*rho*S)); % stall speed (EAS) [m/s]
V_S = 33.0;
V_A = sqrt((n_max*2*MTOW)./(rho*Cl_max*S));  %design maneuver speed at max loading factor (EAS) [m/s] 
% V_A = V_S*sqrt(n_max);
V_A = 0:0.1:V_A; % vector for the first curve up until V_A
n = (V_A.^(2)*rho*Cl_max*S)./(2*MTOW);   % vector for load factors corresponging to V_A_eas


plot(V_A,n,'Color','blue')
hold on
line([0,270],[0,0],'LineStyle','--','Color','r') 

%% V-n Diagram - Horizontal upper part

V_C = sqrt((2*MTOW*n_max)/(rho*S*Cl)); %cruise speed at max loading factor (EAS)[m/s]

V_D = V_C*1.45 %design dive speed (CS23.333), (EAS) [m/s]

line([V_A(end), V_D],[n_max,n_max],'Color','blue') %horizontal max laod
line([V_D,V_D],[0,n_max],'Color','blue')


%% V-n Diagram - NEGATIVE PART

V_S_neg = sqrt((2*abs(n_min)*MTOW)./(rho*Cl_max*S)); %negative stall speed (TAS) [m/s]

line([V_D,V_C],[0,n_min],'Color','blue'); %negative slope line for negative load
line([V_C,V_S_neg],[n_min,n_min],'Color','blue'); %negative horizontal line for constant minimum load factor

V_S_neg = V_S_neg:-0.01:0;
n_neg = -(V_S_neg.^(2)*rho*Cl_max*S)./(2*MTOW);
plot(V_S_neg,n_neg,'Color','blue'); hold on

xlim([0,V_D+10])              
ylim([-6,7.0])            
ylabel('n [-]','fontweight','bold','fontsize',14)
xlabel('V_{EAS} [m/s]','fontweight','bold','fontsize',14)
% title('V-n Maneuver and Gust Diagram','fontweight','bold','fontsize',14)

%% Plot speeds and lines

line([V_S,V_S],[0,1],'LineStyle','--','Color','Magenta') %line for stall speed
line([0,V_S],[1,1],'LineStyle','--','Color','Magenta') %line for stall speed
plot(V_S,0,'.','MarkerSize',20,'Color','Magenta'); %dot for stall speed
x1 = V_S+1; %x location of text for stall speed
y1 = 0.3; %y location of text for stall speed
txt1 ='V-S'; %for stall speed
text(x1,y1,txt1,'FontSize',18) %for stall speed

line([V_A(end),V_A(end)],[0,n_max],'LineStyle','--','Color','blue') %line for maneuver speed
plot(V_A(end),0,'Marker','.','MarkerSize',20,'Color','blue');  %dot for maneuver speed
x1 = V_A(end)+1; %x location of text for maneuver speed
y1 = 0.3; %y location of text for maneuver speed
txt1 = 'V-A'; %for maneuver speed
text(x1,y1,txt1,'FontSize',18) %for maneuver speed

line([V_C,V_C],[n_min,n_max],'LineStyle','--','Color','green') %line for cruise speed
plot(V_C,0,'Marker','.','MarkerSize',20,'Color','green'); %dot for cruise speed
x1 = V_C+1; %x location of text for cruise speed
y1 = 0.3; %y location of text for cruise speed
txt1 = 'V-C'; %for cruise speed
text(x1,y1,txt1,'FontSize',18) %for cruise speed

plot(V_D,0,'Marker','.','MarkerSize',20,'Color','cyan'); %dot for dive speed
x1 = V_D+1; %x location of text for dive speed
y1 = 0.3; %y location of text for dive speed
txt1 = 'V-D'; %for dive speed
text(x1,y1,txt1,'FontSize',18) %for dive speed

%% Addition for Gust

WS = MTOW/S; %wing loading [N/m^2]
Cla = 6.76; %lift curve slope [1/rad] 4.76
mac = 1.036; %mean aerodynamic chord [m]
miu = (2*WS)/(rho_alt*mac*Cla*g); %aircraft mass ratio [-]
Kg = (0.88*miu)/(5.3+miu); %gust alleviation factor [-]
Ub = 20.12; %gust velocity for V_B from CS23.333 [m/s]
Uc = 15.24; %gust velocity for V_C from CS23.333 [m/s]
Ud = 7.62; %gust velocity for V_D from CS23.333 [m/s]
V_B = V_S*sqrt(1+((Kg*Ub*V_C*Cla)/(WS))); %bad weather velocity [m/s]

ng_b = 1+(Kg*rho*Ub*V_B*Cla)/(2*WS); %positive change in load factor due to gust in V_B
ng_b_neg = 1-(Kg*rho*Ub*V_B*Cla)/(2*WS); %negative change in load factor due to gust in V_B
ng_c = 1+(Kg*rho*Uc*V_C*Cla)/(2*WS); %positive change in load factor due to gust in V_C
ng_c_neg = 1-(Kg*rho*Uc*V_C*Cla)/(2*WS); %negative change in load factor due to gust in V_C
ng_d = 1+(Kg*rho*Ud*V_D*Cla)/(2*WS); %positive change in load factor due to gust in V_D
ng_d_neg = 1-(Kg*rho*Ud*V_D*Cla)/(2*WS); %negative change in load factor due to gust in V_D

line([0,V_B],[1,ng_b],'LineStyle','--','Color','black') % positive gust line for dive speed at Ub = 20.12 and V_B=74.32
line([0,V_B],[1,ng_b_neg],'LineStyle','--','Color','black') % negative gust line for dive speed at Ub = 20.12 and V_B=104.0585
x1 = V_B+1-10*5.5; %x location of text for dive speed
y1 = ng_b-0.2; %y location of text for dive speed
txt1 = 'U_{b} = 20.12 [m/s]'; %for dive speed
text(x1,y1,txt1,'FontSize',18) %for dive speed
line([0,V_C],[1,ng_c],'LineStyle','--','Color','black') % positive gust line for cruise speed at Uc = 15.24 and V_C=74.32
line([0,V_C],[1,ng_c_neg],'LineStyle','--','Color','black') % negative gust line for cruise speed at Uc = 15.24 and V_C=74.32
x1 = V_C+1 - 10*3; %x location of text for dive speed
y1 = ng_c-0.2; %y location of text for dive speed
txt1 = 'U_{c} = 15.24 [m/s]'; %for dive speed
text(x1,y1,txt1,'FontSize',18) %for dive speed
line([0,V_D],[1,ng_d],'LineStyle','--','Color','black') % positive gust line for dive speed at Ud = 7.62 and V_D=74.32
line([0,V_D],[1,ng_d_neg],'LineStyle','--','Color','black') % negative gust line for dive speed at Ud = 7.62 and V_D=104.0585
x1 = V_D-1-25; %x location of text for dive speed
y1 = ng_d+0.5+0.2; %y location of text for dive speed
txt1 = 'U_{d} = 7.62 [m/s]'; % for dive speed
text(x1,y1,txt1,'FontSize',18) % for dive speed

set(gca, "FontSize", 20)