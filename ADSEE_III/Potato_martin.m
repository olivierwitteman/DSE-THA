function xyz = Potato_martin()

clc; clear variables
mass_passenger=94.56;    %kg INPUT
mass_cargo=1000;         %kg INPUT
mass_fuel=4000;          %kg INPUT
OEW=11501;               %kg
g=9.81;
X_lemac=12.594;          %TIM INPUT
MAC=2.865;               %TIM INPUT


row1=6.515;              %Pos row 1                     %INPUT ot kartinkite
norm_seat=0.7874;        %spacing between normal rows   %INPUT (wdieh go u saita)
spec_seat=0.8636 ;       %spacing between 8 and 9       %INPUT ot sait

rows_y=[];
rows_y=[rows_y, row1];
r1=1:1:7;               %length till special seat
rows_y=[rows_y, row1+r1*norm_seat];

rows_y=[rows_y, rows_y(end)+spec_seat]; %special seat position

r2=1:1:4;
rows_y=[rows_y, rows_y(end)+r2*norm_seat];
rows_y_left=rows_y;                        %final rows LEFT SIDE Passengers
rows_y_right=rows_y([2:end]);              %final rows RIGHT SIDE Passengers

ypos_total=[];
ypos_total=[ypos_total, rows_y_left]    %all passengers poston in y

ypos_mat=zeros(13,3);                    %create matrix for overwriting
ypos_mat([1:end],1)=ypos_total;          %left column
ypos_mat([2:end],2)=ypos_total([2:end]); %middle column
ypos_mat([1:end],3)=ypos_total;          %right column
%% Moments
arm_OEW=X_lemac+0.3*MAC;              %ASSUMED from page 13 30%MAC from table   %%INPUT
arm_OEW1=100*(arm_OEW-X_lemac)/MAC;
moment_OEW=arm_OEW*OEW*g

arm_cargo=(12.43+24.39)./2; %From embreared documents average Page 13           %%INPUT sredno arimeti4no ot kartinkite mejdu tail i posledniq stol
moment_cargo=arm_cargo*mass_cargo*g

ZF_moment=moment_cargo+moment_OEW;
ZF_weight=(OEW+mass_cargo)*g;
ZF_cg=ZF_moment./ZF_weight
ZF_cg1=100*(ZF_cg-X_lemac)/MAC;

%% momenz FRONT TO BACK
ypos_window=zeros(13,2);                    %create matrix for overwriting
ypos_window([1:end],1)=ypos_total;          %left column
ypos_window([1:end],2)=ypos_total;          %right column

weight_rows_windows=zeros(13,2);
weight_rows_windows([1:end],1)=mass_passenger*g
weight_rows_windows([1:end],2)=mass_passenger*g
weight_rows_windows([1:1],2)=0; %zero spots
weight_rows_windows2=[]; 
for i=1:1:13 %Sum two columns
    weight_rows_windows2=[weight_rows_windows2,...
        weight_rows_windows([i],1)+weight_rows_windows([i],2)];
end
for i=2:1:13  %sum rows
    weight_rows_windows2([i])=...
        weight_rows_windows2([i])+weight_rows_windows2([i-1])
end
weight_rows_windows2=weight_rows_windows2+ZF_weight;

moment_frontback_window=zeros(2,13);
moment_frontback_window=weight_rows_windows.*ypos_window; %arm times weight
moment_frontback_window=moment_frontback_window; %PLUS ZF moment
moment_frontback_window2=[];
for i=1:1:13
    moment_frontback_window2=[moment_frontback_window2,...
        moment_frontback_window([i],1)+moment_frontback_window([i],2)];
end
for i=2:1:13
    moment_frontback_window2([i])=...
        moment_frontback_window2([i])+moment_frontback_window2([i-1]);
end
moment_frontback_window2=moment_frontback_window2+ZF_moment
x_cgs=moment_frontback_window2./weight_rows_windows2;
x_cgs=100*(x_cgs-X_lemac)/MAC;
scatter(x_cgs, weight_rows_windows2/g)
hold on
plot(x_cgs, weight_rows_windows2/g,'Color','[0.5843    0.8157    0.9882]')
line([arm_OEW1, ZF_cg1],[OEW*g/g,ZF_weight/g],'Color','blue')
line([ZF_cg1, x_cgs(1) ],[ZF_weight/g,weight_rows_windows2(1)/g ],'Color','[0.5843    0.8157    0.9882]')

%% Moements ISLE FRON TO BACK
ypos_isle=zeros(13,1);                   %create matrix for overwriting
ypos_isle([1:end],1)=ypos_total;         %left column

weight_rows_isle=zeros(13,1);
weight_rows_isle([1:end],1)=mass_passenger*g
weight_rows_isle([1:1],1)=0; %zero spots
weight_rows_isle2=[]; 
for i=1:1:13 %Sum two columns
    weight_rows_isle2=[weight_rows_isle2,...
        weight_rows_isle([i],1)];
end
for i=2:1:13  %sum rows
    weight_rows_isle2([i])=...
        weight_rows_isle2([i])+weight_rows_isle2([i-1])
end

weight_rows_isle2=weight_rows_isle2+weight_rows_windows2(end);

moment_frontback_isle=zeros(1,13);
moment_frontback_isle=weight_rows_isle.*ypos_isle; %arm times weight
moment_frontback_isle2=[];

for i=2:1:13
    moment_frontback_isle([i])=...
        moment_frontback_isle([i])+moment_frontback_isle([i-1]);
end
moment_frontback_isle2=moment_frontback_isle+moment_frontback_window2(end);
x_cgs=moment_frontback_isle2.'./weight_rows_isle2;
x_cgs=100*(x_cgs-X_lemac)/MAC;
scatter(x_cgs, weight_rows_isle2/g)
plot(x_cgs, weight_rows_isle2/g,'Color','yellow')
hold on

%% BACK TO FRONT
% back to fron windows
ypos_window=zeros(13,2);                   %create matrix for overwriting
ypos_window([1:end],1)=ypos_total;         %left column
ypos_window([1:end],2)=ypos_total;          %right column
ypos_window=flipud(ypos_window) %%%%%%%%%%% to mirror plot

weight_rows_windows=zeros(13,2);
weight_rows_windows([1:end],1)=mass_passenger*g
weight_rows_windows([1:end],2)=mass_passenger*g
weight_rows_windows([1:1],2)=0; %zero spots
weight_rows_windows=flipud(weight_rows_windows) %%%%%%%%

weight_rows_windows2=[]; 
for i=1:1:13 %Sum two columns
    weight_rows_windows2=[weight_rows_windows2,...
        weight_rows_windows([i],1)+weight_rows_windows([i],2)];
end
for i=2:1:13  %sum rows
    weight_rows_windows2([i])=...
        weight_rows_windows2([i])+weight_rows_windows2([i-1])
end
weight_rows_windows2=weight_rows_windows2+ZF_weight;

moment_frontback_window=zeros(2,13);
moment_frontback_window=weight_rows_windows.*ypos_window; %arm times weight
moment_frontback_window=moment_frontback_window; %PLUS ZF moment
moment_frontback_window2=[];
for i=1:1:13
    moment_frontback_window2=[moment_frontback_window2,...
        moment_frontback_window([i],1)+moment_frontback_window([i],2)];
end
for i=2:1:13
    moment_frontback_window2([i])=...
        moment_frontback_window2([i])+moment_frontback_window2([i-1]);
end
moment_frontback_window2=moment_frontback_window2+ZF_moment
x_cgs=moment_frontback_window2./weight_rows_windows2;
x_cgs=100*(x_cgs-X_lemac)/MAC;
scatter(x_cgs, weight_rows_windows2/g)
plot(x_cgs, weight_rows_windows2/g,'Color','green')
line([ZF_cg1, x_cgs(1) ],[ZF_weight/g,weight_rows_windows2(1)/g ],'Color','green')
hold on

%% BAck to Front ISLE
ypos_isle=zeros(13,1);                   %create matrix for overwriting
ypos_isle([1:end],1)=ypos_total;         %left column
ypos_isle=flipud(ypos_isle);

weight_rows_isle=zeros(13,1);
weight_rows_isle([1:end],1)=mass_passenger*g
weight_rows_isle([1:1],1)=0; %zero spots
weight_rows_isle=flipud(weight_rows_isle)

weight_rows_isle2=[]; %TRQBWA POSLEDNIQ RED DA SE DOBAWQ KYM GORNIQ
for i=1:1:13 %Sum two columns
    weight_rows_isle2=[weight_rows_isle2,...
        weight_rows_isle([i],1)];
end
for i=2:1:13  %sum rows
    weight_rows_isle2([i])=...
        weight_rows_isle2([i])+weight_rows_isle2([i-1])
end
weight_rows_isle2=weight_rows_isle2+weight_rows_windows2(end);

moment_frontback_isle=zeros(1,13);
moment_frontback_isle=weight_rows_isle.*ypos_isle; %arm times weight
moment_frontback_isle2=[];

for i=2:1:13
    moment_frontback_isle([i])=...
        moment_frontback_isle([i])+moment_frontback_isle([i-1]);
end
moment_frontback_isle2=moment_frontback_isle+moment_frontback_window2(end);
x_cgs=moment_frontback_isle2.'./weight_rows_isle2;
x_cgs=100*(x_cgs-X_lemac)/MAC;
scatter(x_cgs, weight_rows_isle2/g)
plot(x_cgs, weight_rows_isle2/g,'Color','red')
line([30.57, 31.13 ],[1.487*10.^(4),1.496*10.^(4) ],'Color','red') %%%REDDDD
hold on

%Fuel Part
arm_fuel=13.6;                   %assumed from page 13 from wings           %INPUT
MTOW=weight_rows_isle2(end)/g+mass_fuel+0; %365                                %correction, not needed
moment_MTOW=MTOW/arm_fuel;
arm_fuel1=100*(arm_fuel-X_lemac)/MAC;
line([x_cgs(end), arm_fuel1],[weight_rows_isle2(end)./g, MTOW],'Color','black')
title('Potato Diagram')
ylabel('mass [kg]')
xlabel('x position as % of MAC')
legend('Window','OEW to ZFW','WIndow','Isle','Isle','Payload to RM','')

xyz = "Martin's potato done";
end
