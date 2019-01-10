%% Calculate the shear and stresses in the booms and the skin of the wing box 
%A few assumptions are made in this program. Only the booms carry he normal
%stresses and the skin panels only carry the shear stress.
%The general inputs are shown below.
I_xx=48*10^6;                %INPUT [m^4]    Moment of inertia around x axis
I_yy=52*10^6;                %INPUT [m^4]    Momtet of inertia around y axis
circum=5.;                   %INPUT [m]      Total of top and bottom side of 
n_booms=26.;                 %INPUT          number of booms
spacing_boom= circum/n_booms;%INPUT [m]      %spacing between each booms.
Boom_area=[5 6 2 4 4; 4 5 4 6 5]

CG=[5 5];                 %Calculated and input [x position and y position]
boom_pos=zeros(26,2);
boom_pos(2,2)=25
y=0.2;                       %INPUT [m]     %y spacing of the 
start=1.;                    %Leave this just for codes sake.

boom_pos(
