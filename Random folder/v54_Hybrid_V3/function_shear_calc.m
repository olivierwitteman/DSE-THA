%% function
%% Calculate the shear and stresses in the booms and the skin of the wing box
%A few assumptions are made in this program. Only the booms carry he normal
%stresses and the skin panels only carry the shear stress.
%The general inputs are shown below.
%S_Y is at the boom 2. the same yields for
%% INPUTS
function[shear_flow, shear_stress, max_shear, A ]= function_shear_calc(c, S_x, S_y, I_xx, I_yy, I_xy, t_box, A, Torsion)
                    %INPUT [m]      Chord length      tip=0.64 root=1.6 bekend
%syms A
%I_xx=48*10^-3;                %INPUT [m^4]    Moment of inertia around x axis
%I_yy=52*10^-3;                %INPUT [m^4]    Moment of inertia around y axis
%I_xy=20*10^-3;                %INPUT ]m^4]    Product moment of inertia.
%s_top=1.141/1.96^c;            %INPUT [m]      distance between each boom
%s_bottom=1.138/1.96*c;         %INPUT [m]      distance at the bottom of each boom

A_encl2=0.1.*c.*0.58.*c;               %0.155
%S_x_pos=[S_x_x 0.0804]*c;       %INPUT [m]    location of the shear force in x-direction
%S_y_pos=[0.295 S_y_y]*c;       %INPUT [m]    location of the shear force in y-direction
boom1=A;                   %INPUT [m^2]   Boom area of boom1
boom2=A;                   %INPUT [m^2]   Boom area of boom2
boom3=A;                   %INPUT [m^2]   Boom area of boom3
boom4=A;                   %INPUT [m^2]   Boom area of boom4
boom5=A;                   %INPUT [m^2]   Boom area of boom5
boom6=A;                   %INPUT [m^2]   Boom area of boom6
boom7=A;                   %INPUT [m^2]   Boom area of boom7
boom8=A;                   %INPUT [m^2]   Boom area of boom8
boom9=A;                   %INPUT [m^2]   Boom area of boom9
boom10=A;                  %INPUT [m^2]   Boom area of boom10
boom1_xy=[0.15 0.065];                  %INPUT [m]     Boom position of boom 1 [x y]
boom2_xy=[0.295 0.0804];                %INPUT [m]     Boom position of boom 2 [x y]
boom3_xy=[0.44 0.0755];                 %INPUT [m]     Boom position of boom 3 [x y]
boom4_xy=[0.585 0.0654];                %INPUT [m]     Boom position of boom 4 [x y]
boom5_xy=[0.73 0.0473];                 %INPUT [m]     Boom position of boom 5 [x y]
boom6_xy=[0.73 -0.0197];                %INPUT [m]     Boom position of boom 6 [x y]
boom7_xy=[0.585 -0.0288];               %INPUT [m]
boom8_xy=[0.44 -0.036];                 %INPUT [m]     Boom position of boom 8 [x y]
boom9_xy=[0.295 -0.0414];               %INPUT [m]     Boom position of boom 9 [x y]
boom10_xy=[0.15 -0.0412];                 %INPUT [m]     Boom position of boom 10 [x y]
im=[1 2 3 4 5 6 7 8 9 10];               %for programs sake, nothing to worry about.
%% Calculation for basic shear flow.
% the cut is made in 1-2  Then the basic shear flow is calculated by going
% from 1-2 and going clockwise to 10-1  inner (index 10)
Boom_area=[ boom1 boom2 boom3 boom4 boom5 boom6 boom7 boom8 boom9 boom10];
boom_pos=[boom1_xy; boom2_xy; boom3_xy; boom4_xy; boom5_xy; boom6_xy; boom7_xy; boom8_xy; boom9_xy; boom10_xy]*c;
c_g= [0.44 0.0267]*c;
q_left=-(S_x*I_xx-S_y*I_xy)/(I_xx*I_yy-I_xy^2);         %first term of the basic shear flow.
q_right=-(S_y*I_yy-S_x*I_xy)/(I_xx*I_yy-I_xy^2);        %second term for the basic shear flow
q_b1=zeros(10,1); %leave it
length(q_b1)
boom_pos_cg=zeros(10,2);
for num=im
    boom_pos_cg(num,1)=boom_pos(num,1)-c_g(1);
    boom_pos_cg(num,2)= boom_pos(num,2)-c_g(2);
end
im=[2 3 4 5 6 7 8 9 10];
s=2;

%disp(q_left)
%disp(q_right)
%disp(boom_pos_cg(2,1)*Boom_area(1))
qvaluemat=[];
for number=im
    qvalue=q_left*Boom_area(number)*boom_pos_cg(number,1)+q_right*Boom_area(number)*boom_pos_cg(number,2);
    qvaluemat=[qvaluemat, qvalue];
    qvalue1= q_b1(s-1,1)+qvalue;
    q_b1(s)=[qvalue1] ;      %the basic shear flow, starting from 10-1 outer to 10-1 innter
    length(q_b1)
    s=s+1;
end
%qvaluemat
disp(500000)
%disp(q_b1)
%disp(q_left*Boom_area(3)*boom_pos_cg(3,1));%+q_right*Boom_area(3)*boom_pos_cg(3,2))
%% Calculate the final shear flow.
%this is done by closing the cut and then adding it to the basic shear flow
%that is calculated above. Take the moment around boom 3. This makes shear
%flow q45 and q56 redundant as they do not have a moment around this point.

j=[1 2 3 4 5 6 7 8 9 ];
m=0;
d_relative=zeros(10,2);
for im=j

   d_x = boom_pos_cg(im+1,1)-boom_pos_cg(im,1);
   d_y = boom_pos_cg(im+1,2)-boom_pos_cg(im,2);
   d_relative(im,1)=d_x;
   d_relative(im,2)=d_y;
   d_relative(10,1)=boom_pos_cg(1,1)-boom_pos_cg(10,1);
   d_relative(10,2)=boom_pos_cg(1,2)-boom_pos_cg(10,2);
end
kut=[1 2 3 4 5 6 7 8 9 10];
d_2_relative=zeros(10,2);
for kut1=kut
    d_2_relative(kut1,1)=+boom_pos_cg(kut1,1)-boom_pos_cg(2,1);
    d_2_relative(kut1,2)= -boom_pos_cg(2,2)+boom_pos_cg(kut1,2);
end
%disp(d_2_relative)
q_b2=q_b1;


%q_b2(1)=[];
kutmoment=zeros(10,1);

for kut2=kut
    moment=-q_b2(kut2)*d_relative(kut2,1)*d_2_relative(kut2,1)+q_b2(kut2)*d_relative(kut2,2)*d_2_relative(kut2,2);
    %moment=-q_b2(kut2)*d_5_relative(kut2,1)+q_b2(kut2)*d_5_relative(kut2,2)
    kutmoment(kut2)=moment;
end
%Shear_moment=-S_x*(boom_pos(5,2)-S_x_pos(2))+S_y*(boom_pos(5,1)-S_y_pos(1));
q01moment=2*A_encl2;
q_s_shear=sum(kutmoment)/q01moment;
q_s_torsion=Torsion/(A_encl2.*2);
shear_flow=q_b1+q_s_torsion+q_s_shear;

shear_stress=shear_flow;
max_shear=max(abs(shear_stress));

shear_stress=shear_flow/t_box;
max_shear=max(abs(shear_stress));
disp(q_b1)

end
%% Rate of twist. Ignore below
%now that the basic shear flow is done, the next step in calculating the
%rate of twist in each section and then generate a system of equations and
%then calcualte the resilliant shear flow.
%s_front=boom_pos_cg(1,2)-boom_pos_cg(10,2);       %INPUT [m]      length of the front spar
%s_back =boom_pos_cg(5,2)-boom_pos_cg(6,2);      %INPUT [m]      length of the rear spar

%d_top=s_top/4.;

%d_bottom=s_bottom/4.;

%twist_semicirc_q01=1/(2*A_encl1)*(semi_circ/t_semicirc+s_front/t_box);
%twist_semicirc_q02=1/(2*A_encl1)*(-s_front/t_box);
%twist_semicric_rest=1/(2*A_encl1)*(-q_b1(10)*s_front)/t_box;

%twist_box_rest=1/(2*A_encl2)*(q_b1(2)+q_b1(3)+q_b1(4)+q_b1(5))*d_top+1/(2*A_encl2)*q_b1(6)*s_back+1/(2*A_encl2)*(q_b1(7)+q_b1(8)+q_b1(9)+q_b1(10))*d_bottom+1/(2*A_encl2)*q_b1(11)*s_front;
%twist_box_q01=1/(2*A_encl2)*(-s_front/t_box);
%twist_box_q02=1/(2*A_encl2)*((s_front+s_top+s_bottom+s_back)/t_box);
%% The final shear flow throughout the cross section.
%solvmatrix1=[q01moment q02moment ; twist_semicirc_q01-twist_box_q01 twist_semicirc_q02-twist_box_q02];
%solvmatrix2=[Shear_moment-sum(kutmoment);  twist_box_rest-twist_semicric_rest];
%q0mat=linsolve(solvmatrix1,solvmatrix2);
%disp(sum(kutmoment))
