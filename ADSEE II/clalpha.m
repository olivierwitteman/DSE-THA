function [CLmax, alpha_stall] = clalpha(A, clmax, sweepLE, CLdes, S_ref)
%input values for Cl-alpha curve for finite wing (clean and flapped
%General inputs

sweepc2 = 0.25;
Beta = 1;
eta = 0.95;
%%Airfoil ans wing

alpha0L = 0; %Given by the airfoil
taper = 0.6; %Given
Sharpfactor = 0.02; %from the airfoil.
%%statistical
%general
C1 = 0.02; %From Table (slide:15)
C2 = 0.02;  %From table slide: 15 raymer.


CLmax_base = 0.2;
% Use the datcom method to get the CL slope from statistics
Datcomtop = 2 * pi * A;
Datcombottom=sqrt(4 + (A * Beta / eta)^2 * (1 + ((tan(sweepc2))^2/Beta^2))) + 2;
CLalpha = Datcomtop / Datcombottom;

%next is the trim angle.this is the angle the aircraft needs to fly at to
%fly at CLdes

alphatrim = CLdes / CLalpha + alpha0L;

%Next is the CLmax, 2 methods are being introduced. The Datcom method and a
%general  one. Datcom is prefferred. CLmax

Datcom_choose_method=((C1+1)*cos(sweepLE));

delta_CLmax = 0.5; %Term for M>0.2. Not for take-off and landing.
delta_alpha_CLmax = 0.02; %From table slide 19 (Raymer)
CLM_clmax = 0.05; %slide 17 use right table
CLmax_base = 0.2; %table from Raymer
alpha_CLmax = 0.04;

if A> 4/Datcom_choose_method
    delta_CLmax_Datcom=0.5; %Term for M>0.2. Not for take-off and landing.
    delta_alpha_CLmax=0.02; %From table slide 19 (Raymer)
    CLM_clmax=0.05; %slide 17 use right table
    %high aspect ratio method
    CLmax=CLM_clmax*clmax+delta_CLmax_Datcom;
    alpha_stall=CLmax/CLalpha+alpha0L+delta_alpha_CLmax;
%otherwise low aspect ratio method.
elseif A< 4/Datcom_choose_method
    CLmax_base=0.2; %table from Raymer
    delta_CLmax_Datcom=0.5;%table raymer for compressibility
    delta_alpha_CLmax=0.02; %table Raymer
    alpha_CLmax=0.04;%table Raymer

    CLmax=CLmax_base+delta_CLmax_Datcom;
    alpha_stall=alpha_CLmax+delta_alpha_CLmax;
end

%% High Lift Devices.
%Can we take-off at CLmax-clean? Hopfully no otherwise redesign wing.
%Propose approiate type and combination of HLD's
%estimate wing area Swf/S needs to be known. This can also be done another
%way by using a Known DeltaCLMAX and calculating the needed area to a
%feasibility check.
%assume chord fractions (cf/c) (slide 52)
%Use the following equations
sweep_hinge = 0.025; %hinge line for the flap/slat

Swf = 2; %Area affected by the HLD might be easier to do this other way.

delta_alpha0L_airfoil = 5;%Get from the LE or TE difficult to say.
LE = 'Slat'; %enter the leading edge flap according to the name
TE ='tripple_slotted'; %enter TE flap according to name;
c_ac_c = 0.5; %c'/c is a something we choose in general.
Delta_Clmax_HLD = HLD(LE,TE,c_ac_c);

HLD_Delta_CLmax = 0.9 * Delta_Clmax_HLD * Swf / S_ref * cos(sweep_hinge); %The difference in CLmax
%due to the high lift devices.
HLD_Delta_alpha0L = delta_alpha0L_airfoil * Swf / S_ref*cos(sweep_hinge); %Same but than
%0L angle of attack.
%%
%CLalpha for flapped configuration
%if HLD increase wing surface
S_ac_S = 1 + Swf/S_ref * (c_ac_c - 1);% S'/S
CL_alpha_flap = 1.; %The value is the same if certain flaps are not used.
if S_ac_S > 1.
    CL_alpha_flap = S_ac_S * CLalpha;
end

%Check if in the end the values make sense. Feasability check is very
%important.
