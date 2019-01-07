<<<<<<< HEAD
function [CLmax_clean, alpha_stall] = clalpha(A, clmax, CLdes, sweep_c2, M, sweepLE, t_c)
=======
function [CLmax_clean, alpha_stall,CL_alpha_clean] = clalpha(A, clmax, CLdes, sweep_c2, M, sweepLE, t_c)
>>>>>>> 8e04a7f2bbf68aa5db112dc28670559be59188f5
%input values for Cl-alpha curve for finite wing (clean and flapped
%General inputs
Beta = sqrt(1-M^2);                  %Compressibility factor for the mach number.
eta = 0.95;             %<-----INPUT, but you cannot change it as we do not have flight data
%AIRFOIL INPUTS
alpha0l = -2.5;         %<-----INPUT, form the airfoil
Sharpfactor = 26*t_c    %<----INPUT Used for some graphs in Raymer.
%INPUTS FROM RAYMER TABLES
C1 = 0.4;               %<-----INPUT: From Table Raymer (slide:15)
C2 = 1.1;               %<-----INPUT  From Table Raymer:(slide:15)
%% Datcom Method for calculating the lift curve slope and CLmax clean and the stall angle.
Datcomtop = 2 * pi * A;         %Top part of the lift curve slope equation 
Datcombottom=sqrt(4 + (A * Beta / eta)^2 * (1 + ((tan(sweep_c2))^2/Beta^2))) + 2; %bottom part of the lift curve slope equation.
CL_alpha_clean = Datcomtop / Datcombottom *pi/180  %The lift curve slope in 1/degrees.


alphatrim = CLdes / CL_alpha_clean + alpha0l %The trim angle, the angle at which cldes is reached.

Datcom_choose_method=((C1+1)*cos(sweepLE));  %WHICH METHOD TO USE: HIGH OR LOW ASPECT RATIO METHOD
if A> 4/Datcom_choose_method
    delta_CLmax_Datcom=-0.2;    %<----INPUT: Term for M>0.2. Not for take-off and landing.
    delta_alpha_CLmax=2.1;       %<----INPUT: Difference for where CLmax,clean occurs. From Raymer table (slide 19).
    CLMax_clmax=0.9;            %<----INPUT: Conversion from 2D to 3D. Input from RAYMER table.
    CLmax_clean=CLMax_clmax*clmax+delta_CLmax_Datcom    %CLmax in clean configuration
    alpha_stall=CLmax_clean/CL_alpha_clean+alpha0l+delta_alpha_CLmax  %stall angle in clean configuration.


%The low aspect ratio method.
elseif A< 4/Datcom_choose_method
    
<<<<<<< HEAD
    disp('something went wrong'
=======
    disp('something went wrong')
>>>>>>> 8e04a7f2bbf68aa5db112dc28670559be59188f5
end

%% High Lift Devices.
%We can go 2 ways with high lift devices, calculate the maximum CLmax by
%assuming a Swf_S and then continue, or by filling in the CLmax and then
%calculating the space necessary for the CLmax.
sweep_hinge = 10*pi/180;            %<----INPUT: The hinge angle of the FLAP
sweep_hinge= sweep_hinge*pi/180.;   %CHANGE DEGREES TO RADIANS.
%% Here make a choice about which way you want to go!!!!!!!
Swf_S = 0.4;                        %<-----(POSSIBLE INPUT): Ratio of the area affected by the HLD
%HLD_Delta_CLmax=0.5;                %<-----(POSSIBLE INPUT): Necessary deltaCLMax needed to be provieded by the HLD.
%% Calculation for the final CLMAX
delta_alpha0L_airfoil = -10;        %<----INPUT: Shift of 0 lift angle. First approxamation: -10 Take-off, -15 Landing
LE = 'slat';                            %<----INPUT: Choose your Leading edge HLD Look in HLD function to check the names.
TE= 'fowler';                       %<----INPUT: Choose your Trailing edge HLD. Look in HLD function to check the names.
c_ac_c = 1.2;                       %<----INPUT: How far is the chord extended in for example a fowler flap. This ratio>1
Delta_Clmax_HLD = HLD(LE,TE,c_ac_c);
%METHOD FOR CALCULATING THE DIFFERENCE IN CLMAX
HLD_Delta_CLmax = 0.9 * Delta_Clmax_HLD * Swf_S * cos(sweep_hinge); %The shift in CLMax due to the HLD, THIS IS NOT THE FINAL VALUE.
%METHOD FOR CALCULATING THE AREA RATIO IF DELTA CLMAX IS KNOWN.
<<<<<<< HEAD
Swf_S= HLD_Delta_CLmax/0.9/Delta_Clmax_HLD/cos(sweep_hinge);
=======
%Swf_S= HLD_Delta_CLmax/0.9/Delta_Clmax_HLD/cos(sweep_hinge);
>>>>>>> 8e04a7f2bbf68aa5db112dc28670559be59188f5

HLD_Delta_alpha0L = delta_alpha0L_airfoil * Swf_S*cos(sweep_hinge); %The shift in alpha 0 Lift, THIS IS NOT THE FINAL VALUE.


TotalCLMAX=HLD_Delta_CLmax+CLmax_clean;
Total_alpha_0L=HLD_Delta_alpha0L+alpha0l ;
%% LIFT CURVE SLOPE CALCULATION DUE TO THE FLAPS OR HLD IN GENERAL.


S_ac_S = 1 + Swf_S * (c_ac_c - 1);          %Increased wing area ratio: Relation from Raymer (slide 59 notes.)
<<<<<<< HEAD
CL_alpha_HLD=S_ac_S*CL*CL_alpha_clean;
=======
CL_alpha_HLD=S_ac_S*CL_alpha_clean;
>>>>>>> 8e04a7f2bbf68aa5db112dc28670559be59188f5

end
