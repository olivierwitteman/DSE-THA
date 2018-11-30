%input values for Cl-alpha curve for finite wing (clean and flapped
%General inputs
A=0.002;
sweepc2 =0.25; 
M=0.3;
Beta=sqrt(1-M^2);
eta=0.95;
%%Airfoil ans wing
sweepLE=0.25; %Calculate
alpha0L=0; %Given by the airfoil
CLdes=0.65; %Given
taper=0.6; %Given
Sharpfactor=0.02; %from the airfoil.
%%statistical
%general
C1=0.02; %From Table (slide:15)
C2=0.02  %From table slide: 15 raymer.


CLmax_base=0.2;
%% Use the datcom method to get the CL clope from statistics
Datcomtop=2*pi*A; 
Datcombottom=sqrt(4+(A*Beta/eta)^2*(1+((tan(sweepc2))^2/Beta^2)))+2;
CLalpha=Datcomtop/Datcombottom;  

%next is the trim angle.this is the angle the aircraft needs to fly at to
%fly at CLdes

alphatrim=CLdes/CLalpha+alpha0L;

%Next is the CLmax, 2 methods are being introduced. The Datcom method and a
%general  one. Datcom is prefferred. CLmax
Datcom_choose_method=((C1+1)*cos(sweepLE));
if A> 4/Datcom_choose_method
    delta_CLmax=0.5; %Term for M>0.2. Not for take-off and landing.
    delta_alpha_CLmax=0.02; %From table slide 19 (Raymer)
    CLM_clmax=0.05; %slide 17 use right table 
    %high aspect ratio method
    CLmax=CLM_clmax*clmax+delta_CLmax;
    alpha_stall=CLmax/CLalpha+alpha0L+delta_alpha_CLmax;
%otherwise low aspect ratio method.    
elseif A< 4/Datcom_choose_method
    CLmax_base=0.2; %table from Raymer
    delta_CLmax=0.5;
    delta_alpha_CLmax=0.02;
    alpha_CLmax=0.04;
    
    CLmax=CLmax_base+delta_CLmax;
    alpha_stall=alpha_CLmax+delta_alpha_CLmax;
end
%Slide 21 verdergaan met alle dingen beschrijven en schrijven waar ze
%vandaan komen. Verder zijn het rare namen voor de variables die je
%misschien aan zou kunnen passen. Dan kun je verder gaan op slide 24 met
%HLD.
