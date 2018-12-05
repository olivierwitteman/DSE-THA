function [Cldes, CLdes] = Airfoilselection(LAMBDA, rho, Wfiml, Vinf, WSsc, WSec)
Ltot = Wfiml;        %Adsee II - Lecture 1 - Slide 48 - Test for Robel
Lwing = 1.1*Ltot; 
CLdes = 1.1/(0.5*rho*Vinf^2)*0.5*(WSsc+WSec); %Adsee II - Lecture 1 - Slide 49
%Vinf is Vcruise i think
Cldes = CLdes/(cos(LAMBDA)^2);
end

%Adsee II

%Airfoil Selection Steps
% 0 - Assume the t/c is given
% 1 - Determine the airfoil design lift coefficient Cldes
% 2 - For the Cldes, look for the airfoil with minium drag at Cldes and the
% widest possible drag bucket around Cldes
% 3 - Select an airfoil with the largest Clmax possible. However avoid
% airfoils with sharp drop in Cl immediately after stall
% 4 - After having identified airfoils with low enough Cdmin and high
% enough Clmax select the one with the lowest Cm possible at Cldes

%INPUTS
%LAMBDA =            %Wingsweep at 0.25MAC
%rho=                %rho at cruise alt
%Wfiml=              %Aircraft weight at fuel intensive mission leg %ADSEEII-LECTURE1-SLIDE48
%Vinf =              %Cruise Speed
%WSsc =              %Wing loading at the start of the cruise
%Wsec =              %Wing loading at the end of the cruise
