%University of Pittsburgh PropLab
%Injector calculations
% Nathan Smith
clc
clear 
close all


%% Injector Geometry
%Set by us
L_star = 1.5; %[m] characteristic length
d_chamber = 0.08; %[m] chamber diameter
A_star = 2.5447e-04; %[m^2] throat area, from Ben CEA

A_chamber = pi*d_chamber^2/4; %[m^2]
%display(A_chamber)
% A_chamber = 0.0050 [m^2]


%% Calculating Orifice Size
%Equation: [A = wdot/(Cd*sqrt(2*g*rho*dP))]


Cd_ox = 0.625;  %Assume L/d at 1.2 and and R/d at 0
Cd_fuel = 0.625; %Assume L/d at 1.2 and and R/d at 0

mdot_ox = 0.4245;  %[kg/s]  From Ben CEA
mdot_fuel = 0.2022; %[kg/s] From Ben CEA

wdot_ox = mdot_ox*2.20462; %[lb/s]
wdot_fuel = mdot_fuel*2.20462; %[lb/s]

dP = 100*144; %[psf]

nfuel = 8; % number of fuel orifices
nox= 16 ; % number of oxidizer orifices

rho_ox = 76.5367; %[lb/ft^3] density of N2O
rho_fuel = 49.0684; %[lb/ft^3] density of Iso

g = 32.2; %[ft/s^2] - gravitational field strength

Afuel = wdot_fuel/(Cd_fuel*sqrt(2*g*dP*rho_fuel)); %[ft^2]
Afuel = Afuel*144; %[in^2]
Afuel = Afuel/nfuel; % indivuidual orifice area Iso

Dfuel = 2*sqrt(Afuel/pi); %[in]
display(Dfuel);

Aox = wdot_ox/(Cd_ox*sqrt(2*g*dP*rho_ox)); %[ft^2]
Aox = Aox*144; %[in^2]
Aox = Aox/nox; % individual orifice area N2O

Dox = 2*sqrt(Aox/pi); %[in]
display(Dox)

%Dox = 0.0451 [in]
%Dfuel = 0.0492 [in]

%% Impingement Distance & Orifice Seperation

%Using L_impingement/D_avg = 6.8
D_avg = (Dox + Dfuel)/2; %[in]
L_impingement = 6.8*D_avg; %[in]
%display(L_impingement)
%L_impingement = 0.3704 [in]


%Using seperation = L_impingement/tand(theta) --> assuming theta is
%taken from horizontal
theta = 60; %[degrees] 
seperation = L_impingement/tand(theta); %[in]
%display(seperation)
%seperation = 0.2139 [in]




