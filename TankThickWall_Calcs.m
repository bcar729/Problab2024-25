% University of Pittsburgh PropLab
% Tank Volume & Dimensioning Calculations
% Thick Wall Pressure calculations
% Nathan Smith


%% Calculating Stresses accounting for Internal and External Pressures 

%Radial Stress given by: Sigma_r1 = (a^*p_i-b^2*p_o)/(b^2-a^2) - a^2*b^2*(p_i-p_o)/((b^2-a^2)*r)
%Maximum Sigma_r1 occurs at r = a with Sigma_r1 = -800 [psi]
%Maximum Sigma_t1 occurs at r = a and is given by the following

a = 3.5; %[in]
p_i = 800; %[psi]
p_o = 14.7; %[psi]
b = 4;


Sigma_t1 = (a^2*p_i-b^2*p_o)/(b^2-a^2) + (p_i-p_o)*b^2*a^2/((b^2-a^2)*a^2); %[psi]
%display(Sigma_t1) % Sigma_t1 = 5901.2 [psi] 

%The axial stress is uniform throughout the wall and given by

Sigma_a = (a^2*p_i-b^2*p_o)/(b^2-a^2); %[psi]
%display(Sigma_a); %Sigma_a = 2550.6 [psi]


%% Calculating Stresses accounting for Internal Pressures Only

%By applying the Maximum Shear Stress Theory of Failure
%F_smax = (F_tmax-F_rmax)/2

F_smax = b^2*p_i/(b^2-a^2);
%display(F_smax) %F_smax = 3413.3 [psi]


%% Tank Volume & Height Calculations

% Tank volume
m_dot_fuel = 0.2002; %[kg/s] From Ben CEA
m_dot_ox = 0.4245; %[kg/s] From Ben CEA
burn_time = 4; %[s]
rho_fuel = 786/61023.7; %[kg/in^3]
rho_ox = 1220/61023.7; %[kg/in^3]

Ox_volume = m_dot_ox*burn_time/rho_ox; %[in^3]
%display(Ox_volume) %Ox_volume = 84.9330 [in^3]
Fuel_volume = m_dot_fuel*burn_time/rho_fuel; %[in^3]
%display(Fuel_volume) %Fuel_volume = 62.1727 %[in^3]

Total_volume = Ox_volume + Fuel_volume; %[in^3]
%display(Total_volume)  %Total_volume = 141.1057 [in^3]


%Section Heights
Ox_Height = Ox_volume/pi/a^2; % [in]
%display(Ox_Height) %Ox_Height = 2.2069 [in]

Fuel_Height = Fuel_volume/pi/a^2; %[in]
%display(Fuel_Height) %Fuel_Height = 1.6115 [in]



