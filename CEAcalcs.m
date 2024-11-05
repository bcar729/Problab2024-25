clear
close all
clc

%% Defined by us
m_dot_ratio = 2.2; %5.7 Stoich; %[-]
d_star = 0.018; %[m]
%m_dot = 2.5; %[kg/s] % note: want m_dot_ox to be about 0.4
P_chamber = 3.447e6; %[Pa]
L_star = 1.5; %[m] characteristic length
d_chamber = 0.08; %[m]
f_film = .2; %[-] fraction of total mass flow that goes towards film cooling
%geometry = flipud(load("nozzle points 230412.txt"))/1000;

%Constants
R = 8.314; %[J/mol-K]
P_0 = 101.3e3; %[Pa] ambient
T_amb = 298.15; %[K] ambient temperature                                       

%From CEA
thermo_lib = open("thermo_lib.mat");
CEA_RESULTS = CEA('prob','rocket','equilibrium', ...    %problem type
    'o/f',m_dot_ratio,'p,bar',(P_chamber)./(1e5), ...   % o/f and chamber pressures
    'pi/p', P_chamber./P_0, ...                         %chamber/exit pressure ratios
    'reac','fu','C3H8O,1propanol','wt%',100.,'t(k)',298.15, ...    %RP-1 data
    'ox','N2O','wt%',100,'t(k)',298.15, ...                %N2O data
    'output','transport','end');

gamma_chamber = CEA_RESULTS.output.eql.gamma(1); %[-]
gamma_star = CEA_RESULTS.output.eql.gamma(2); %[-]
gamma_e = CEA_RESULTS.output.eql.gamma(3); %[-]
T_chamber = CEA_RESULTS.output.eql.temperature(1); %[K]
MW_e = CEA_RESULTS.output.eql.mw(3)/1000; %[kg/mol]
viscdyn_chamber = CEA_RESULTS.output.eql.viscosity(1)*1e-4; %[Pa-s] dynamic viscosity
Cp_chamber = CEA_RESULTS.output.eql.cp(1)*1000; %[J/kg-k]
Pr_chamber = CEA_RESULTS.output.eql.prandtl.eql(1); %[-]
C_star_chamber = CEA_RESULTS.output.eql.cstar(1); %[m/s] characteristic velocity


%% Thrust Equations
A_star = (pi/4)*(d_star)^2; %[m^2]
P_e = P_0; %optimize for sea level
M_e = sqrt(2*((P_e/P_chamber)^(-(gamma_e-1)/gamma_e)-1)/(gamma_e-1));
A_e = ((gamma_e+1)/2)^(-(gamma_e+1)/(2*(gamma_e-1)))*(1+(gamma_e-1)/2*M_e^2)^((gamma_e+1)/(2*(gamma_e-1)))*(1/M_e)*A_star; %[m^2] exit area
T_e = (1+(gamma_e-1)/2*M_e^2)^(-1)*T_chamber; %[K]
V_e = M_e*sqrt(gamma_e*R*T_e/MW_e); %[m/s] might need to divide by average molar mass
m_dot = (A_star*P_chamber)/sqrt(T_chamber)*sqrt(gamma_e/(R/MW_e))*((gamma_e+1)/2)^(-(gamma_e+1)/(2*(gamma_e-1)));
%P_chamber = 1/((A_star/m_dot)/sqrt(T_chamber)*sqrt(gamma_e/(R/MW_e))*((gamma_e+1)/2)^(-(gamma_e+1)/(2*(gamma_e-1))));

thrust = (m_dot)*(V_e) + (P_e-P_0)*A_e; %[N]
Isp = thrust/(m_dot*9.81);
thrust_coeff = (thrust)/(A_star*P_chamber);

%Mass flows
m_dot_fuel = m_dot/(m_dot_ratio+1);
m_dot_ox = m_dot_fuel * m_dot_ratio;

%Chamber geometry
A_ratio = A_e/A_star;
d_e = d_star*sqrt(A_ratio);

A_chamber = .25*d_chamber^2*pi;
V_chamber = L_star*A_star;
L_chamber = V_chamber/A_chamber;

%Throat Conditions
T_star = T_chamber*(1+(gamma_e-1)/2)^(-1); %note M = 1 at throat

%% Thermal Performance

%geometry
% A_inner = 0.056801; %[m^2] area of inner surface of whole thrust chamber and nozzle
% n_ch = 40; %[-] number of regen channels
% t_w_in = .001; %[m] thickness of inner wall before channel
% height_ch = .008; %[m] height of cooling channel
% w_ch = .001; %[m] width of cooling channel
% t_w_out = .001; %[m] thickness of outer wall after channel
% 
% t_w_tot = t_w_in + height_ch + t_w_out; %[m] total thickness of thrust chamber;
% w_w_chamber = (((d_chamber+2*t_w_in)*pi)/n_ch)-w_ch; %[m] width of ribs between channels in thrust chamber
% w_w_star = (((d_star+2*t_w_in)*pi)/n_ch)-w_ch; %[m] width of ribs between channels at throat
% w_w_e = (((d_e+2*t_w_in)*pi)/n_ch)-w_ch; %[m] width of ribs between channels at exit
% w_ribs = [w_w_chamber w_w_star w_w_e];
% 
% ratio_rib = w_ribs./(w_ribs + w_ch);
% ratio_ch = [1 1 1];%current model considers only channel w_ch./(w_ribs + w_ch);
% 
% %heat transfer calcs
% m_dot_fuel_reg = m_dot_fuel/(1-f_film); %adjust mass flow used in combustion calcs to include film cooling mass
% 
% %regenerative coolant properties
% 
% k_inco = 11; %[W/m-K] thermal conductivity of inconel 718
% 
% % Re_g = u_g.*rho_g.*D_h_g./viscdyn_g; %[-] reynolds number gases
% % Re_reg = u_reg.*rho_reg.*D_h_reg./viscdyn_reg; %[-] reynolds number coolant
% 
% %discrete model
% n_regions = 4; %div, throat, conv, chamber
% n_cells = [20, 20, 20, 20]; %number of cells per region
% S = [.0475 .0156 .0730 .0562]; %[m] arc length of each segment
% dS = [];
% for i = 1:n_regions
%     dS = horzcat(dS,repmat(S(i)/n_cells(i),1,n_cells(i))); %size of each element
% end
% 
% geom = geometry(1:end-2,:);
% 
% A_ratio = ((geometry(:,2)).^2)/(d_star/2)^2;
% [~, i_star] = min(geometry(:,2));
% 
% % for i = 1:length(geom)
% %      if i > i_star %index of the throat
% %          flowType = "sub";
% %      end
% % %      gamma_g(i) = gamma_g_func(geom(i,1),gamma_e,gamma_star,gamma_chamber);
% % %      [M_g(i), T_g(i), P_g(i), rho_g(i)] = flowisentropic(gamma_g(i),A_ratio(i),flowType);
% % %      T_g(i) = T_g(i)*T_chamber;
% % %      P_g(i) = P_g(i)*T_chamber;
% % %      rho_g(i) = rho_g(i)*T_chamber;
% % end
% 
% %initial conditions for first element
% RP1 = py.rocketprops.rocket_prop.get_prop('RP1');
% T_co = zeros(1,length(dS)) + 300; %all coolant initially at 300 [K]
% D_h_co = repmat((4*w_ch*height_ch)/(2*(w_ch + height_ch)),1,length(dS));%[m] regen hydraulic diameter as function of axial position
% rho_co = repmat(RP1.SGLiqAtTdegR(KtoR(T_co(1))).*1000,1,length(dS)); %[kg/m^3] regen density as function of axial position
% viscdyn_co = repmat(RP1.ViscAtTdegR(KtoR(T_co(1))).*0.1,1,length(dS)); %[Pa-s] regn dynamic viscosity as function of axial position
% k_co = repmat(CONDtoMETRIC(RP1.CondAtTdegR(KtoR(T_co(1)))),1,length(dS)); %[W/m-K] regen thermal conductivity
% Cp_co = repmat(CptoMETRIC(RP1.CpAtTdegR(KtoR(T_co(1)))),1,length(dS));
% Pr_co = (Cp_co.*viscdyn_co)./(k_co);%[-] regen prandtl number as function of axial position
% u_co = (m_dot_fuel_reg/n_ch)./(rho_co.*(w_ch * height_ch));  %[m/s] regen speed as function of axial position
% Re_co = u_co.*rho_co.*D_h_co./viscdyn_co; %[-]
% 
% viscdyn_g = zeros(1,size(geometry,1));
% cp_g = zeros(1,size(geometry,1));
% k_g = zeros(1,size(geometry,1));
% Pr_g = zeros(1,size(geometry,1));
% u_g = zeros(1,size(geometry,1));
% gamma_g = zeros(1,size(geometry,1));
% M_g = zeros(1,size(geometry,1));
% T_g = zeros(1,size(geometry,1));
% P_g = zeros(1,size(geometry,1));
% rho_g = zeros(1,size(geometry,1));

% %run CEA for every area ratio to get transport properties
% flowType = 'sup';
% for i = 1:size(geometry,1)
%     if i > i_star %index of the throat, change mode here to subsonic flow
%         flowType = 'sub';
%     end
%     CEA_TRANS = CEA('prob','rocket','equilibrium', ...  %problem type
%         'o/f',m_dot_ratio,'p,bar',(P_chamber)./(1e5), ...   % o/f and chamber pressures
%         flowType, A_ratio(i), ...                         %chamber/exit pressure ratios
%         'reac','fu','RP-1','wt%',100.,'t(k)',298.15, ...    %RP-1 data
%         'ox','N2O','wt%',100,'t(k)',300, ...                %N2O data
%         'output','transport','end');
% 
%     viscdyn_g(i) = CEA_TRANS.output.eql.viscosity(3)*1e-4;
%     cp_g(i) = CEA_TRANS.output.eql.cp_tran.eql(3)*1e3;
%     k_g(i) = CEA_TRANS.output.eql.conduct.eql(3)*1e-1;
%     Pr_g(i) = CEA_TRANS.output.eql.prandtl.eql(3);
%     u_g(i) = CEA_TRANS.output.eql.mach(3)*CEA_TRANS.output.eql.sonvel(3);
%     gamma_g(i) = CEA_TRANS.output.eql.gamma(3);
%     M_g(i) = CEA_TRANS.output.eql.mach(3);
%     T_g(i) = CEA_TRANS.output.eql.temperature(3);
%     P_g(i) = CEA_TRANS.output.eql.pressure(3)*1e5;
%     rho_g(i) = CEA_TRANS.output.eql.density(3);  
% end
% D_h_g = 2*squeeze(geometry(:,2))';
% 
% Re_g = (rho_g.*u_g.*D_h_g)./viscdyn_g;
% 
% figure(1)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Mach Number vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), M_g)
% ylabel("Mach Number [-]")
% hold off
% 
% figure(2)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Gas Temperature vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), T_g)
% ylabel("Temperature [K]")
% hold off
% 
% figure(3)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Pressure vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), P_g)
% ylabel("Pressure [Pa]")
% hold off
% 
% figure(4)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Gas Velocity vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), u_g)
% ylabel("Velocity [m/s]")
% hold off
% 
% figure(5)
% subplot(2,2,1)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Dynamic Viscosity vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), viscdyn_g)
% ylabel("Dynamic Viscosity [Pa-s]")
% hold off
% 
% subplot(2,2,2)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Constant Pressure Heat Capacity vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), cp_g)
% ylabel("Heat Capcity [J/kg-K]")
% hold off
% 
% subplot(2,2,3)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Thermal Conductivity vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), k_g)
% ylabel("Thermal Conductivity [W/m-K]")
% hold off
% 
% subplot(2,2,4)
% hold on
% plot(geometry(:,1),geometry(:,2), 'Color','k','LineWidth',2)
% ylim([0, .1])
% title("Prandtl Number vs. Axial Position")
% ylabel("Radius [m]")
% xlabel("x [m]")
% yyaxis right
% plot(geometry(:,1), Pr_g)
% ylabel("Pradtl Number [-]")
% hold off


% q = zeros(1,length(dS));
% h_g = zeros(1,length(dS));
% q_guess = 1e5;
% q_calc = q_guess + 1e2;
% T_co_guess = T_co(1);
% R_w_in = t_w_in./(k_inco.*dS.*w_ch); %dS * w_ch = area of channel segment
% perimeter_reg = height_ch*2 + w_ch*2;
 
% for i = 1:length(dS)
%     iterations = 0;
% 
%     %resistances at each segment, REYNOLDS NUMBER IS TOO HIGH FOR THIS
%     %CORRELATOIN, SWITCH TO BARTZ OR SOME OTHER METHOD
%     f_gnielinski_g = (0.79*log(Re_g(i))-1.64)^-2; %factor f for gnielinski correlation
%     Nu_gnielinski_g = ((f_gnielinski_g/8)*(Re_g(i)-1000)*Pr_g(i))/(1+12.7*sqrt(f_gnielinski_g/8)*(Pr_g(i)^(2/3)-1));
%     R_g_iter = (Nu_gnielinski_g * k_g(i)/D_h_g(i))^-1;
%     
%     %coolant properties
%     while abs(q_guess - q_calc) > 1e-5
%         q_guess = abs(q_guess + q_calc)/2;
%         
%         T_reg_w_i = T_g(i) - q_guess*(R_g_iter + R_w_in(i)); %update inner wall temp of regen channel
% 
%         %update coolant properties
%         rho_co_iter = RP1.SGLiqAtTdegR(KtoR(T_co_guess)).*1000; %[kg/m^3] regen density as function of axial position
%         viscdyn_co_iter = RP1.ViscAtTdegR(KtoR(T_co_guess)).*0.1; %[Pa-s] regn dynamic viscosity as function of axial position
%         k_co_iter = CONDtoMETRIC(RP1.CondAtTdegR(KtoR(T_co_guess))); %[W/m-K] regen thermal conductivity
%         Cp_co_iter = CptoMETRIC(RP1.CpAtTdegR(KtoR(T_co_guess)));
%         Pr_co_iter = (Cp_co_iter.*viscdyn_co_iter)./(k_co_iter);%[-] regen prandtl number as function of axial position
%         u_co_iter = (m_dot_fuel_reg/n_ch)./(rho_co_iter.*(w_ch * height_ch));  %[m/s] regen speed as function of axial position
%         Re_co_iter = u_co_iter*rho_co_iter.*D_h_co(i)./viscdyn_co_iter; %[-]
% 
%         %calculate coolant thermal resistance
%         R_reg_in_iter =  ((k_co_iter./D_h_co(i)).*.023*(Re_co_iter).^(.8).*Pr_co_iter.^(.4)).^(-1);
%         
%         R_tot = R_g_iter + R_w_in(i) + R_reg_in_iter;
% 
%         %calculate coolant temperature
%         T_co_guess = fzero(@(T_co_iter) (T_reg_w_i-T_co_iter)*exp((-dS(i)*perimeter_reg/2)/((m_dot_fuel_reg/n_ch)*Cp_co_iter*(R_tot)))-(T_reg_w_i-T_co_iter),T_co_guess);
%         q_calc = (T_g(i) - T_co_guess)/(R_tot);
% 
%         iterations = iterations + 1;
%         if mod(iterations,1000) == 0
%             fprintf("q_guess = %.2f \t q_calc = %.2f \t T_co_guess = %.2f \n",q_guess,q_calc,T_co_guess)
%         end
%     end
% 
%     q(i) = q_calc;
%     h_g(i) = (R_g_iter)^-1;
%     T_co(i) = T_co_guess;
%     rho_co(i) = rho_co_iter; %[kg/m^3] regen density as function of axial position
%     viscdyn_co(i) = viscdyn_co_iter; %[Pa-s] regn dynamic viscosity as function of axial position
%     k_co(i) = k_co_iter; %[W/m-K] regen thermal conductivity
%     Cp_co(i) = Cp_co_iter;
%     Pr_co(i) = Pr_co_iter;%[-] regen prandtl number as function of axial position
%     u_co(i) = u_co_iter;  %[m/s] regen speed as function of axial position
%     Re_co(i) = Re_co_iter; %[-]
% 
%     q_guess = q(i);
%     q_calc = q(i) + 100;
%     T_co_guess = T_co(i);
%     fprintf('T_co = %.2f [K]\t q_calc = %.2f [W/m^2] \t iter = %d\n',T_co(i),q(i), iterations)
%     
% 
% end
% 
% plot(geometry(:,1),geometry(:,2))
% ylim([0 55])
% 
% fprintf("mdot ox = %.2f [kg/s] \nmdot fuel = %.2f [kg/s] \n", m_dot_ox, m_dot_fuel_reg)
% fprintf("u_fuel = %.2f [m/s] \n", u_co(1))
% fprintf("chamber pressure = %.4e [Pa] \n", P_chamber)
%% CONVERSION FUNCTIONS
% function gamma_g = gamma_g_func(x,gamma_e,gamma_star,gamma_chamber)
%     if x > 132.45
%         gamma_g = interp1([180 132.45],[gamma_e gamma_star],x);
%     else
%         gamma_g = interp1([132.45 0],[gamma_star gamma_chamber],x);
%     end
% end


