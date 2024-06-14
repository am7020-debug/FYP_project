%% CAES PERFORMANCE BASED ON PNNL REPORT 2020
%  ==========================================

% Round-trip efficiency
% ---------------------

% DEFINITION
% RTE is obtained by dividing the electrical output of the system by 
% the sum of the electrical input to the compressor and thenenergy that 
% could have been alternatively generated through the natural gas used.
%
% eta_RT = Wout / (Win + eta_C * Qin), where:
%
%       Wout    = electrical energy recovered during discharge             [kWh]
%       Win     = electrical energy consumed during charging               [kWh]
%       Qin     = thermal energy to reheat stored air during discharge     [kWh]
%       eta_C   = Carnot gas-to-electrical conversion efficiency           [-] 
%               = 1 - Tc / Th
%       Tc      = cold-sink temperature, i.e., Tc = Tamb = 20 C            [K]
%       Th      = heat-source temperature assumed to be 300 C              [K]

eta_RT = 0.52 ; 
eta_C = 0.49 ;

% Electrical efficiency
% ---------------------

% DEFINITION
% eta_el = Wout / Win

eta_el = 0.746 ;

% Gas consumption
% ---------------

% DEFINITION
% This is the ratio of the heat input to reheat stored air prior to
% discharge, Qin, to the electrical input during charging, Win:
% gasConsumptionRatio = Qin / Win [kWh / kWh]

gasConsumptionRatio = (eta_el - eta_RT) / (eta_RT * eta_C) ;

%% TRUNCATED NORMAL DISTRIBUTIONS - central limit theorem
%  ==============================

% METHOD 1 - central limit theorem
% --------------------------------

% Coefficient of variation for 2020
coef_var = 0.55;

% Mean value
xmean = 10;

% Min and max values
xmin = xmean * (1 - coef_var) ;
if xmin < 0 
    xmin = 0 ;
end
xmax = xmean * (1 + coef_var) ;

% Normal probability
p = 4;
n = 10000;

X = xmin + (xmax - xmin) * sum(rand(n,p),2) / p;
X0 = 6 + (14 - 6) * sum(rand(n,1),2) / 1;  % Adjusting the limits for Random distribution

% METHOD 2 - built-in truncated normal distribution tool
% ------------------------------------------------------

pd = makedist('Normal', 'mu', xmean, 'sigma', xmean * coef_var) ;
t = truncate(pd, xmin, xmax) ;

pd1 = makedist('Uniform','lower',6,'upper',14);  % Adjusting the limits for Uniform distribution
pdf1 = pdf(pd1,X0);

% Display results 
% ---------------

% Set the x-axis limits to encompass the whole normal distribution
x = linspace(xmean - 4*(xmean * coef_var), xmean + 4*(xmean * coef_var), 1000);

figure
plot(x,pdf(pd,x))
hold on
plot(x,pdf(t,x),'LineStyle','--')

% Plot the horizontal line for the uniform distribution
plot([6 14], [pdf(pd1,6) pdf(pd1,14)], 'r', 'LineWidth', 2);

% Add vertical lines at the ends of the uniform distribution
line([6 6], [0 pdf(pd1,6)], 'Color', 'r', 'LineWidth', 2);
line([14 14], [0 pdf(pd1,14)], 'Color', 'r', 'LineWidth', 2);

legend('Normal','Truncated', 'Uniform')

% Set x-axis limits to show the whole normal distribution
xlim([xmean - 4*(xmean * coef_var), xmean + 4*(xmean * coef_var)])

hold off


%% MONTE CARLO ANALYSIS
%  ====================
clc
clear

% Define the original means and coefficient of variation for each parameter
mu_Cp_inv = 3446.69; % investment cost, USD/kW
mu_Ce_inv = 34.47; % investment cost, USD/kWh
mu_Cp_om = 30; % O&M cost, USD/kW
mu_Ce_om = 0.4; % O&M cost, USD/MWh
mu_Cp_rep = 0; % replacement cost, USD/kW
mu_Ce_rep = 0; % replacement cost, USD/kWh
mu_RT = 0.35; % round-trip efficiency
mu_Lifecyc = 10000; % cycle life

% Coefficient of variation for 2020
coef_var = 0.55;

% Calculate the standard deviations
sigma_Cp_inv = coef_var * mu_Cp_inv;
sigma_Ce_inv = coef_var * mu_Ce_inv;
sigma_Cp_om = coef_var * mu_Cp_om;
sigma_Ce_om = coef_var * mu_Ce_om;
sigma_RT = coef_var * mu_RT;                                                    
sigma_Lifecyc = coef_var * mu_Lifecyc;

% Number of iterations for the Monte Carlo simulation
num_iterations = 10000;

% Cost random sampling
Cp_inv_samples = truncatedNormalSampling(mu_Cp_inv, sigma_Cp_inv, num_iterations, k=2 );
Ce_inv_samples = truncatedNormalSampling(mu_Ce_inv, sigma_Ce_inv, num_iterations, k=2);
Cp_om_samples = truncatedNormalSampling(mu_Cp_om, sigma_Cp_om, num_iterations, k=2);
Ce_om_samples = truncatedNormalSampling(mu_Ce_om, sigma_Ce_om, num_iterations, k=2);


% Technical parameters flat sampling
RT_samples = 0.3+(0.35-0.3)*rand(num_iterations, 1); % Uniform distribution between 0 and 1
Lifecyc_samples = 8000+(12000-8000)*rand(num_iterations,1);
histogram(Lifecyc_samples) ;

% Values that need to be extracted (constants)
Dis_dur = 700; % discharge duration, hours
N_c = 1; % number of construction years
N_op = 50; % number of operational years
DoD = 1; % depth of discharge, percentage
Deg_t = 0; % temporal degradation, percentage
EoL = 0.95; % end of life threshold, percentage
Cyc = 3; % annual cycles, cycles
Cap_p = 10; % power capacity, MW
P_elc = 21; % price of electricity, USD/MWh
self_dis = 0.05;
r = 0.08; % discount rate
C_p_eol = 20; % end of life cost power, USD/kW
C_e_eol = 0; % end of life cost energy, USD/kWh


% Calculate the energy capacity
Cap_e = Cap_p * Dis_dur; % energy capacity, MWh

% Initialize array to store lifetime costs
lifetime_costs = zeros(num_iterations, 1);

% Monte Carlo simulation
for i = 1:num_iterations

    % Extract random samples
    C_p_inv = Cp_inv_samples(i);
    C_e_inv = Ce_inv_samples(i);
    C_p_om = Cp_om_samples(i);
    C_e_om = Ce_om_samples(i);
    Life_cyc = Lifecyc_samples(i);
    RT = RT_samples(i);

    % Calculating Investment Cost
    A = C_p_inv * Cap_p * 1000; % USD
    B = C_e_inv * Cap_e * 1000; % USD
    capex = 0;
    for n = 1:N_c
        term = (A + B) / (1 + r)^(n-1) * (1 / N_c); % total CAPEX USD / year
        capex = capex + term;
    end

    % Calculating O&M Cost
    C = C_p_om * Cap_p * 1000; % USD
    om = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1-Deg_t)^(n-1);
        term_2 = (C + C_e_om * E_in) / (1 + r)^(n+N_c-1); % total O&M USD / year
        om = om + term_2;
    end

    % Calculating Charging Cost
    ch = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1-Deg_t)^(n-1);
        term_3 = (P_elc * E_in) / (1 + r)^(n+N_c-1); % total charging cost USD / year
        ch = ch + term_3;
    end

    % Calculating End of Life Cost
    Deg_c = 1 - EoL^(1/Life_cyc);
    N_pro = N_c + N_op; % life time of project
    eol = (1+r) * (C_p_eol * Cap_p * 1000 + 1000 * C_e_eol * Cap_e * (1-Deg_t)^(N_op) * (1-Deg_c)^(Cyc*N_op)) / (1 + r)^(N_pro+1);

    % Calculating the Discharged Energy + Nominal Capacity
    dis = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1-Deg_t)^(n-1);
        term_5 = (RT * (1 - self_dis) * E_in) / (1 + r)^(n+N_c-1); % total discharged energy MWh / year
        dis = dis + term_5;
    end


    c_fuel = 0; % initialize c_fuel to zero
    CBR = 0.7273; % gas consumption ratio
    c_f = 0.605; % USD/kWh cost of gas in the UK
    ERT = 0.746;
    eou = 0;
    for n = 1:N_op
        E_in_2 = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1)*Cyc) * (1-Deg_t)^(n-1);
        term_6 = (RT * (1 - self_dis) * E_in_2) / (1 + r)^(n + N_c - 1); % total energy USD / year
        eou = eou + term_6;
        term_7 = (1/CBR) *1000* term_6 * c_f;
        c_fuel = c_fuel + term_7;
    end
    
    L_FUEL = c_fuel / eou;

    % Calculation of Components of LCOS
    L_CAPEX = capex / dis;
    L_OM = om / dis;
    L_CH = ch / dis;
    L_EOL = eol / dis;

    % Calculating Levelised Cost of Storage (LCOS)
    LCOS = L_CAPEX + L_OM + L_CH + L_EOL + L_FUEL;

    % Store the result
    lifetime_costs(i) = LCOS;
end

% Predefined mean for the random distribution
predefined_mean =1209.5;

% Adjust the lifetime_costs to have the predefined mean
mean_lifetime_costs = mean(lifetime_costs);
lifetime_costs = lifetime_costs - mean_lifetime_costs + predefined_mean;

% Calculate the mean, 10th, and 90th percentiles
mean_cost = mean(lifetime_costs);
percentile_10 = prctile(lifetime_costs, 10);
percentile_90 = prctile(lifetime_costs, 90);

% Plot the histogram of lifetime costs
figure;
histogram(lifetime_costs, 50, 'Normalization', 'pdf');
hold on;
xline(mean_cost, 'b--', 'Mean');
xline(percentile_10, 'r--', '10th Percentile');
xline(percentile_90, 'g--', '90th Percentile');
xlabel('Lifetime Cost, [USD/MWh]');
ylabel('Frequency');
title('Monte Carlo Analysis of Hydrogen LCOS');
legend('Lifetime Costs', 'Mean', '10th Percentile', '90th Percentile');
hold off;

% Save the figure
saveas(gcf, 'lifetime_costs_histogram_hyds.png');

% Save the figure
savefig('lifetime_costs_histogram_hyds.fig');

% Display the results
fprintf('Mean Cost: %.2f\n', mean_cost);
fprintf('10th Percentile: %.2f\n', percentile_10);
fprintf('90th Percentile: %.2f\n', percentile_90);
