clear
clc

% Define the original means and coefficient of variation for each parameter
mu_Cp_inv = 700; % investment cost, USD/kW
mu_Ce_inv = 450; % investment cost, USD/kWh
mu_Cp_om = 10; % O&M cost, USD/kW
mu_Ce_om = 2; % O&M cost, USD/MWh
mu_Cp_rep = 0; % replacement cost, USD/kW
mu_Ce_rep = 0; % replacement cost, USD/kWh
mu_RT = 0.68; % round-trip efficiency
mu_Lifecyc = 20000; % cycle life

% Coefficient of variation for 2020
coef_var = 0.55;

% Calculate the standard deviations
sigma_Cp_inv = coef_var * mu_Cp_inv;
sigma_Ce_inv = coef_var * mu_Ce_inv;
sigma_Cp_om = coef_var * mu_Cp_om;
sigma_Ce_om = coef_var * mu_Ce_om;
sigma_Cp_rep = coef_var * mu_Cp_rep;
sigma_Ce_rep = coef_var * mu_Ce_rep;
sigma_RT = coef_var * mu_RT;
sigma_Lifecyc = coef_var * mu_Lifecyc;

% Number of iterations for the Monte Carlo simulation
num_iterations = 1000;

% Generate random samples for each parameter
Cp_inv_samples = normrnd(mu_Cp_inv, sigma_Cp_inv, num_iterations, 1);
Ce_inv_samples = normrnd(mu_Ce_inv, sigma_Ce_inv, num_iterations, 1);
Cp_om_samples = normrnd(mu_Cp_om, sigma_Cp_om, num_iterations, 1);
Ce_om_samples = normrnd(mu_Ce_om, sigma_Ce_om, num_iterations, 1);
Cp_rep_samples = normrnd(mu_Cp_rep, sigma_Cp_rep, num_iterations, 1);
Ce_rep_samples = normrnd(mu_Ce_rep, sigma_Ce_rep, num_iterations, 1);
RT_samples = normrnd(mu_RT, sigma_RT, num_iterations, 1);
Lifecyc_samples = normrnd(mu_Lifecyc, sigma_Lifecyc, num_iterations, 1);

% Values that need to be extracted (constants)
Dis_dur = 700; % discharge duration, hours
N_c = 1; % number of construction years
N_op = 34; % number of operational years
DoD = 1; % depth of discharge, percentage
Deg_t = 0.0015; % temporal degradation, percentage
EoL = 0.95; % end of life threshold, percentage
Cyc = 3; % annual cycles, cycles
Cap_p = 10; % power capacity, MW
P_elc = 50; % price of electricity, USD/MWh
self_dis = 0;
r = 0.08; % discount rate
C_p_eol = 20; % end of life cost power, USD/kW
C_e_eol = -100; % end of life cost energy, USD/kWh

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
    RT = RT_samples(i);
    Life_cyc = Lifecyc_samples(i);

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

    % Calculation of Components of LCOS
    L_CAPEX = capex / dis;
    L_OM = om / dis;
    L_CH = ch / dis;
    L_EOL = eol / dis;

    % Calculating Levelised Cost of Storage (LCOS)
    LCOS = L_CAPEX + L_OM + L_CH + L_EOL;

    % Store the result
    lifetime_costs(i) = LCOS;
end

% Predefined mean for the random distribution
predefined_mean = 13060.7885;

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
xlabel('Lifetime Cost');
ylabel('Frequency');
title('Monte Carlo Analysis of Flow Battery LCOS');
legend('Lifetime Costs', 'Mean', '10th Percentile', '90th Percentile');
hold off;


% Save the figure
saveas(gcf, 'lifetime_costs_histogram_red.png');

% Save the plot as a .fig file for future editing
savefig('lifetime_costs_histogram_red.fig');

% Display the results
fprintf('Mean Cost: %.2f\n', mean_cost);
fprintf('10th Percentile: %.2f\n', percentile_10);
fprintf('90th Percentile: %.2f\n', percentile_90);
