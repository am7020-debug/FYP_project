clear
clc

% Load the Excel file
filename = 'tech_cost_specification.xlsx';

% List of sheet names
sheet_names = {'Redox flow battery', 'Lithium ion battery', 'Pumped hydro', 'CAES', 'Hydrogen'};

% Variables to be extracted
variables_to_extract = {'Cp_inv', 'Ce_inv', 'Cp_om', 'Ce_om', 'Cp_rep', 'Ce_rep', 'RT', 'Lifecyc'};

% Preallocate storage for results
lifetime_costs_all = zeros(1000, length(sheet_names));

% Loop through each sheet
for s = 1:length(sheet_names)
    sheet_name = sheet_names{s};
    
    % Read data from the corresponding sheet
    data = readtable(filename, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
    
    % Normalize column names
    var_name_column = 'Variable name in code';
    
    % Extract values from the data table
    values = struct();
    for j = 1:length(variables_to_extract)
        var_name = variables_to_extract{j};
        idx = strcmp(data.(var_name_column), var_name);
        if any(idx)
            values.(var_name) = data{idx, 'value'};
        else
            values.(var_name) = NaN;
        end
    end
    
    % Assign values to variables
    mu_Cp_inv = values.Cp_inv;
    mu_Ce_inv = values.Ce_inv;
    mu_Cp_om = values.Cp_om;
    mu_Ce_om = values.Ce_om;
    mu_Cp_rep = values.Cp_rep;
    mu_Ce_rep = values.Ce_rep;
    mu_RT = values.RT;
    mu_Lifecyc = values.Lifecyc;

    % Coefficient of variation for 2020
    coef_var = 0.25;

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
    constant_vars = {'Dis_dur', 'N_c', 'N_op', 'DoD', 'Deg_t', 'EoL', 'Cyc', 'Cap_p', 'P_elc', 'self_dis', 'r', 'C_p_eol', 'C_e_eol'};
    constant_values = struct();
    for j = 1:length(constant_vars)
        var_name = constant_vars{j};
        idx = strcmp(data.(var_name_column), var_name);
        if any(idx)
            constant_values.(var_name) = data{idx, 'value'};
        else
            constant_values.(var_name) = NaN;
        end
    end
    
    % Assign constant values to variables
    Dis_dur = constant_values.Dis_dur;
    N_c = constant_values.N_c;
    N_op = constant_values.N_op;
    DoD = constant_values.DoD;
    Deg_t = constant_values.Deg_t;
    EoL = constant_values.EoL;
    Cyc = constant_values.Cyc;
    Cap_p = constant_values.Cap_p;
    P_elc = constant_values.P_elc;
    self_dis = constant_values.self_dis;
    r = constant_values.r;
    C_p_eol = constant_values.C_p_eol;
    C_e_eol = constant_values.C_e_eol;

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
    
    % Store the lifetime costs for this technology
    lifetime_costs_all(:, s) = lifetime_costs;
end

% Calculate the 10th and 90th percentiles for each technology
percentile_10 = prctile(lifetime_costs_all, 10);
percentile_90 = prctile(lifetime_costs_all, 90);

% Plot the histogram of lifetime costs for the first technology as an example
figure;
histogram(lifetime_costs_all(:, 1), 50, 'Normalization', 'pdf');
hold on;
xline(percentile_10(1), 'r--', '10th Percentile');
xline(percentile_90(1), 'g--', '90th Percentile');
xlabel('Lifetime Cost');
ylabel('Frequency');
title(['Monte Carlo Simulation of Lifetime Costs for ', sheet_names{1}]);
legend('Lifetime Costs', '10th Percentile', '90th Percentile');
hold off;

% Calculate the probability of the first technology exhibiting the lowest lifetime cost compared to the second one as an example
lifetime_costs_A = lifetime_costs_all(:, 1);  % First technology
lifetime_costs_B = lifetime_costs_all(:, 2);  % Second technology

prob_A_min = mean(lifetime_costs_A < lifetime_costs_B);

% Display the results for the first technology as an example
fprintf('Technology: %s\n', sheet_names{1});
fprintf('10th Percentile: %.2f\n', percentile_10(1));
fprintf('90th Percentile: %.2f\n', percentile_90(1));
fprintf('Probability of %s having the lowest lifetime cost compared to %s: %.2f\n', sheet_names{1}, sheet_names{2}, prob_A_min);
