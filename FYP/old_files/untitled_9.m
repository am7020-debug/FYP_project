clear
clc

% Load the Excel file
filename = 'tech_cost_specification.xlsx';

% List of technologies
technologies = {'Redox flow battery', 'Lithium ion battery', 'Pumped hydro', 'CAES', 'Hydrogen'};

% Variables to be extracted
required_vars = {'Dis_dur', 'C_p_inv', 'C_e_inv', 'N_c', 'C_p_om', 'C_e_om', 'N_op', 'DoD', 'RT', 'Deg_t', 'EoL', 'Life_cyc', 'C_p_eol', 'C_e_eol', 'self_dis', 'cof_heat'};

% Preallocate storage for LCOS components and results
L_CAPEX = zeros(length(technologies), 1);
L_OM = zeros(length(technologies), 1);
L_CH = zeros(length(technologies), 1);
L_EOL = zeros(length(technologies), 1);
L_FUEL = zeros(length(technologies), 1);  % New cost component for fuel
LCOS_values = zeros(length(technologies), 1);

% Loop through each technology
for i = 1:length(technologies)
    tech = technologies{i};
    
    % Read data from the corresponding sheet
    data = readtable(filename, 'Sheet', tech, 'VariableNamingRule', 'preserve');
    
    % Normalize column names
    var_name_column = 'Variable name in code';
    
    % Extract values from the data table
    values = struct();
    for j = 1:length(required_vars)
        var_name = required_vars{j};
        idx = strcmp(data.(var_name_column), var_name);
        if any(idx)
            values.(var_name) = data{idx, 'value'};
        else
            values.(var_name) = NaN;
        end
    end
    
    % Assign values to variables
    Dis_dur = values.Dis_dur;
    C_p_inv = values.C_p_inv;
    C_e_inv = values.C_e_inv;
    N_c = values.N_c;
    C_p_om = values.C_p_om;
    C_e_om = values.C_e_om;
    N_op = values.N_op;
    DoD = values.DoD;
    RT = values.RT;
    Deg_t = values.Deg_t;
    EoL = values.EoL;
    Life_cyc = values.Life_cyc;
    C_p_eol = values.C_p_eol;
    C_e_eol = values.C_e_eol;
    self_dis = values.self_dis;
    cof_heat = values.cof_heat;  % Cost of fuel per unit energy
    Cap_p = 10; % power capacity, MW
  
    % Fixed values
    r = 0.08; % discount rate
    Cyc = 3; % annual cycles, cycles
    P_elc = 50; % price of electricity, USD/MWh
    
    % Calculations
    Cap_e = Cap_p * Dis_dur; % energy capacity, MWh
    A = C_p_inv * Cap_p * 1000; % USD
    B = C_e_inv * Cap_e * 1000; % USD

    % Calculate CAPEX
    capex = 0;
    for n = 1:N_c
        term = (A + B) / (1 + r)^(n-1) * (1 / N_c); % total CAPEX USD / year
        capex = capex + term;
    end
    display (capex);

    % Calculate O&M cost
    C = C_p_om * Cap_p * 1000; % USD
    om = 0;
    E_out = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1)*Cyc) * (1-Deg_t)^(n-1);
        term_6 = (RT * (1 - self_dis) * E_in) / (1 + r)^(n + N_c - 1); % total energy USD / year
        E_out = E_out + term_6;
        
        term_2 = (C + C_e_om * E_in + E_out *cof_heat) / (1 + r)^(n + N_c - 1); % total O&M USD / year
        om = om + term_2;
    end

    % Calculate charging cost
    ch = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1)*Cyc) * (1-Deg_t)^(n-1);
        term_3 = (P_elc * E_in) / (1 + r)^(n + N_c - 1); % total charging USD / year
        ch = ch + term_3;
    end

    % Calculate end of life cost
    Deg_c = 1 - EoL^(1/Life_cyc);
    N_pro = N_c + N_op; % lifetime of project
    eol = (1+r) * (C_p_eol * Cap_p * 1000 + 1000 * C_e_eol * Cap_e * (1-Deg_t)^N_op * (1-Deg_c)^(Cyc*N_op)) / (1 + r)^(N_pro + 1);

    % Calculate discharged energy
    dis = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1)*Cyc) * (1-Deg_t)^(n-1);
        term_5 = (RT * (1 - self_dis) * E_in) / (1 + r)^(n + N_c - 1); % total energy USD / year
        dis = dis + term_5;
    end



    % Calculate LCOS components
    L_CAPEX(i) = capex / dis;
    L_OM(i) = om / dis;
    L_CH(i) = ch / dis;
    L_EOL(i) = eol / dis;


    % Final LCOS calculation
    LCOS_values(i) = L_CAPEX(i) + L_OM(i) + L_CH(i) + L_EOL(i) ;
    
    % Display results with comma formatting
    fprintf('Technology: %s\n', tech);
    formatted_LCOS = sprintf('%.2f', LCOS_values(i));
    fprintf('Levelized cost of storage (LCOS): $%s / MWh\n', insert_commas(formatted_LCOS));
    fprintf('               ~~~\n');
end

% Create a table with results
results_table = table(technologies', L_CAPEX, L_OM, L_CH, L_EOL, LCOS_values, ...
    'VariableNames', {'Technology', 'CAPEX', 'O&M', 'Charging', 'End_of_Life', 'Fuel', 'LCOS'});

% Save the table to an Excel file
writetable(results_table, 'LCOS_results_with_fuel.xlsx');

% Plotting the stacked bar chart for LCOS breakdown
figure;
bar_data = [L_CAPEX, L_OM, L_CH, L_EOL];
bar(bar_data, 'stacked');
set(gca, 'XTickLabel', technologies);
legend({'CAPEX', 'O&M', 'Charging', 'End of Life', 'Fuel'}, 'Location', 'bestoutside');
ylabel('LCOS ($/MWh)');
title('Breakdown of Levelized Cost of Storage (LCOS) by Technology');
grid on;

function str_with_commas = insert_commas(str)
    % Helper function to insert commas into a number string
    [int_part, dec_part] = strtok(str, '.');
    int_with_commas = regexprep(int_part, '\B(?=(\d{3})+(?!\d))', ',');
    str_with_commas = [int_with_commas, dec_part];
end
