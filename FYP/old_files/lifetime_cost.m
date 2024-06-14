clear
clc

% Define the filename and technology details
filename = 'tech_cost_specification.xlsx';

% List of technologies
technologies = {'Redox flow battery', 'Lithium ion battery', 'Pumped hydro', 'CAES', 'Hydrogen'};

% Define colors for each technology
colors = [0.6, 0.8, 0.6; 1, 0.8, 0.8; 0.8, 0.6, 0.8; 1, 1, 0.8; 0.8, 0.8, 1];

% Initialize the figure for plotting
figure;
hold on;

% Loop through each technology to calculate and plot LCOS
for t = 1:length(technologies)
    sheet = technologies{t};
    data = readtable(filename, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    
    % Extract values from the table using exact column names
    Dis_dur = data.value(strcmp(data.('value'), 'Dis_dur'));
    C_p_inv = data.value(strcmp(data.('value'), 'C_p_inv'));
    C_e_inv = data.value(strcmp(data.('value'), 'C_e_inv'));
    N_c = data.value(strcmp(data.('value'), 'N_c'));
    C_p_om = data.value(strcmp(data.('value'), 'C_p_om'));
    C_e_om = data.value(strcmp(data.('value'), 'C_e_om'));
    N_op = data.value(strcmp(data.('value'), 'N_op'));
    DoD = data.value(strcmp(data.('value'), 'DoD'));
    RT = data.value(strcmp(data.('value'), 'RT'));
    Deg_t = data.value(strcmp(data.('value'), 'Deg_t'));
    EoL = data.value(strcmp(data.('value'), 'EoL'));
    Life_cyc = data.value(strcmp(data.('value'), 'Life_cyc'));
    C_p_eol = data.value(strcmp(data.('value'), 'C_p_eol'));
    C_e_eol = data.value(strcmp(data.('value'), 'C_e_eol'));
    self_dis = data.value(strcmp(data.('value'), 'self_dis'));
    
    % Ensure values are scalars
    Dis_dur = Dis_dur(1);
    C_p_inv = C_p_inv(1);
    C_e_inv = C_e_inv(1);
    N_c = N_c(1);
    C_p_om = C_p_om(1);
    C_e_om = C_e_om(1);
    N_op = N_op(1);
    DoD = DoD(1);
    RT = RT(1);
    Deg_t = Deg_t(1);
    EoL = EoL(1);
    Life_cyc = Life_cyc(1);
    C_p_eol = C_p_eol(1);
    C_e_eol = C_e_eol(1);
    self_dis = self_dis(1);
    
    % Calculation Constants
    Cap_p = 10; % power capacity, MW
    Cap_e = Cap_p * Dis_dur; % energy capacity, MWh
    A = C_p_inv * Cap_p * 1000; % USD
    B = C_e_inv * Cap_e * 1000; % USD
    r = 0.08; % discount rate
    Cyc = 3; % annual cycles, cycles
    P_elc = 50; % price of electricity, USD/MWh

    %% CALCULATING INVESTMENT COST
    capex = 0;
    for n = 1:N_c
        term = (A + B) / (1 + r)^(n-1) * (1 / N_c); % total CAPEX USD / year
        capex = capex + term;
    end

    %% CALCULATING O&M COST
    C = C_p_om * Cap_p * 1000; % USD
    om = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1 - Deg_t)^(n-1);
        term_2 = (C + C_e_om * E_in) / (1 + r)^(n + N_c - 1); % total O&M cost USD / year
        om = om + term_2;
    end

    %% CALCULATING CHARGING COST
    ch = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1 - Deg_t)^(n-1);
        term_3 = (P_elc * E_in) / (1 + r)^(n + N_c - 1); % total charging cost USD / year
        ch = ch + term_3;
    end

    %% CALCULATING END OF LIFE COST 
    Deg_c = 1 - EoL^(1/Life_cyc);
    N_pro = N_c + N_op; % life time of project
    eol = (1 + r) * (C_p_eol * Cap_p * 1000 + 1000 * C_e_eol * Cap_e * (1 - Deg_t)^(N_op) * (1 - Deg_c)^(Cyc * N_op)) / (1 + r)^(N_pro + 1);

    %% CALCULATING THE DISCHARGED ENERGY + NOMINAL CAPACITY
    dis = 0;
    for n = 1:N_op
        E_in = ((Cap_e * DoD * Cyc) / RT) * (EoL^(1/Life_cyc))^((n-1) * Cyc) * (1 - Deg_t)^(n-1);
        term_5 = (RT * (1 - self_dis) * E_in) / (1 + r)^(n + N_c - 1); % total discharged energy MWh / year
        dis = dis + term_5;
    end

    %% CALCULATION OF COMPONENTS OF LCOS
    L_CAPEX = capex / dis;
    L_OM = om / dis;
    L_CH = ch / dis;
    L_EOL = eol / dis;

    %% CALCULATING LEVELISED COST OF STORAGE
    LCOS = L_CAPEX + L_OM + L_CH + L_EOL;

    % Plotting the LCOS breakdown
    values = [L_CAPEX, L_OM, L_CH, L_EOL];
    b = bar(t, values, 'stacked');
    for k = 1:length(b)
        b(k).FaceColor = colors(k, :);
    end

    % Adding value labels to each bar component
    y_offset = 0; % Initial offset for the label
    for i = 1:length(values)
        y_offset = y_offset + values(i); % Increment the offset
        text(t, y_offset - values(i) / 2, sprintf('$%.2f', values(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% Label the y-axis
ylabel('LCOS in $/MWh');

% Set x-axis labels
set(gca, 'XTick', 1:length(technologies), 'XTickLabel', technologies);

% Add title to the plot
title('LCOS Components for Different Technologies');

% Add a legend
legend({'Investment', 'O&M', 'Charging', 'End-of-life'}, 'Location', 'eastoutside');

% Set the x-axis limit to provide space for the labels
xlim([0.5 length(technologies) + 0.5]);

hold off;
