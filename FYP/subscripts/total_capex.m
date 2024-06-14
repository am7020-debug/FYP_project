clear
clc

%% VALUES THAT NEED TO BE EXTRACTED

Dis_dur = 700; % discharge duration, hours
C_p_inv = 250; % investment cost, USD/kW
C_e_inv = 300; % investment cost, USD/kWh
N_c = 1; % number of construction years
C_p_om = 5; % O&M cost, USD/kW
C_e_om = 0.4; % O&M cost, USD/MWh
N_op = 21.8; % number of operational years
DoD = 0.8; % depth of discharge, percentage
RT = 0.86 ; % round trip efficiency, percentage
Deg_t = 0.01; % temporal degradation, percentage
EoL = 0.8; % end of life threshold, percentage
Life_cyc = 3500; %life cycles
C_p_eol = 20; % end of life cost power, USD/kW
C_e_eol = 0; % end of life cost energy, USD/kWh
self_dis = 0.01;

%% CALCULATING INVESTMENT COST

r = 0.08; % discount rate
Cyc = 3; % annual cycles, cycles
Cap_p = 10; % power capacity, MW
P_elc = 50; % price of electricity, USD/MWh
Cap_e = Cap_p * Dis_dur; % energy capacity, MWh
A = C_p_inv * Cap_p*1000; % USD
B = C_e_inv * Cap_e*1000; % USD

% Initialize the sum
capex = 0;

% Loop through the terms of the sum
for n = 1:N_c
    term = (A + B) / (1 + r)^(n-1) * (1 / N_c); % total CAPEX USD / year
    capex = capex + term;
end

% To display with commas, you can use the following approach
formatted_cost = num2str(capex, '%.2f');
formatted_cost_with_commas = regexprep(formatted_cost, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Total CAPEX: $%s\n', formatted_cost_with_commas);

%% CALCULATING O&M COST

C = C_p_om * Cap_p*1000; % USD

% Initialize the sum
om = 0;

% Loop through the terms of the sum
for n = 1:N_op
    E_in = ((Cap_e * DoD * Cyc)/RT)*(EoL^(1/Life_cyc))^((n-1)*Cyc)*(1-Deg_t)^(n-1);
    term_2 = (C + C_e_om * E_in) / (1 + r)^(n+N_c-1); % total CAPEX USD / year
    om = om + term_2;
end

% To display with commas, you can use the following approach
formatted_cost_2 = num2str(om, '%.2f');
formatted_cost_with_commas_2 = regexprep(formatted_cost_2, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Total O&M cost: $%s\n', formatted_cost_with_commas_2);

%% CALCULATING CHARGING COST

% Initialize the sum
ch = 0;

% Loop through the terms of the sum
for n = 1:N_op
    E_in = ((Cap_e * DoD * Cyc)/RT)*(EoL^(1/Life_cyc))^((n-1)*Cyc)*(1-Deg_t)^(n-1);
    term_3 = (P_elc * E_in) / (1 + r)^(n+N_c-1); % total CAPEX USD / year
    ch = ch + term_3;
end

% To display with commas, you can use the following approach
formatted_cost_3 = num2str(ch, '%.2f');
formatted_cost_with_commas_3 = regexprep(formatted_cost_3, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Total charging cost: $%s\n', formatted_cost_with_commas_3);

%% CALCULATING END OF LIFE COST 

Deg_c = 1 - EoL^(1/Life_cyc);
N_pro = N_c + N_op; % life time of project

% initialise
eol = (1+r) * (C_p_eol * Cap_p * 1000 + 1000 * C_e_eol * Cap_e * (1-Deg_t)^(N_op) * (1-Deg_c)^(Cyc*N_op)) / (1 + r)^(N_pro+1);

% To display with commas, you can use the following approach
formatted_cost_4 = num2str(eol, '%.2f');
formatted_cost_with_commas_4 = regexprep(formatted_cost_4, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Total end-of-life cost: $%s\n', formatted_cost_with_commas_4);

%% CALCULATING THE DISCHARGED ENERGY + NOMINAL CAPACITY

% Initialize the sum
dis = 0;

% Loop through the terms of the sum
for n = 1:N_op
    E_in = ((Cap_e * DoD * Cyc)/RT)*(EoL^(1/Life_cyc))^((n-1)*Cyc)*(1-Deg_t)^(n-1);
    term_5 = (RT * (1 - self_dis) * E_in) / (1 + r)^(n+N_c-1); % total CAPEX USD / year
    dis = dis + term_5;
end

fprintf('               ~~~\n');

% To display with commas, you can use the following approach
formatted_cost_5 = num2str(dis, '%.2f');
formatted_cost_with_commas_5 = regexprep(formatted_cost_5, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Total energy discharge: %s MWh\n', formatted_cost_with_commas_5);
formatted_cost_6 = num2str(Cap_e, '%.2f');
formatted_cost_with_commas_6 = regexprep(formatted_cost_6, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Nominal energy capacity: %s MWh\n', formatted_cost_with_commas_6);

%% CALCULATION OF COMPONENTS OF LCOS

L_CAPEX = capex / dis;
L_OM = om / dis;
L_CH = ch / dis;
L_EOL = eol / dis;

fprintf('               ~~~\n');

formatted_cost_7 = num2str(L_CAPEX, '%.2f');
formatted_cost_with_commas_7 = regexprep(formatted_cost_7, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Levelised investment cost: $%s / MWh\n', formatted_cost_with_commas_7);

formatted_cost_8 = num2str(L_OM, '%.2f');
formatted_cost_with_commas_8 = regexprep(formatted_cost_8, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Levelised O&M cost: $%s / MWh\n', formatted_cost_with_commas_8);

formatted_cost_9 = num2str(L_CH, '%.2f');
formatted_cost_with_commas_9 = regexprep(formatted_cost_9, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Levelised charging cost: $%s / MWh\n', formatted_cost_with_commas_9);

formatted_cost_10 = num2str(L_EOL, '%.2f');
formatted_cost_with_commas_10 = regexprep(formatted_cost_10, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Levelised end-of-life cost: $%s / MWh\n', formatted_cost_with_commas_10);

fprintf('               ~~~\n');

%% CALCULATING LEVELISED COST OF STORAGE
LCOS = L_CAPEX + L_OM + L_CH + L_EOL ;

formatted_cost_11 = num2str(LCOS, '%.2f');
formatted_cost_with_commas_11 = regexprep(formatted_cost_11, '\d(?=(\d{3})+\.)', '$&,');
fprintf('Levelised cost of storage (LCOS): $%s / MWh\n', formatted_cost_with_commas_11);

% Define the four values
value1 = L_CAPEX;
value2 = L_OM;
value3 = L_CH;
value4 = L_EOL;

% Create an array of the values
values = [value1, value2, value3, value4];
labels = {'Investment', 'O&M', 'Charging', 'End-of-life'};
colors = [0.6, 0.8, 0.6; 1, 0.8, 0.8; 0.8, 0.6, 0.8; 1, 1, 0.8]; % Define custom colors

% Create a figure and plot the stacked bar chart
figure;
b = bar(1, values, 'stacked');

% Apply colors to the bars
for k = 1:length(b)
    b(k).FaceColor = colors(k, :);
end

% Label the y-axis
ylabel('LCOS in $/MWh');

% Remove x-axis labels
set(gca, 'XTick', []);

% Add title to the plot
title('LCOS Components');

% Add a legend
legend(labels, 'Location', 'eastoutside');

% Optionally, add value labels to each bar component on the side of the bar
hold on;
y_offset = 0; % Initial offset for the label
for i = 1:length(values)
    y_offset = y_offset + values(i); % Increment the offset
    text(1.1, y_offset - values(i)/2, sprintf('$%.2f', values(i)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
hold off;

% Set the x-axis limit to provide space for the labels
xlim([0.5 1.5]);