clc;

% Defining system features
system_capacity = 100; % System capacity in GW
avg_CF = 0.6; % Average capacity factor, 60%
red_CF = 0.1; % Reduced capacity factor, 10%
hours = 720; % Hours in a 30-day period

% Calculating expected and reduced energy generation in GWh
exp_energy = system_capacity * avg_CF * hours;
red_energy = system_capacity * red_CF * hours;

% Calculating overall shortage in GWh and converting to TWh
shortage_GWh = exp_energy - red_energy;

% Define the parameters
a_base = 0.1; % GWh in 2021
year = 2012;
target_year = 2050;
max_period = target_year - year;
growth_rate = 0.41; % 0.26 to 0.57 based on literatures from WoodMac and BNEF

% Target cumulative market capacity in 2050
target_cumulative_capacity = shortage_GWh; % GWh

% Define a function to calculate the cumulative market capacity for a given a_sat
function cumulative_capacity = calculate_cumulative_capacity(a_base, growth_rate, max_period, a_sat)
    annual_market_capacity = zeros(max_period + 1, 1);
    for n_period = 0:max_period
        annual_market_capacity(n_period + 1) = a_sat / (1 + ((a_sat - a_base) / a_base) * exp(-growth_rate * n_period));
    end
    cumulative_capacity = sum(annual_market_capacity);
end

% Use fminsearch to find the optimal a_sat
options = optimset('Display', 'iter');
optimal_a_sat = fminsearch(@(a_sat) abs(calculate_cumulative_capacity(a_base, growth_rate, max_period, a_sat) - target_cumulative_capacity), 100, options);

% Calculate the annual and cumulative market capacities with the optimal a_sat
annual_market_capacity = zeros(max_period + 1, 1);
for n_period = 0:max_period
    annual_market_capacity(n_period + 1) = optimal_a_sat / (1 + ((optimal_a_sat - a_base) / a_base) * exp(-growth_rate * n_period));
end
cumulative_market_capacity = cumsum(annual_market_capacity);

% Define the years from 2012 to 2050
years = year:(year + max_period);

% Create a table for the results
results_table = table(years', annual_market_capacity, cumulative_market_capacity, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});

% Save the table to an Excel file
writetable(results_table, 'experience_and_market_growth_curves_.xlsx');

% Plot the results
figure;
subplot(2, 1, 1);
plot(years, annual_market_capacity, '-o');
xlabel('Year');
ylabel('Annual Market Capacity (GWh)');
title('Annual Market Capacity from 2012 to 2050');
grid on;

subplot(2, 1, 2);
plot(years, cumulative_market_capacity, '-o');
xlabel('Year');
ylabel('Cumulative Market Capacity (GWh)');
title('Cumulative Market Capacity from 2012 to 2050');
grid on;

% Display the results
disp('Optimal a_sat:');
disp(optimal_a_sat);
disp('Annual and Cumulative Market Capacity from 2012 to 2050:');
disp(results_table);

% Load the Excel file
filename = 'technology_spec.xlsx'; % Specify your Excel file name

% Get sheet names
[~, sheet_names] = xlsfinfo(filename);

% Initialize a cell array to store the results
results = cell(numel(sheet_names), 3); % Preallocate for sheet name, capacity, and experience curve

% Loop through each sheet
for i = 1:numel(sheet_names)
    sheetname = sheet_names{i};
    
    % Load data from the current sheet
    data = readtable(filename, 'Sheet', sheetname);
    
    % Display the content of the current sheet to inspect its structure
    disp(['Content of sheet: ', sheetname]);
    disp(data);
    
    % Extract p_1 and ER dynamically based on parameter names
    % Find the indices of 'p_1' and 'experience rate'
    p_1_idx = find(strcmp(data.parameter, 'p_1'));
    ER_idx = find(strcmp(data.parameter, 'experience rate'));
    
    % Extract p_1 and ER
    p_1 = data.value(p_1_idx); % Assuming p_1 value is in the 'value' column
    ER = data.value(ER_idx) / 100; % Assuming ER value is in the 'value' column and converting percentage to a fraction
    
    % Check if p_1 and ER were successfully extracted
    if isempty(p_1) || isempty(ER)
        error('Unable to extract p_1 or ER from sheet: %s', sheetname);
    end
    
    % Define the parameters
    a = p_1;
    b = log(1 - ER) / log(2); % Experience curve exponent
    
    % Define the capacity range
    x_min = 0.1; % GWh
    x_max = 100000; % GWh
    
    % Create a vector of capacities
    x = logspace(log10(x_min), log10(x_max), 1000); % Logarithmic scale for better visualization
    
    % Calculate the experience curve values
    exp_curve = a * x .^ b;
    
    % Store results in the cell array
    results{i, 1} = sheetname;
    results{i, 2} = x';
    results{i, 3} = exp_curve';
    
    % Plot the experience curve for the current sheet
    figure;
    loglog(x, exp_curve, '-o');
    xlabel('Capacity (GWh)');
    ylabel('Cost (or similar metric)');
    title(['Experience Curve for ', sheetname]);
    grid on;
end

% Save the results to a new Excel file
output_filename = 'experience_and_market_growth_curves_.xlsx';
for i = 1:numel(sheet_names)
    sheetname = results{i, 1};
    T = table(results{i, 2}, results{i, 3}, 'VariableNames', {'Capacity', 'ExperienceCurve'});
    writetable(T, output_filename, 'Sheet', sheetname);
end

% Display the results
disp('Experience Curve Data for each technology has been saved to experience_curves.xlsx');
