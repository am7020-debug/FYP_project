clc;

%% DEFINING THE SHORTAGE OF ENERGY DURING 30 DAY PERIOD

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
a_base = 1; % GWh in 2021
year = 2016;
target_year = 2050;
max_period = target_year - year;
growth_rate = 0.41; % 0.26 to 0.57 based on literatures from WoodMac and BNEF

%% READING DATA IN EXCEL TO CALCULATE MARKET GROWTH OF EACH TECHNOLOGY

% Load the Excel file
filename = 'technology_spec.xlsx'; % Specify your Excel file name

% Get sheet names
[~, sheet_names] = xlsfinfo(filename);

% Initialize a cell array to store the results
results_1 = cell(numel(sheet_names), 4); % Preallocate for sheet name, capacity, and product price
product_prices_1 = cell(numel(sheet_names), 1); % Preallocate for product prices

for i = 1:numel(sheet_names)
    sheetname = sheet_names{i};
    
    % Load data from the current sheet
    data = readtable(filename, 'Sheet', sheetname);
    
    % Display the content of the current sheet to inspect its structure
    disp(['Content of sheet: ', sheetname]);
    disp(data);
    
    % Extract a_base_tech and year_tech dynamically based on parameter names
    a_base_tech_idx_1 = find(strcmp(data.parameter, 'installed capacity in reference year'));
    year_tech_idx_1 = find(strcmp(data.parameter, 'reference year'));
    
    % Extract values
    a_base_tech_1 = data.value(a_base_tech_idx_1); % GWh in base year
    year_tech_1 = data.value(year_tech_idx_1); % Base year

    % Check if a_base_tech_1 and year_tech_1 were successfully extracted
    if isempty(a_base_tech_1) || isempty(year_tech_1)
        error('Unable to extract necessary values from sheet: %s', sheetname);
    end

    % Define the parameters
    target_year = 2050;
    max_period_tech_1 = target_year - year_tech_1;
    growth_rate = 0.41; % 0.26 to 0.57 based on literatures from WoodMac and BNEF
    
    %% DEFINING CALCULATION FOR SPECIFIED MARKET GROWTH
    
    % Target cumulative market capacity in 2050
    target_cumulative_capacity = shortage_GWh; % GWh
    
    % Use fminsearch to find the optimal a_sat_1
    options_1 = optimset('Display', 'iter');
    optimal_a_sat_1 = fminsearch(@(a_sat_1) abs(calculate_cumulative_capacity_1(a_base_tech_1, growth_rate, max_period_tech_1, a_sat_1) - target_cumulative_capacity), 100, options_1);
    
    % Calculate the annual and cumulative market capacities with the optimal a_sat_1
    annual_market_capacity_1 = zeros(max_period_tech_1 + 1, 1);
    for n_period_1 = 0:max_period_tech_1
        annual_market_capacity_1(n_period_1 + 1) = optimal_a_sat_1 / (1 + ((optimal_a_sat_1 - a_base_tech_1) / a_base_tech_1) * exp(-growth_rate * n_period_1));
    end
    cumulative_market_capacity_1 = cumsum(annual_market_capacity_1);
    
    % Define the years from base year to 2050
    years_1 = year_tech_1:(year_tech_1 + max_period_tech_1);
    
    %% SAVING DATA FROM MARKET GROWTH TO AN EXCEL FILE
    
    % Create a table for the results
    results_table_1 = table(years_1', annual_market_capacity_1, cumulative_market_capacity_1, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});
    
    % Save the table to an Excel file
    writetable(results_table_1, ['market_growth_curves_', sheetname, '.xlsx']);
    
    %% PLOTTING DATA FROM MARKET GROWTH 
    
    % Plot the market growth results
    figure;
    subplot(2, 1, 1);
    plot(years_1, annual_market_capacity_1, '-o');
    xlabel('Year');
    ylabel('Annual Market Capacity (GWh)');
    title(['Annual Market Capacity from ', num2str(year_tech_1), ' to 2050 for ', sheetname]);
    grid on;
    
    subplot(2, 1, 2);
    plot(years_1, cumulative_market_capacity_1, '-o');
    xlabel('Year');
    ylabel('Cumulative Market Capacity (GWh)');
    title(['Cumulative Market Capacity from ', num2str(year_tech_1), ' to 2050 for ', sheetname]);
    grid on;
    
    % Display the results
    disp(['Optimal a_sat_1 for ', sheetname, ':']);
    disp(optimal_a_sat_1);
    disp(['Annual and Cumulative Market Capacity from ', num2str(year_tech_1), ' to 2050 for ', sheetname, ':']);
    disp(results_table_1);
    
    %% READING DATA ABOUT EACH TECHNOLOGY TO DETERMINE EXPERIENCE CURVE

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
    b = log(1 - (ER)) / log(2); % Experience curve exponent
    a = p_1 / ( a_base_tech_1 ^ b);

    % Define the capacity range
    x_min = 0.1; % GWh
    x_max = 100000; % GWh

    % Create a vector of capacities
    x = logspace(log10(x_min), log10(x_max), 1000); % Logarithmic scale for better visualization

    % Calculate the experience curve values
    exp_curve = a * x .^ b;

    % Store results in the cell array
    results_1{i, 1} = sheetname;
    results_1{i, 2} = x';
    results_1{i, 3} = exp_curve';

    % Plot the experience curve for the current sheet
    figure;
    loglog(x, exp_curve, '-o');
    xlabel('Capacity (GWh)');
    ylabel('Cost (or similar metric)');
    title(['Experience Curve for ', sheetname]);
    grid on;

    % Define the experience curve function
    experience_curve = @(C, a, b) a .* C .^ b;

    % Initialize a vector to store the product price over the years
    product_price = zeros(length(years_1), 1);

    % Loop through each year and calculate the product price
    for j = 1:length(years_1)
        C = cumulative_market_capacity_1(j); % Cumulative market capacity at year j
        product_price(j) = experience_curve(C, a, b);
    end

    % Store results in the cell array
    results_1{i, 4} = [years_1', product_price];
end

%% FUNCTION DEFINITIONS

function cumulative_capacity_1 = calculate_cumulative_capacity_1(a_base_1, growth_rate, max_period_1, a_sat_1)
    annual_market_capacity_1 = zeros(max_period_1 + 1, 1);
    for n_period_1 = 0:max_period_1
        annual_market_capacity_1(n_period_1 + 1) = a_sat_1 / (1 + ((a_sat_1 - a_base_1) / a_base_1) * exp(-growth_rate * n_period_1));
    end
    cumulative_capacity_1 = sum(annual_market_capacity_1);
end

%% SAVING THE RESULTS TO EXCEL FILE

% Save the results to a new
