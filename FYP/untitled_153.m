clc;

%% DEFINING THE SHORTAGE OF ENERGY DURING 30 DAY PERIOD

% Defining system features
system_capacity = 100; % Wind generation system capacity in GW
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

%% DEFINING CALCULATION FOR MARKET GROWTH GIVEN THE SHORTAGE

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

% Define the years from 2016 to 2050
years = year:(year + max_period);

% Initialize a table to store the general market growth data
general_market_growth_data = table();

%% SAVING DATA FROM MARKET GROWTH TO AN EXCEL FILE

% Create a table for the results
results_table = table(years', annual_market_capacity, cumulative_market_capacity, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});

% Save the table to an Excel file
writetable(results_table, 'data_collection.xlsx');

%% PLOTTING DATA FROM MARKET GROWTH 

% Plot the market growth results
figure;
subplot(2, 1, 1);
plot(years, annual_market_capacity, '-o');
xlabel('Year');
ylabel('Annual Market Capacity (GWh)');
title('Annual Market Capacity from 2016 to 2050');
grid on;

subplot(2, 1, 2);
plot(years, cumulative_market_capacity, '-o');
xlabel('Year');
ylabel('Cumulative Market Capacity (GWh)');
title('Cumulative Market Capacity from 2016 to 2050');
grid on;

% Display the results
disp('Optimal a_sat:');
disp(optimal_a_sat);
disp('Annual and Cumulative Market Capacity from 2016 to 2050:');
disp(results_table);

%% READING DATA IN EXCEL TO CALCULATE MARKET GROWTH OF EACH TECHNOLOGY

% Load the Excel file
filename = 'technology_spec_2.xlsx'; % Specify your Excel file name

% Get sheet names
[~, sheet_names] = xlsfinfo(filename);

% Initialize a cell array to store the results
results_1 = cell(numel(sheet_names), 4); % Preallocate for sheet name, capacity, and product price
product_prices_1 = cell(numel(sheet_names), 1); % Preallocate for product prices

for i = 1:numel(sheet_names)
    sheetname = sheet_names{i};
    
    % Load data from the current sheet
    data = readtable(filename, 'Sheet', sheetname);
    
    % Display the column names to inspect their structure
    disp(['Column names in sheet ', sheetname, ':']);
    disp(data.Properties.VariableNames);

    % Convert the 'parameter' column to string for easier comparison
    param_col_name = data.Properties.VariableNames{1}; % Assume 'parameter' is the first column
    value_col_name = data.Properties.VariableNames{3}; % Assume 'value' is the third column

    data.(param_col_name) = string(data.(param_col_name));

    % Extract a_base_tech, year_tech, and final energy storage required dynamically based on parameter names
    a_base_tech_idx_1 = find(strcmpi(data.(param_col_name), 'installed capacity in reference year'));
    year_tech_idx_1 = find(strcmpi(data.(param_col_name), 'reference year'));
    final_storage_idx = find(strcmpi(data.(param_col_name), 'final energy storage required'));
    
    % Extract values
    a_base_tech_1 = data.(value_col_name)(a_base_tech_idx_1); % GWh in base year
    year_tech_1 = data.(value_col_name)(year_tech_idx_1); % Base year
    target_cumulative_capacity = data.(value_col_name)(final_storage_idx); % Final energy storage required

    % Check if a_base_tech_1, year_tech_1, and final_storage_idx were successfully extracted
    if isempty(a_base_tech_1) || isempty(year_tech_1) || isempty(target_cumulative_capacity)
        error('Unable to extract necessary values from sheet: %s', sheetname);
    end

    % Define the parameters
    target_year = 2050;
    max_period_tech_1 = target_year - year_tech_1;
    growth_rate = 0.41; % 0.26 to 0.57 based on literatures from WoodMac and BNEF
    
    %% DEFINING CALCULATION FOR SPECIFIED MARKET GROWTH
    
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
    years_1 = (year_tech_1:(year_tech_1 + max_period_tech_1))';

    %% SAVING DATA FROM MARKET GROWTH TO AN EXCEL FILE
    
    % Create a table for the results
    results_table_1 = table(years_1, annual_market_capacity_1, cumulative_market_capacity_1, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});
    
    % Save the table to an Excel file
    writetable(results_table_1, 'data_collection.xlsx', 'Sheet', [sheetname, '_MarketGrowth']);
    
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
    p_1_idx = find(strcmpi(data.(param_col_name), 'p_1'));
    ER_idx = find(strcmpi(data.(param_col_name), 'experience rate'));

    % Extract p_1 and ER
    p_1 = data.(value_col_name)(p_1_idx); % Assuming p_1 value is in the 'value' column
    ER = data.(value_col_name)(ER_idx) / 100; % Assuming ER value is in the 'value' column and converting percentage to a fraction

    % Check if p_1 and ER were successfully extracted
    if isempty(p_1) || isempty(ER)
        error('Unable to extract p_1 or ER from sheet: %s', sheetname);
    end

    % Define the parameters
    b = log(1 - (ER)) / log(2); % Experience curve exponent
    a = p_1 / ( a_base_tech_1 ^ b);
    display (a);
    display (b);
    display (a_base_tech_1);

    % Define the capacity range
    x_min = 0.1; % GWh
    x_max = 100000; % GWh

    % Create a vector of capacities
    x = logspace(log10(x_min), log10(x_max), 1000); % Logarithmic scale for better visualization

    % Calculate the experience curve values
    exp_curve = a * x .^ b;

    % Ensure x and exp_curve are column vectors
    x = x(:);
    exp_curve = exp_curve(:);

    % Store results in the cell array
    results_1{i, 1} = sheetname;
    results_1{i, 2} = x;
    results_1{i, 3} = exp_curve;

    % Plot the experience curve for the current sheet
    figure;
    loglog(x, exp_curve, '-o');
    xlabel('Capacity (GWh)');
    ylabel('Cost (or similar metric)');
    title(['Experience Curve for ', sheetname]);
    grid on;

    % Define the experience curve function
    experience_curve = @(C, a, b) a .* C .^ b;

    %% CONVOLUTION OF MARKET GROWTH AND EXPERIENCE CURVE

    % Initialize a vector to store the product price over the years
    product_price = zeros(length(years_1), 1);

    % Loop through each year and calculate the product price
    for j = 1:length(years_1)
        C = cumulative_market_capacity_1(j); % Cumulative market capacity at year j
        product_price(j) = experience_curve(C, a, b);
    end

    % Store results in the cell array
    results_1{i, 4} = [years_1, product_price];

    %% PLOTTING PRODUCT PRICE PROJECTION

    figure;
    plot(years_1, product_price, '-o');
    xlabel('Year');
    ylabel('Product Price (USD/kWh)');
    ylim([0 1.5 * max(product_price)]); % 50% more than the max value
    title(['Product Price Projection for ', sheetname]);
    grid on;

    %% SAVING DATA TO EXCEL FILE

    % Append data to general market growth data
    temp_table = table(years_1, annual_market_capacity_1, cumulative_market_capacity_1, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});
    general_market_growth_data = [general_market_growth_data; temp_table];
    
    % Save experience rate, cost projection, and market growth for each technology
    exp_rate_table = table(x, exp_curve, 'VariableNames', {'Capacity', 'ExperienceCurve'});
    writetable(exp_rate_table, 'data_collection.xlsx', 'Sheet', [sheetname, '_ExperienceCurve']);
    
    cost_proj_table = table(years_1, product_price, 'VariableNames', {'Year', 'ProductPrice'});
    writetable(cost_proj_table, 'data_collection.xlsx', 'Sheet', [sheetname, '_CostProjection']);
    
    market_growth_table = table(years_1, annual_market_capacity_1, cumulative_market_capacity_1, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'});
    writetable(market_growth_table, 'data_collection.xlsx', 'Sheet', [sheetname, '_MarketGrowth']);
end

%% FUNCTION DEFINITIONS

function cumulative_capacity_1 = calculate_cumulative_capacity_1(a_base_1, growth_rate, max_period_1, a_sat_1)
    annual_market_capacity_1 = zeros(max_period_1 + 1, 1);
    for n_period_1 = 0:max_period_1
        annual_market_capacity_1(n_period_1 + 1) = a_sat_1 / (1 + ((a_sat_1 - a_base_1) / a_base_1) * exp(-growth_rate * n_period_1));
    end
    cumulative_capacity_1 = sum(annual_market_capacity_1);
end

%% PLOTTING ALL COST PROJECTION CURVES ONTO ONE GRAPH WITH A LOGARITHMIC SCALE

figure;
hold on;
for i = 1:numel(sheet_names)
    sheetname = sheet_names{i};
    price_data = results_1{i, 4};
    years_1 = price_data(:, 1);
    product_price = price_data(:, 2);
    plot(years_1, product_price, '-o', 'DisplayName', sheetname);
end
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
xlabel('Year');
ylabel('Product Price (USD/kWh)');
title('Product Price Projections for All Technologies');
legend show;
grid on;
hold off;

