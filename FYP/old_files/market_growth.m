clc

% Define the parameters
a_base = 0.1; % GWh in 2012
year = 2012;
target_year = 2050;
max_period = target_year - year;
growth_rate = 0.41; % 0.26 to 0.57 based on literatures from WoodMac and BNEF

% Target cumulative market capacity in 2050
target_cumulative_capacity = 36000; % GWh

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

% Plot the results
figure;
subplot(2, 1, 1);
plot(years, annual_market_capacity, '-o');
xlabel('Year');
ylabel('Annual Market Capacity (GWh)');
title('Annual Market Capacity');
grid on;

subplot(2, 1, 2);
plot(years, cumulative_market_capacity, '-o');
xlabel('Year');
ylabel('Cumulative Market Capacity (GWh)');
title('Cumulative Market Capacity');
grid on;

% Display the results
disp('Optimal a_sat:');
disp(optimal_a_sat);
disp('Annual and Cumulative Market Capacity:');
disp(table(years', annual_market_capacity, cumulative_market_capacity, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'}));

% Display the results
disp('Annual Market Capacity');
disp(table(years', annual_market_capacity, cumulative_market_capacity, 'VariableNames', {'Year', 'AnnualMarketCapacity', 'CumulativeMarketCapacity'}));
