clc;
%% 

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
output_filename = 'experience_curves.xlsx';
for i = 1:numel(sheet_names)
    sheetname = results{i, 1};
    T = table(results{i, 2}, results{i, 3}, 'VariableNames', {'Capacity', 'ExperienceCurve'});
    writetable(T, output_filename, 'Sheet', sheetname);
end

% Display the results
disp('Experience Curve Data for each technology has been saved to experience_curves.xlsx');
