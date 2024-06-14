clc;
clear;

% Load the Excel file
filename = 'data_collection_1.xlsx'; % Specify your Excel file name

% Check if the file exists
if ~isfile(filename)
    error('File %s not found.', filename);
end

% Initialize a variable to store sheet names
sheet_names = {};

% Try to get sheet names
try
    [~, sheet_names] = xlsfinfo(filename);
catch ME
    disp('Error reading the Excel file:');
    disp(ME.message);
    return;
end

% Display the sheet names to verify correct file access
disp('Sheet names in the file:');
disp(sheet_names);

% Ensure sheet_names is a cell array
if ~iscell(sheet_names)
    sheet_names = cellstr(sheet_names);
end

% Initialize a cell array to store the results
cost_projection_2050 = {};

% Loop through each sheet
for i = 1:numel(sheet_names)
    sheetname = sheet_names{i};
    
    % Check if the sheet name ends with '_CostProjection'
    if endsWith(sheetname, '_CostProjection')
        disp(['Processing sheet: ', sheetname]);
        
        % Load data from the current sheet
        try
            data = readtable(filename, 'Sheet', sheetname);
        catch ME
            disp(['Error reading sheet: ', sheetname]);
            disp(ME.message);
            continue;
        end
        
        % Display the first few rows of the sheet to inspect its structure
        disp(['First few rows of sheet: ', sheetname]);
        disp(head(data));
        
        % Check if the required columns exist
        if ismember('Year', data.Properties.VariableNames) && ismember('ProductPrice', data.Properties.VariableNames)
            % Find the row where the year is 2050
            row_2050 = find(data.Year == 2050);
        
            % Check if year 2050 is found
            if ~isempty(row_2050)
                % Extract the cost projection for the year 2050
                cost_2050 = data.ProductPrice(row_2050);
            
                % Append the technology name and cost projection to the results
                tech_name = erase(sheetname, '_CostProjection');
                cost_projection_2050 = [cost_projection_2050; {tech_name, cost_2050}];
            else
                disp(['Year 2050 not found in sheet: ', sheetname]);
            end
        else
            disp(['Sheet ', sheetname, ' does not contain the required columns.']);
        end
    end
end

% Ensure there are results to display and save
if isempty(cost_projection_2050)
    disp('No cost projection data for 2050 found.');
    return;
end

% Convert the results to a table
try
    cost_projection_2050_table = cell2table(cost_projection_2050, 'VariableNames', {'Technology', 'CostProjection2050'});
catch ME
    disp('Error creating the table:');
    disp(ME.message);
    return;
end

% Display the results
disp('Cost Projection for the Year 2050:');
disp(cost_projection_2050_table);

% Save the results to a new sheet in the Excel file
try
    writetable(cost_projection_2050_table, filename, 'Sheet', 'CostProjection2050Summary');
catch ME
    disp('Error writing to the Excel file:');
    disp(ME.message);
end

% Save the table to a new file for user to download
output_filename = 'CostProjection2050Summary.xlsx';
try
    writetable(cost_projection_2050_table, output_filename);
    disp(['Cost projection summary saved to ', output_filename]);
catch ME
    disp('Error saving the summary table to a new file:');
    disp(ME.message);
end

% Convert the cost projections to numeric array if they are not already
if iscell(cost_projection_2050_table.CostProjection2050)
    costs = cellfun(@(x) str2double(x), cost_projection_2050_table.CostProjection2050);
else
    costs = cost_projection_2050_table.CostProjection2050;
end

% Extract costs for 'Electroliser' and 'Fuel Cell'
electrolyser_idx = strcmp(cost_projection_2050_table.Technology, 'Electroliser');
fuel_cell_idx = strcmp(cost_projection_2050_table.Technology, 'Fuel Cell');

if any(electrolyser_idx)
    electrolyser_cost = costs(electrolyser_idx);
else
    electrolyser_cost = 0;
    disp('"Electroliser" cost not found.');
end

if any(fuel_cell_idx)
    fuel_cell_cost = costs(fuel_cell_idx);
else
    fuel_cell_cost = 0;
    disp('"Fuel Cell" cost not found.');
end

% Sum of the costs for Hydrogen
hydrogen_cost = [electrolyser_cost, fuel_cell_cost];

% Extract costs for other technologies
other_technologies = {'Pumped Hydro', 'Lithium Ion', 'Flow Battery'};
other_costs = zeros(length(other_technologies), 1);

for i = 1:length(other_technologies)
    idx = strcmp(cost_projection_2050_table.Technology, other_technologies{i});
    if any(idx)
        other_costs(i) = costs(idx);
    else
        disp([other_technologies{i}, ' cost not found.']);
    end
end

% Prepare data for plotting
all_technologies = [{'Hydrogen'}; other_technologies'];
all_costs = [sum(hydrogen_cost); other_costs];

% Plotting the bar chart with stacked bar for Hydrogen and simple bars for other technologies
figure;
hold on;
bar(1, hydrogen_cost, 'stacked'); % Stacked bar for Hydrogen

% Simple bars for other technologies
for i = 2:length(all_costs)
    bar(i, all_costs(i), 'FaceColor', [0.5 0.5 0.5]);
end
hold off;

% Adding legend
legend({'Electroliser', 'Fuel Cell', 'Other Technologies'}, 'Location', 'northwest');

xlabel('Technology');
ylabel('Cost Projection for 2050 (USD/kWh)');
title('Cost Projections for 2050 by Technology');
grid on;
xticks(1:length(all_technologies));
xticklabels(all_technologies);
xtickangle(45); % Rotate the x-axis labels for better readability

% Save the plot as an image
saveas(gcf, 'CostProjection2050Plot.png');

% Save the plot as a .fig file for future editing
savefig('CostProjection2050Plot.fig');
