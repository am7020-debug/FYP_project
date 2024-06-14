clear;
clc;

% Define the folder containing the data
dataFolder = '/MATLAB Drive/FYP';

% List of data files and corresponding technology names
dataFiles = {'lifetime_costs_histogram_red.fig', 'lifetime_costs_histogram_pum.fig', 'lifetime_costs_histogram_lith.fig', 'lifetime_costs_histogram_hyd.fig', 'lifetime_costs_histogram_caes.fig'};
techNames = {'Flow Battery', 'Pumped Hydro', 'Lithium Ion Battery', 'Hydrogen', 'CAES'};

% Initialize variables to store data
data = cell(length(dataFiles), 1);

% Load the data from each figure file
for i = 1:length(dataFiles)
    dataPath = fullfile(dataFolder, dataFiles{i});
    if isfile(dataPath)
        fig = openfig(dataPath, 'invisible'); % Open the figure invisibly
        ax = gca; % Get the current axes of the figure
        
        % Try to find histogram objects
        dataObjs = findobj(ax, 'Type', 'histogram'); 
        
        % Check if any histogram objects are found
        if ~isempty(dataObjs)
            % Assuming the first histogram object contains the data we need
            data{i}.x = dataObjs(1).BinEdges;
            data{i}.y = dataObjs(1).Values;
        else
            error('No suitable graphical objects found in figure %s', dataFiles{i});
        end
        
        close(fig); % Close the figure
    else
        error('File %s not found', dataFiles{i});
    end
end

% Function to convert histogram data to pdf
function [x, pdf] = histogram_to_pdf(edges, values)
    bin_widths = diff(edges);
    pdf = values ./ sum(values) ./ bin_widths;
    x = edges(1:end-1) + bin_widths / 2;
end

% Convert the histogram data to PDFs
pdfData = cell(size(data));
for i = 1:length(data)
    [pdfData{i}.x, pdfData{i}.y] = histogram_to_pdf(data{i}.x, data{i}.y);
end


% Plot the data for Lithium-Ion Battery and Flow Battery
figure1 = figure;
hold on;
colors = {'b', 'r'}; % Blue, Red for the two datasets
indices = [3, 1]; % Indices for Lithium-Ion Battery and Flow Battery
for i = 1:length(indices)
    idx = indices(i);
    if ~isempty(data{idx})
        % Interpolate and smooth the data
        fine_x = linspace(min(pdfData{idx}.x), max(pdfData{idx}.x), 1000);
        fine_y = interp1(pdfData{idx}.x, pdfData{idx}.y, fine_x, 'pchip');
        % Add zero points to close the area to x-axis
        fine_x = [fine_x(1) fine_x fine_x(end)];
        fine_y = [0 fine_y 0];
        % Plot the filled area under the curve
        fill(fine_x, fine_y, colors{i}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', techNames{idx});
    end
end
hold off;
legend show;
title('Lithium-Ion Battery and Flow Battery LCOS Comparison');
xlabel('Cost Projections');
ylabel('Probability');
grid on;
xlim([min(cellfun(@(x) min(x.x), pdfData(indices))), max(cellfun(@(x) max(x.x), pdfData(indices)))]);
ylim([0, max(cellfun(@(x) max(x.y), pdfData(indices)))]);

% Save the first figure
saveas(figure1, 'Lithium_Ion_vs_Flow.fig');
saveas(figure1, 'Lithium_Ion_vs_Flow.png');

% Plot the data for Pumped Hydro, Hydrogen, and CAES
figure2 = figure;
hold on;
colors = {'g', 'k', 'm'}; % Green, Black, Magenta for the three datasets
indices = [2, 4, 5]; % Indices for Pumped Hydro, Hydrogen, and CAES
for i = 1:length(indices)
    idx = indices(i);
    if ~isempty(data{idx})
        % Interpolate and smooth the data
        fine_x = linspace(min(pdfData{idx}.x), max(pdfData{idx}.x), 1000);
        fine_y = interp1(pdfData{idx}.x, pdfData{idx}.y, fine_x, 'pchip');
        % Add zero points to close the area to x-axis
        fine_x = [fine_x(1) fine_x fine_x(end)];
        fine_y = [0 fine_y 0];
        % Plot the filled area under the curve
        fill(fine_x, fine_y, colors{i}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', techNames{idx});
    end
end
hold off;
legend show;
title('Pumped Hydro, Hydrogen, and CAES LCOS Comparison');
xlabel('Cost Projections [USD/MWh]');
ylabel('Probability');
grid on;
xlim([min(cellfun(@(x) min(x.x), pdfData(indices))), max(cellfun(@(x) max(x.x), pdfData(indices)))]);
ylim([0, max(cellfun(@(x) max(x.y), pdfData(indices)))]);

% Save the second figure
saveas(figure2, 'Pumped_Hydro_Hydrogen_CAES.fig');
saveas(figure2, 'Pumped_Hydro_Hydrogen_CAES.png');

% Plot the data for all five technologies together
figure3 = figure;
hold on;
colors = {'r', 'g', 'b', 'k', 'm'}; % Red, Green, Blue, Black, Magenta for the datasets
indices = 1:5; % Indices for all five technologies
for i = 1:length(indices)
    idx = indices(i);
    if ~isempty(data{idx})
        % Interpolate and smooth the data
        fine_x = linspace(min(pdfData{idx}.x), max(pdfData{idx}.x), 1000);
        fine_y = interp1(pdfData{idx}.x, pdfData{idx}.y, fine_x, 'pchip');
        % Add zero points to close the area to x-axis
        fine_x = [fine_x(1) fine_x fine_x(end)];
        fine_y = [0 fine_y 0];
        % Plot the filled area under the curve
        fill(fine_x, fine_y, colors{i}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', techNames{idx});
    end
end
hold off;
legend show;
title('All Technologies LCOS Comparison');
xlabel('Cost Projections [USD/MWh]');
ylabel('Probability');
grid on;
xlim([min(cellfun(@(x) min(x.x), pdfData(indices))), max(cellfun(@(x) max(x.x), pdfData(indices)))]);
ylim([0, max(cellfun(@(x) max(x.y), pdfData(indices)))]);

% Save the third figure
saveas(figure3, 'All_Technologies.fig');
saveas(figure3, 'All_Technologies.png');

% Function to calculate the likelihood that one distribution is lower than another
function likelihood = calc_likelihood(pdf1, x1, pdf2, x2)
    common_x = linspace(min([x1, x2]), max([x1, x2]), 1000);
    interp_pdf1 = interp1(x1, pdf1, common_x, 'linear', 0);
    interp_pdf2 = interp1(x2, pdf2, common_x, 'linear', 0);
    non_overlap_area = trapz(common_x, max(0, interp_pdf1 - interp_pdf2));
    total_area1 = trapz(common_x, interp_pdf1);
    likelihood = non_overlap_area / total_area1;
end

% Calculate the likelihoods
likelihood_lithium_vs_flow = 1-calc_likelihood(pdfData{3}.y, pdfData{3}.x, pdfData{1}.y, pdfData{1}.x);
likelihood_flow_vs_pumped = 1-calc_likelihood(pdfData{1}.y, pdfData{1}.x, pdfData{2}.y, pdfData{2}.x);
likelihood_pumped_vs_hydrogen = 1-calc_likelihood(pdfData{2}.y, pdfData{2}.x, pdfData{4}.y, pdfData{4}.x);
likelihood_hydrogen_vs_caes = 1-calc_likelihood(pdfData{4}.y, pdfData{4}.x, pdfData{5}.y, pdfData{5}.x);

% Display the likelihoods as percentages
fprintf('The likelihood that the Lithium-Ion Battery distribution is more expensive than the Flow Battery distribution is: %.2f%%\n', (1 - likelihood_lithium_vs_flow) * 100);
fprintf('The likelihood that the Flow Battery distribution is more expensive than the Pumped Hydro distribution is: %.2f%%\n', (1 - likelihood_flow_vs_pumped) * 100);
fprintf('The likelihood that the Pumped Hydro distribution is more expensive than the Hydrogen distribution is: %.2f%%\n', (1 - likelihood_pumped_vs_hydrogen) * 100);
fprintf('The likelihood that the Hydrogen distribution is more expensive than the CAES distribution is: %.2f%%\n', (1 - likelihood_hydrogen_vs_caes) * 100);
