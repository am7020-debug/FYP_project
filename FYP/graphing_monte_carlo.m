clear;
clc;

% Define the folder containing the data
dataFolder = '/MATLAB Drive/FYP';

% List of data files and corresponding technology names
dataFiles = {'lifetime_costs_histogram_red.fig', 'lifetime_costs_histogram_pum.fig', 'lifetime_costs_histogram_lith.fig', 'lifetime_costs_histogram_hyd.fig', };
techNames = {'Flow Battery', 'Pumped Hydro', 'Lithium Ion Battery', 'Hydrogen'};

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
likelihood_hydro_vs_pumped = calc_likelihood(pdfData{4}.y, pdfData{4}.x, pdfData{2}.y, pdfData{2}.x);
likelihood_hydro_vs_flow = calc_likelihood(pdfData{4}.y, pdfData{4}.x, pdfData{1}.y, pdfData{1}.x);
likelihood_flow_vs_lithium = calc_likelihood(pdfData{1}.y, pdfData{1}.x, pdfData{3}.y, pdfData{3}.x);

% Display the likelihoods as percentages
fprintf('The likelihood that the Hydrogen distribution is lower than the Pumped Hydro distribution is: %.2f%%\n', likelihood_hydro_vs_pumped * 100);
fprintf('The likelihood that the Hydrogen distribution is lower than the Flow Battery distribution is: %.2f%%\n', likelihood_hydro_vs_flow * 100);
fprintf('The likelihood that the Flow Battery distribution is lower than the Lithium Ion Battery distribution is: %.2f%%\n', likelihood_flow_vs_lithium * 100);

% Plot the data on one graph
figure;
hold on;
colors = {'r', 'g', 'b', 'k'}; % Red, Green, Blue, Black for different datasets
for i = 1:length(data)
    if ~isempty(data{i})
        % Plot histogram data
        stairs(data{i}.x, [data{i}.y 0], 'Color', colors{i}, 'DisplayName', techNames{i});
    end
end
hold off;
legend show;
title('Comparing Monte Carlo Analysis of each Technology LCOS');
xlabel('Cost Projections');
ylabel('Probability');
grid on;
xlim([-15000 40000]);
