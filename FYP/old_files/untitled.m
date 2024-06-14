

% Plot the distributions
figure;
hold on;
histogram(dist1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'b');
histogram(dist2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('Distribution 1', 'Distribution 2');
title('Random Distributions');
hold off;

% Calculate the percentage chance that one distribution is lower than the other
count_dist1_lower = sum(dist1 < dist2);
percentage_dist1_lower = (count_dist1_lower / length(dist1)) * 100;

% Display the result
fprintf('Percentage chance that Distribution 1 is lower than Distribution 2: %.2f%%\n', percentage_dist1_lower);
