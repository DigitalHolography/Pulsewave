function figHistogram(R, level, mask, color1, color2, name, ToolBox)
%FIG_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
% Set the threshold
threshold = level(3);

m = min(R(mask));
M = max(R(mask));

% Bin the data and count occurrences
edges = linspace(m, M, 50); % Set bin edges (modify as needed)
[counts, centers] = histcounts(R(mask), edges);
counts = [counts 0];

% Separate bins based on threshold
below_thresh = centers < threshold; % Logical array for bins below threshold
above_thresh = centers >= threshold; % Logical array for bins above threshold

% Plot the histogram with different colors based on threshold
figure;
hold on;

% Bars below threshold
if strcmp(color1,'white') == 1
    bar(centers(below_thresh), counts(below_thresh), 'FaceColor', color1, 'EdgeColor', 'black');
else
    bar(centers(below_thresh), counts(below_thresh), 'FaceColor', color1, 'EdgeColor', 'none');
end

% Bars above threshold
if strcmp(color2,'white') == 1
    bar(centers(above_thresh), counts(above_thresh), 'FaceColor', color2, 'EdgeColor', 'black');
else
    bar(centers(above_thresh), counts(above_thresh), 'FaceColor', color2, 'EdgeColor', 'none');
end

% thresholds
for i = 2:size(level,2)
xline(level(i), 'k--', 'LineWidth', 2)
end

% Add labels and title
xlabel('Data Value');
ylabel('Frequency');
title('Histogram with Threshold Coloring');
axis tight
set(gca, 'Linewidth', 2)
pbaspect([1.68 1 1])
box on
hold off;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s.png', name))))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s.eps', name))))

end

