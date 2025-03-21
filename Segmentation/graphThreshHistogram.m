function graphThreshHistogram(R, thresholds, mask, colors, name)
%FIG_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
% Set the threshold

ToolBox = getGlobalToolBox;
numLevel = size(thresholds, 2);
numColors = size(colors, 1);

if contains(name, "arter")
    cmap = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
elseif contains(name, "vein")
    cmap = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
else
    cmap = cmapLAB(256, [0 0 0], 0, [1 1 1], 1);
end

figure('Visible', 'off')
imagesc(R .* mask)
title('Map with Colorbar');
colormap(cmap)
colorbar
axis off
axis image
exportgraphics(gca, fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s_MapColorbar.png', name))))

if size(colors, 1) ~= numLevel + 1
    error('Wrong size for colors, there should be N+1 colors for N levels\n numLevel = %d and numColors = %d', numLevel, numColors)
end

m = min(R(mask));
M = max(R(mask));
thresholds = [m thresholds M];

% Bin the data and count occurrences
edges = linspace(m, M, 50); % Set bin edges (modify as needed)
[counts, centers] = histcounts(R(mask & R ~= 0), edges);
counts = [counts 0];

% Plot the histogram with different colors based on threshold
figure (Visible = "off");
hold on;

for ii = 1:(numLevel + 1)

    % Separate bins based on threshold
    inside_thresholds = centers > thresholds(ii) & centers <= thresholds(ii + 1); % Logical array for bins inside thresholds

    % Bars below threshold

    hsvColor = rgb2hsv(colors(ii, :));

    if hsvColor(2) < 0.5
        bar(centers(inside_thresholds), counts(inside_thresholds), 'FaceColor', colors(ii, :), 'EdgeColor', 'black');
    else
        bar(centers(inside_thresholds), counts(inside_thresholds), 'FaceColor', colors(ii, :), 'EdgeColor', 'none');
    end

end

% thresholds
for i = 2:(numLevel + 1)
    xline(thresholds(i), 'k--', 'LineWidth', 2)
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
exportgraphics(gca, fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s_Histo.png', name))))
exportgraphics(gca, fullfile(ToolBox.path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s_Histo.eps', name))))

end
