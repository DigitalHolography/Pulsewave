function [Q_t, dQ_t] = plotRadius(Q_mat, dQ_mat, fullTime, rad, idx_start, idx_end, name)
% Get global toolbox and parameters
ToolBox = getGlobalToolBox;
numCircles = size(Q_mat, 1); % Number of circles (radii)

% Define color for shaded regions
Color_std = [0.7 0.7 0.7];

% Sum across sections for each circle and compute uncertainties
Q_rt = squeeze(sum(Q_mat, 2, "omitnan")); % Sum over sections
dQ_rt = squeeze(sqrt(sum(dQ_mat .^ 2, 2, "omitnan"))); % Quadratic sum of uncertainties

% Compute time-averaged mean and standard deviation
Q_r = squeeze(mean(Q_rt(:, idx_start:idx_end), 2))'; % Mean over time
N = size(dQ_rt, 2);
dQ_r = squeeze(sqrt(sum(dQ_rt(:, idx_start:idx_end) .^ 2, 2)))' / N; % RMS of uncertainties

% Create shaded region for uncertainty
curve1 = Q_r + dQ_r; % Upper bound
curve2 = Q_r - dQ_r; % Lower bound
rad2 = [rad, fliplr(rad)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

% Plot time-averaged blood volume rate vs. radius
figure("Visible", "off");
fill(rad2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(rad, Q_r, '-k', 'LineWidth', 2); % Mean curve
yline(mean(Q_r), '--k', 'LineWidth', 2); % Horizontal mean line
legend({'', '', '', '', sprintf('Mean = %0.2f µL/min', mean(Q_r))});

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
box on;
ylabel('Blood Volume Rate (µL/min)');
xlabel('Radius (pixels)');
title("Time-Averaged Blood Volume Rate");
set(gca, 'PlotBoxAspectRatio', [1.618 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_mean_%s_radius.png", ToolBox.main_foldername, name)));

% Plot radial variations of blood volume rate over time
figure("Visible", "off");
hold on;

for circleIdx = 1:numCircles
    plot(fullTime, squeeze(mean(Q_mat(circleIdx, :, :), 2, "omitnan")), 'LineWidth', 2);
end

axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
box on;
ylabel('Blood Volume Rate (µL/min)');
xlabel('Time (s)');
title("Radial Variations of Blood Volume Rate");
set(gca, 'PlotBoxAspectRatio', [1.618 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_variance_%s_time.png", ToolBox.main_foldername, name)));

% Compute total blood volume rate over time
Q_t = squeeze(mean(Q_rt, 1)); % Mean over time
N = size(dQ_rt, 1);
dQ_t = squeeze(sqrt(sum(dQ_rt .^ 2, 1))) / N; % RMS of uncertainties

% Compute statistics for the time range
mean_Q = mean(Q_t(idx_start:idx_end)); % Time-averaged mean
max_Q = max(Q_t(idx_start:idx_end)); % Maximum value in the range
N = size(dQ_t(idx_start:idx_end), 1);
mean_dQ = sqrt(sum(dQ_t(idx_start:idx_end) .^ 2)) / N; % RMS of the uncertainty

% Plot total blood volume rate over time
figure("Visible", "off");
curve1 = Q_t + dQ_t; % Upper bound
curve2 = Q_t - dQ_t; % Lower bound
ft2 = [fullTime, fliplr(fullTime)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

hold on;
fill(ft2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
yline(0, 'k-', 'LineWidth', 2); % Zero line
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(fullTime, Q_t, '-k', 'LineWidth', 2); % Mean curve
yline(mean_Q, '--k', 'LineWidth', 2); % Horizontal mean line

% Mark the time range used for averaging
plot(fullTime(idx_start), 1.07 * max_Q, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_end), 1.07 * max_Q, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_start:idx_end), repmat(1.07 * max_Q, 1, idx_end - idx_start + 1), '-k', 'LineWidth', 2);

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);
box on;
ylabel('Blood Volume Rate (µL/min)');
xlabel('Time (s)');
title(sprintf("Total Blood Volume Rate (Avg. %0.2f µL/min)", mean_Q));
set(gca, 'PlotBoxAspectRatio', [1.618 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_allrad_%s_time.png", ToolBox.main_foldername, name)));

% Write results to a text file
fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_main_outputs', '.txt')), 'a');
fprintf(fileID, 'Flow Rate %s : %f (µL/min) \r\n', name, mean_Q);
fprintf(fileID, 'Flow Rate Standard Deviation %s : %f (µL/min) \r\n', name, mean_dQ);
fclose(fileID);
end
