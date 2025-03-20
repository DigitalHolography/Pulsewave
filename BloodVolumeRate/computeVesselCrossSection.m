function [D, D_std, A, A_std, c1, c2, rsquare] = computeVesselCrossSection(subImg, figName, TB)

% Parameters
params = TB.getParams;
px_size = params.cropSection_pixelSize / (2 ^ params.k);

% Compute velocity profile
profile = mean(subImg, 1);
L = length(profile);

if mean(profile) > 0
    central_range = find(profile > 0.1 * max(profile));
else
    central_range = find(profile < 0.1 * min(profile));
end

centt = mean(central_range);

r_range = (central_range - centt) * px_size;
[p1, p2, p3, rsquare, p1_err, p2_err, p3_err] = customPoly2Fit(r_range', profile(central_range)');
[r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err);

if isnan(r1)
    c1 = 2;
    c2 = L - 1;
else
    c1 = round(centt + (r1 / px_size));
    c2 = round(centt + (r2 / px_size));
    c3 = max(c1, c2);
    c1 = max(min(c1, c2), 2);
    c2 = min(c3, L - 1);
end

% Determine cross-section width
if rsquare < 0.6 || isnan(r1) || isnan(r2)
    D = mean(sum(subImg ~= 0, 2));
    D_std = std(sum(subImg ~= 0, 2));
else
    D = abs(r1 - r2) / px_size;
    D_std = sqrt(r1_err ^ 2 + r2_err ^ 2) / px_size;
end

% Compute cross-sectional area
A = pi * (px_size / 2) ^ 2 * D ^ 2;
A_std = pi * (px_size / 2) ^ 2 * sqrt(D_std ^ 4 + 2 * D_std ^ 2 * D ^ 2);

% Calculate x-axis values (position in µm)
r_ = ((1:L) - centt) * px_size * 1000;

% Calculate standard deviation and confidence interval
stdprofile = std(subImg, [], 1);
curve1 = profile + 0.5 * stdprofile;
curve2 = profile - 0.5 * stdprofile;

% Create figure
f = figure('Visible', 'off');

% Plot confidence interval
Color_std = [0.7, 0.7, 0.7]; % Gray color for confidence interval
fill([r_, fliplr(r_)], [curve1, fliplr(curve2)], Color_std, 'EdgeColor', 'none');
hold on;

% Plot upper and lower bounds of confidence interval
plot(r_, curve1, "Color", Color_std, 'LineWidth', 2);
plot(r_, curve2, "Color", Color_std, 'LineWidth', 2);

% Plot measured data points used for fitting
plot(r_(central_range), profile(central_range), 'xk', 'MarkerSize', 10, 'LineWidth', 1.5);

% Plot zero line
yline(0, '--k', 'LineWidth', 1);

% Plot Poiseuille fit
x_fit = linspace(min(r_), max(r_), 100) / 1000; % Interpolate for smooth fit
y_fit = p1 * x_fit .^ 2 + p2 * x_fit + p3;
plot(x_fit * 1000, y_fit, 'k', 'LineWidth', 1.5);

% Plot vessel boundaries
plot(r1 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(r2 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(linspace(r1 * 1000, r2 * 1000, 10), repmat(-2, 10), '-k', 'LineWidth', 1.5);

% Adjust axes
axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), - 5, axP(4) * 1.07])
box on
set(gca, 'LineWidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])

% Add labels and title
fontsize(gca, 12, "points");
xlabel('Position (µm)');
ylabel('Velocity (mm/s)');
title('velocity profile and laminar flow model fit');

% Save figure

exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', 'projection', ...
    sprintf('%s_proj_poiseuille_%s.png', TB.main_foldername, figName)))

% Close figure
close(f);
end
