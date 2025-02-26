function poiseuilleProfileFigure(subImg, profile, centt, central_range, p1, p2, p3, r1, r2, rsquare, insert, name_section, ToolBox)
% Generate a figure showing the velocity profile and Poiseuille fit.
%
% Inputs:
%   subImg          - 2D array, the sub-image of the blood vessel.
%   profile         - 1D array, the velocity profile across the vessel.
%   centt           - Scalar, the center of the profile.
%   central_range   - Indices of the central range used for fitting.
%   p1, p2, p3      - Coefficients of the quadratic fit.
%   r1, r2          - Roots of the quadratic fit (vessel boundaries).
%   rsquare         - Coefficient of determination.
%   insert          - String, additional identifier for the filename.
%   name_section    - String, name of the section.
%   ToolBox         - Struct, contains parameters and paths.
%   sectionIdx      - Scalar, index of the current section.

% Get parameters
PW_params = ToolBox.getParams;
k = PW_params.k;

% Calculate x-axis values (position in µm)
r_ = ((1:length(profile)) - centt) * (PW_params.cropSection_pixelSize / 2^k) * 1000;

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
y_fit = p1 * x_fit.^2 + p2 * x_fit + p3;
plot(x_fit * 1000, y_fit, 'k', 'LineWidth', 1.5);

% Plot vessel boundaries
plot(r1 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(r2 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(linspace(r1 * 1000, r2 * 1000, 10), repmat(-2, 10), '-k', 'LineWidth', 1.5);

% Adjust axes
axis tight;
ax = axis;
ax(3) = -5; % Set lower y-axis limit
axis(ax);
box on
set(gca, 'LineWidth', 2)

% Add labels and title
xlabel('Position (µm)');
ylabel('Velocity (mm/s)');
title('Velocity Profile and Poiseuille Fit');

% Add legend
legend('', '', '', 'Measured Data', '', sprintf('fit R² = %d', rsquare), ...
    'Location', 'northeast');

% Save figure
saveas(f, fullfile(ToolBox.PW_path_png, 'projection', ...
    sprintf('%s_%s_proj_poiseuille_%s.png', ToolBox.main_foldername, insert, name_section)));

% Close figure
close(f);

end