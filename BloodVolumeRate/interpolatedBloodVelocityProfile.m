function interpolatedBloodVelocityProfile(v_cell, dv_cell, sysIdx, diasIdx, numSections, name, rad)
% interpolatedBloodVelocityProfile - Combines and optimizes the computation and plotting of velocity profiles.
% Inputs:
%   v_cell: Cell array containing mean velocity profiles for each circle, section, and frame.
%   dv_cell: Cell array containing velocity profile uncertainties.
%   sysIdx: Indices of systolic frames.
%   diasIdx: Indices of diastolic frames.
%   numSections: Number of sections for each circle.
%   name: Name identifier for saving files.
%   rad: Radius values for each circle.

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = ToolBox.exportVideos;

Color_err_sys = [0.8 0.6 0.6];
Color_err_dias = [0.6 0.6 0.8];
Color_err = [0.7 0.7 0.7];

% Get sizes
numCircles = params.json.BloodVolumeRateAnalysis.NumberOfCircles;
numFrames = size(v_cell{1}, 2);
numInterp = params.json.BloodVolumeRateFigures.numInterp;

% Preallocate interpolated profiles
v_interp = cell(1, numCircles);
dv_interp = cell(1, numCircles);

% Interpolate profiles for all circles, sections, and frames
parfor cIdx = 1:numCircles
    v_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp); % Preallocate for mean velocity
    dv_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp); % Preallocate for uncertainties

    for sectionIdx = 1:numSections(cIdx)

        for frameIdx = 1:numFrames
            v = v_cell{cIdx}{sectionIdx, frameIdx}; % Mean velocity profile
            dv = dv_cell{cIdx}{sectionIdx, frameIdx}; % Velocity uncertainty profile

            % Handle edge case with negative velocities
            if any(v < 0)
                [~, locs] = findpeaks(-v); % Find peaks in negative velocities

                if isempty(locs)
                    indx = find(v > 0); % Use positive velocities if no peaks
                elseif length(locs) > 1
                    indx = locs(1):locs(end); % Use range between first and last peak
                else

                    if locs(1) > length(v) / 2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(v);
                    end

                end

            else
                indx = find(v > 0); % Use positive velocities
            end

            % Interpolate profiles
            L = length(indx);
            v_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, v(indx), linspace(1, L, numInterp)); % Interpolate mean velocity
            dv_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, dv(indx), linspace(1, L, numInterp)); % Interpolate uncertainty
        end

    end

end

% Preallocate variables for systolic and diastolic profiles
v_sys = zeros(numInterp, 1);
dv_sys = zeros(numInterp, 1);
v_dias = zeros(numInterp, 1);
dv_dias = zeros(numInterp, 1);

for cIdx = 1:numCircles

    for sectionIdx = 1:numSections(cIdx)
        % Average profiles for systolic and diastolic frames
        v_sys = v_sys + squeeze(sum(v_interp{cIdx}(sectionIdx, sysIdx, :), 2) / numSections(cIdx));
        dv_sys = dv_sys + squeeze(sum(dv_interp{cIdx}(sectionIdx, sysIdx, :) .^ 2, 2) / numSections(cIdx));
        v_dias = v_dias + squeeze(sum(v_interp{cIdx}(sectionIdx, diasIdx, :), 2) / numSections(cIdx));
        dv_dias = dv_dias + squeeze(sum(dv_interp{cIdx}(sectionIdx, diasIdx, :) .^ 2, 2) / numSections(cIdx));
    end

end

numSys = length(sysIdx);
numDias = length(diasIdx);
v_sys = (v_sys / (numSys * numCircles))';
dv_sys = (sqrt(dv_sys) / (numSys * numCircles))';
v_dias = (v_dias / (numDias * numCircles))';
dv_dias = (sqrt(dv_dias) / (numDias * numCircles))';

% Create curves for plotting
curve1_sys = v_sys + dv_sys;
curve2_sys = v_sys - dv_sys;
inBetween_sys = [curve1_sys, fliplr(curve2_sys)];
ft2_sys = [(1:numInterp), fliplr(1:numInterp)];

curve1_dias = v_dias + dv_dias;
curve2_dias = v_dias - dv_dias;
inBetween_dias = [curve1_dias, fliplr(curve2_dias)];
ft2_dias = [(1:numInterp), fliplr(1:numInterp)];

% Plot systolic and diastolic profiles
figure("Visible", "off");
ax = gca;
hold on;

fill(ax, ft2_sys, inBetween_sys, Color_err_sys, 'EdgeColor', 'none');
plot(ax, 1:numInterp, curve1_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, curve2_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, v_sys, '--r', 'LineWidth', 2);

fill(ax, ft2_dias, inBetween_dias, Color_err_dias, 'EdgeColor', 'none');
plot(ax, 1:numInterp, curve1_dias, 'Color', Color_err_dias, 'LineWidth', 2);
plot(ax, 1:numInterp, curve2_dias, 'Color', Color_err_dias, 'LineWidth', 2);
plot(ax, 1:numInterp, v_dias, '--b', 'LineWidth', 2);

% Poiseuille fit
warning("off");
[~, centt_sys] = max(v_sys);
r_range_sys = (1:numInterp) - centt_sys;
f_sys = fit(r_range_sys', v_sys', 'poly2');
poiseuille_fit_sys = f_sys.p1 * r_range_sys .^ 2 + f_sys.p2 * r_range_sys + f_sys.p3;
poiseuille_fit_sys(poiseuille_fit_sys < 0) = 0;

[~, centt_dias] = max(v_dias);
r_range_dias = (1:numInterp) - centt_dias;
f_dias = fit(r_range_dias', v_dias', 'poly2');
poiseuille_fit_dias = f_dias.p1 * r_range_dias .^ 2 + f_dias.p2 * r_range_dias + f_dias.p3;
poiseuille_fit_dias(poiseuille_fit_dias < 0) = 0;

plot(ax, 1:numInterp, poiseuille_fit_sys, '-r', 'LineWidth', 2);
plot(ax, 1:numInterp, poiseuille_fit_dias, '-b', 'LineWidth', 2);
warning("on");

% Adjust axes and labels
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);
hold off;

box on;
set(gca, 'Linewidth', 2);
xlabel('Profile (AU)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);
pbaspect([1.618 1 1]);

% Export figure
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_diasys_%s.png", ToolBox.main_foldername, name)));

% Interpolated blood velocity profile and video export
if exportVideos
    fig = figure("Visible", "off");
    ax = axes(fig);
    hold(ax, "on");

    fillPlot = fill(ax, NaN, NaN, Color_err, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    curve1Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    curve2Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    meanPlot = plot(ax, NaN, NaN, '-k', 'LineWidth', 2);
    poiseuillePlot = plot(ax, NaN, NaN, '-r', 'LineWidth', 2);

    axis tight;
    ax.YLim = [-10; 30];
    box on;
    set(gca, 'Linewidth', 2);

    parfor cIdx = 1:numCircles
        video = zeros(420, 560, 3, numFrames, 'single'); % Preallocate video array
        title(['interpolated time-averaged velocity profile at radius = ', num2str(rad(cIdx)), ' pix']);

        for frameIdx = 1:numFrames
            % Compute mean and RMS profiles
            mean_v_interp_profile = mean(v_interp{cIdx}(:, frameIdx, :), 1);
            rms_dv_interp_profile = sqrt(mean(dv_interp{cIdx}(:, frameIdx, :) .^ 2, 1));

            curve1 = mean_v_interp_profile + rms_dv_interp_profile;
            curve2 = mean_v_interp_profile - rms_dv_interp_profile;
            inBetween = [curve1, fliplr(curve2)]';
            ft2 = [(1:numInterp), fliplr(1:numInterp)];

            % Update plot elements
            set(fillPlot, 'XData', ft2, 'YData', inBetween);
            set(curve1Plot, 'XData', 1:numInterp, 'YData', curve1);
            set(curve2Plot, 'XData', 1:numInterp, 'YData', curve2);
            set(meanPlot, 'XData', 1:numInterp, 'YData', mean_v_interp_profile);

            % Poiseuille fit
            warning("off");
            [~, centt] = max(mean_v_interp_profile);
            r_range = (1:numInterp) - centt;
            f = fit(r_range', mean_v_interp_profile(:), 'poly2');
            poiseuille_fit = f.p1 * r_range .^ 2 + f.p2 * r_range + f.p3;
            poiseuille_fit(poiseuille_fit < 0) = 0;
            set(poiseuillePlot, 'XData', 1:numInterp, 'YData', poiseuille_fit);
            warning("on");

            % Capture frame
            video(:, :, :, frameIdx) = rescale(frame2im(getframe(fig)));
        end

        % Write video to disk
        writeGifOnDisc(video, sprintf("interp_profile_%s%d", name, cIdx), "ToolBox", ToolBox);
    end

end

end
