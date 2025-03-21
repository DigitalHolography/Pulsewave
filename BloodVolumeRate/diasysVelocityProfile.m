function diasysVelocityProfile(v_cell, dv_cell, sysIdx, diasIdx, numSections, name, numInterp)

return

ToolBox = getGlobalToolBox;
numCircles = size(v_cell, 2);

% Preallocate interpolated profiles
v_interp = cell(1, numCircles);
dv_interp = cell(1, numCircles);
v_sys = zeros(numInterp, 1, 'double');
dv_sys = zeros(numInterp, 1, 'double');
v_dias = zeros(numInterp, 1, 'double');
dv_dias = zeros(numInterp, 1, 'double');

parfor cIdx = 1:numCircles
    v_interp{cIdx} = zeros(numSections(cIdx), numInterp, 'double');
    dv_interp{cIdx} = zeros(numSections(cIdx), numInterp, 'double');

    for sectionIdx = 1:numSections(cIdx)
        v = v_cell{cIdx}{sectionIdx}; % Mean velocity profile
        dv = dv_cell{cIdx}{sectionIdx};

        % Handle edge case with negative velocities
        if any(v < 0)
            [~, locs] = findpeaks(-v);

            if isempty(locs)
                indx = find(v > 0);
            elseif length(locs) > 1
                indx = locs(1):locs(end);
            else

                if locs(1) > length(v) / 2
                    indx = 1:locs(1);
                else
                    indx = locs(1):length(v);
                end

            end

        else
            indx = find(v > 0);
        end

        % Interpolate profiles
        L = length(indx);
        v_interp{cIdx}(sectionIdx, :) = interp1(1:L, v(indx), linspace(1, L, numInterp));
        dv_interp{cIdx}(sectionIdx, :) = interp1(1:L, dv(indx), linspace(1, L, numInterp));
    end

end

% Average over the cicles and the systolic or diastolic frames
for cIdx = 1:numCircles

    for sectionIdx = 1:numSections(cIdx)
        v_sys = v_sys + sum(v_interp{cIdx}(sectionIdx, sysIdx), 2);
        dv_sys = dv_sys + sum(dv_interp{cIdx}(sectionIdx, sysIdx) .^ 2, 2);
        v_dias = v_dias + sum(v_interp{cIdx}(sectionIdx, diasIdx), 2);
        dv_dias = dv_dias + sum(dv_interp{cIdx}(sectionIdx, diasIdx) .^ 2, 2);
    end

end

numSys = length(sysIdx);
numDias = length(diasIdx);
v_sys = v_sys / (numSys * numCircles);
dv_sys = sqrt(dv_sys) / (numSys * numCircles);
v_dias = v_dias / (numDias * numCircles);
dv_dias = sqrt(dv_dias) / (numDias * numCircles);

% Create curves for plotting
curve1_sys = v_sys + dv_sys;
curve2_sys = v_sys - dv_sys;
inBetween_sys = [curve1_sys; flipud(curve2_sys)];
ft2_sys = [(1:numInterp)'; flipud((1:numInterp)')];

curve1_dias = v_dias + dv_dias;
curve2_dias = v_dias - dv_dias;
inBetween_dias = [curve1_dias; flipud(curve2_dias)];
ft2_dias = [(1:numInterp)'; flipud((1:numInterp)')];

figure("Visible", "off")

% Plot systolic and diastolic profiles

Color_err_sys = [0.8 0.6 0.6];
ax = gca;
fill(ax, ft2_sys, inBetween_sys, Color_err_sys, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(ax, 1:numInterp, curve1_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, curve2_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, v_sys, '--r', 'LineWidth', 2);


Color_err_dias = [0.6 0.6 0.8];
fill(ax, ft2_dias, inBetween_dias, Color_err_dias, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(ax, 1:numInterp, curve1_dias, 'Color', Color_err_dias, 'LineWidth', 2);
plot(ax, 1:numInterp, curve2_dias, 'Color', Color_err_dias, 'LineWidth', 2);
plot(ax, 1:numInterp, v_dias, '--b', 'LineWidth', 2);

% Poiseuille fit
warning("off");
[~, centt_sys] = max(v_sys);
r_range_sys = (1:numInterp)' - centt_sys;
f_sys = fit(r_range_sys, v_sys, 'poly2');
poiseuille_fit_sys = f_sys.p1 * r_range_sys .^ 2 + f_sys.p2 * r_range_sys + f_sys.p3;
poiseuille_fit_sys(poiseuille_fit_sys < 0) = 0;

[~, centt_dias] = max(v_dias);
r_range_dias = (1:numInterp)' - centt_dias;
f_dias = fit(r_range_dias, v_dias, 'poly2');
poiseuille_fit_dias = f_dias.p1 * r_range_dias .^ 2 + f_dias.p2 * r_range_dias + f_dias.p3;
poiseuille_fit_dias(poiseuille_fit_dias < 0) = 0;

plot(ax, 1:numInterp, poiseuille_fit_sys, '-r', 'LineWidth', 2);
plot(ax, 1:numInterp, poiseuille_fit_dias, '-b', 'LineWidth', 2);
warning("on");

exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_diasys_%s.png", ToolBox.main_foldername, name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, 'volumeRate', 'velocityProfiles', sprintf("%s_diasys_%s.eps", ToolBox.main_foldername, name)))

end
