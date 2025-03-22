function diasysVelocityProfile(v_cell, dv_cell, sysIdx, diasIdx, numSections, name, numInterp)

ToolBox = getGlobalToolBox;
numCircles = size(v_cell, 2);
numFrames = size(v_cell{1}, 2);

% Preallocate interpolated profiles
v_interp = cell(1, numCircles);
dv_interp = cell(1, numCircles);

parfor cIdx = 1:numCircles
    v_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp);
    dv_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp);

    for sectionIdx = 1:numSections(cIdx)
        for frameIdx = 1:numFrames
            v = v_cell{cIdx}{sectionIdx, frameIdx}; % Mean velocity profile
            dv = dv_cell{cIdx}{sectionIdx, frameIdx};

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
            v_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, v(indx), linspace(1, L, numInterp));
            dv_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, dv(indx), linspace(1, L, numInterp));
        end
    end

end

% Average over the cicles and the systolic or diastolic frames

v_sys = zeros(numInterp, 1);
dv_sys = zeros(numInterp, 1);
v_dias = zeros(numInterp, 1);
dv_dias = zeros(numInterp, 1);

for cIdx = 1:numCircles

    for sectionIdx = 1:numSections(cIdx)
        v_sys = v_sys + squeeze(sum(v_interp{cIdx}(sectionIdx, sysIdx, :), 2) / numSections(cIdx) );
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
ft2_sys = [(1:numInterp), fliplr((1:numInterp))];

curve1_dias = v_dias + dv_dias;
curve2_dias = v_dias - dv_dias;
inBetween_dias = [curve1_dias, fliplr(curve2_dias)];
ft2_dias = [(1:numInterp), fliplr((1:numInterp))];

figure("Visible", "off")
hold on
% Plot systolic and diastolic profiles

Color_err_sys = [0.8 0.6 0.6];
ax = gca;
fill(ax, ft2_sys, inBetween_sys, Color_err_sys, 'EdgeColor', 'none');
plot(ax, 1:numInterp, curve1_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, curve2_sys, 'Color', Color_err_sys, 'LineWidth', 2);
plot(ax, 1:numInterp, v_sys, '--r', 'LineWidth', 2);

Color_err_dias = [0.6 0.6 0.8];
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

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)])
hold off

axis tight;
box on
set(gca, 'Linewidth', 2)

xlabel('Profile (AU)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);
pbaspect([1.618 1 1]);

exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_diasys_%s.png", ToolBox.main_foldername, name)))

end
