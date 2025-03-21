function interpolatedBloodVelocityProfile(v_profiles_cell, dv_profiles_cell, numSections, name, rad, numInterp)

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
main_folder = ToolBox.main_foldername;
params = ToolBox.getParams;
numCircles = size(v_profiles_cell, 2);
Color_err = [0.7 0.7 0.7];

parfor cIdx = 1:numCircles

    % bloodVelocityProfiles Figure
    figure("Visible", "off");

    for sectionIdx = 1:numSections(cIdx)
        plot(v_profiles_cell{cIdx}{sectionIdx}, Linewidth = 2)
        hold on
    end

    title(sprintf('measured time-averaged velocity profiles at radius = %d pix', rad(cIdx)))
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_bloodVelocity_profiles_%s%d.png", main_folder, name, cIdx)))
    % interpolatedBloodVelocityProfile Figure

    figure("Visible", "off");
    v_profile_interp = zeros([numSections(cIdx), numInterp], 'double');
    dv_interp_profile = zeros([numSections(cIdx), numInterp], 'double');

    for sectionIdx = 1:numSections(cIdx)

        v_profile = v_profiles_cell{cIdx}{sectionIdx}; % mean velocity profile
        dv_profile = dv_profiles_cell{cIdx}{sectionIdx};

        if any(v_profile < 0) % edge case when there is negative velocities
            [~, locs] = findpeaks(-v_profile);
            % we find the minimums and set them as the borders of the
            % vessel profile
            if length(locs) > 1
                indx = locs(1):locs(end);
            else

                if isempty(locs)
                    indx = find(v_profile > 0);
                elseif locs(1) > length(v_profile) / 2
                    indx = 1:locs(1);
                else
                    indx = locs(1):length(v_profile);
                end

            end

        else % main case
            indx = find(v_profile > 0);
        end

        L = length(indx);
        v_profile_interp(sectionIdx, :) = interp1(1:L, v_profile(indx), linspace(1, L, numInterp));
        dv_interp_profile(sectionIdx, :) = interp1(1:L, dv_profile(indx), linspace(1, L, numInterp));
    end

    mean_v_interp_profile = sum(v_profile_interp, 1);
    rms_dv_interp_profile = sqrt(sum(dv_interp_profile .^ 2, 1));

    curve1 = mean_v_interp_profile + rms_dv_interp_profile;
    curve2 = mean_v_interp_profile - rms_dv_interp_profile;
    ft2 = [(1:numInterp), fliplr(1:numInterp)];
    inBetween = [curve1, fliplr(curve2)]';

    fill(ft2, inBetween, Color_err);
    hold on;
    plot(1:numInterp, curve1, "Color", Color_err, 'LineWidth', 2);
    plot(1:numInterp, curve2, "Color", Color_err, 'LineWidth', 2);
    plot(1:numInterp, mean_v_interp_profile, '-k', 'LineWidth', 2);
    axis tight;

    % adding a poiseuille fiting (poly2)
    [~, centt] = max(mean_v_interp_profile);
    central_range = 1:numInterp; %max(1,centt-round(Ninterp/6)):min(Ninterp,centt+round(Ninterp/6));
    r_range = (central_range - centt);
    f = fit(r_range', mean_v_interp_profile(central_range)', 'poly2');
    poiseuille_fit = f.p1 * ((1:numInterp) -centt) .^ 2 + f.p2 * ((1:numInterp) -centt) + f.p3;
    poiseuille_fit(poiseuille_fit < 0) = 0;
    plot(poiseuille_fit, '-r', 'LineWidth', 2);

    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), axP(3), 1.07 * axP(4)])
    hold off

    title(['interpolated time-averaged velocity profile at radius = ', num2str(rad(cIdx)), ' pix'])
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_interp_profile_%s%d.png", ToolBox.main_foldername, name, cIdx)))

end

numFrames = size(v_profiles_cell{1}, 2);

if params.exportVideos

    fig = figure("Visible", "off");
    ax = axes(fig);
    hold(ax, "on");

    fillPlot = fill(ax, NaN, NaN, Color_err, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Preallocate fill area
    curve1Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    curve2Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    meanPlot = plot(ax, NaN, NaN, '-k', 'LineWidth', 2);
    poiseuillePlot = plot(ax, NaN, NaN, '-r', 'LineWidth', 2);

    axis tight;
    ax.YLim = [-10; 30];
    box on

    set(gca, 'Linewidth', 2)

    parfor cIdx = 1:numCircles
        video = zeros(420, 560, 3, numFrames, 'single'); % Preallocate video array
        title(['interpolated time-averaged velocity profile at radius = ', num2str(rad(cIdx)), ' pix'])

        % Precompute profiles
        v_profile_interp = zeros([numSections(cIdx), numInterp], 'double');
        dv_interp_profile = zeros([numSections(cIdx), numInterp], 'double');

        for frameIdx = 1:numFrames

            for sectionIdx = 1:numSections(cIdx)
                v_profile = v_profiles_cell{cIdx}{sectionIdx, frameIdx};
                dv_profile = dv_profiles_cell{cIdx}{sectionIdx, frameIdx};

                if any(v_profile < 0) % edge case when there is negative velocities
                    [~, locs] = findpeaks(-v_profile);
                    % we find the minimums and set them as the borders of the
                    % vessel profile
                    if length(locs) > 1
                        indx = locs(1):locs(end);
                    else

                        if isempty(locs)
                            indx = find(v_profile > 0);
                        elseif locs(1) > length(v_profile) / 2
                            indx = 1:locs(1);
                        else
                            indx = locs(1):length(v_profile);
                        end

                    end

                else % main case
                    indx = find(v_profile > 0);
                end

                L = length(indx);
                v_profile_interp(sectionIdx, :) = interp1(1:L, v_profile(indx), linspace(1, L, numInterp));
                dv_interp_profile(sectionIdx, :) = interp1(1:L, dv_profile(indx), linspace(1, L, numInterp));
            end

            % Compute new plot data
            mean_v_interp_profile = mean(v_profile_interp, 1);
            rms_dv_interp_profile = sqrt(sum(dv_interp_profile .^ 2, 1)) / size(dv_interp_profile, 1);
            curve1 = mean_v_interp_profile + rms_dv_interp_profile;
            curve2 = mean_v_interp_profile - rms_dv_interp_profile;
            inBetween = [curve1, fliplr(curve2)]';
            ft2 = [(1:numInterp), fliplr(1:numInterp)];

            % Update plot elements instead of re-creating them
            set(fillPlot, 'XData', ft2, 'YData', inBetween);
            set(curve1Plot, 'XData', 1:numInterp, 'YData', curve1);
            set(curve2Plot, 'XData', 1:numInterp, 'YData', curve2);
            set(meanPlot, 'XData', 1:numInterp, 'YData', mean_v_interp_profile);

            % Poiseuille fit
            warning("off")
            [~, centt] = max(mean_v_interp_profile);
            r_range = (1:numInterp) - centt;
            f = fit(r_range', mean_v_interp_profile(:), 'poly2');
            poiseuille_fit = f.p1 * r_range .^ 2 + f.p2 * r_range + f.p3;
            poiseuille_fit(poiseuille_fit < 0) = 0;
            set(poiseuillePlot, 'XData', 1:numInterp, 'YData', poiseuille_fit);
            warning("on")

            % Capture frame
            video(:, :, :, frameIdx) = rescale(frame2im(getframe(fig)));
        end

        % Write only once per circle
        writeGifOnDisc(video, sprintf("interp_profile_%s%d", name, cIdx), "ToolBox", ToolBox);
    end

end

end
