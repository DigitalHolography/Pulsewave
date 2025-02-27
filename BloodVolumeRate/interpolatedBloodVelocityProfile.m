function interpolatedBloodVelocityProfile(v_profiles_avg_r, v_profiles_std_r, numSections, name, rad, numInterp)

ToolBox = getGlobalToolBox;
PW_params = ToolBox.getParams;
numCircles = size(v_profiles_avg_r, 2);
Color_std = [0.7 0.7 0.7];

for circleIdx = 1:numCircles

    % bloodVelocityProfiles Figure

    figure("Visible", "off");

    for sectionIdx = 1:numSections(circleIdx)
        plot(mean(v_profiles_avg_r{circleIdx}{sectionIdx}, 2), Linewidth = 2)
        hold on
    end

    title(['measured time-averaged velocity profiles at radius = ', num2str(rad(circleIdx)), ' pix'])
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s_bloodVelocityProfiles.png", ToolBox.main_foldername, circleIdx, name)))

    % interpolatedBloodVelocityProfile Figure

    figure("Visible", "off");
    interp_profile = zeros([numSections(circleIdx), numInterp], 'single');
    interp_profile_std = zeros([numSections(circleIdx), numInterp], 'single');

    parfor sectionIdx = 1:numSections(circleIdx)

        profile_avg = mean(v_profiles_avg_r{circleIdx}{sectionIdx}, 2); % mean velocity profile
        profile_std = mean(v_profiles_std_r{circleIdx}{sectionIdx}, 2);

        if any(profile_avg < 0) % edge case when there is negative velocities
            [~, locs] = findpeaks(-profile_avg);
            % we find the minimums and set them as the borders of the
            % vessel profile
            if length(locs) > 1
                indx = locs(1):locs(end);
            else

                if isempty(locs)
                    indx = find(profile_avg > 0);
                elseif locs(1) > length(profile_avg) / 2
                    indx = 1:locs(1);
                else
                    indx = locs(1):length(profile_avg);
                end

            end

        else % main case
            indx = find(profile_avg > 0);
        end

        interp_profile(sectionIdx, :) = interp1(1:length(indx), profile_avg(indx), linspace(1, length(indx), numInterp));
        interp_profile_std(sectionIdx, :) = interp1(1:length(indx), profile_std(indx), linspace(1, length(indx), numInterp));
    end

    mean_interp_profile = mean(interp_profile, 1);
    std_interp_profile = mean(interp_profile_std, 1);
    curve1 = mean_interp_profile + 0.5 * std_interp_profile;
    curve2 = mean_interp_profile - 0.5 * std_interp_profile;
    ft2 = [(1:numInterp), fliplr(1:numInterp)];
    inBetween = [curve1, fliplr(curve2)]';

    fill(ft2, inBetween, Color_std);
    hold on;
    plot(1:numInterp, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(1:numInterp, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(1:numInterp, mean_interp_profile, '-k', 'LineWidth', 2);
    axis tight;

    % adding a poiseuille fiting (poly2)
    [~, centt] = max(mean_interp_profile);
    central_range = 1:numInterp; %max(1,centt-round(Ninterp/6)):min(Ninterp,centt+round(Ninterp/6));
    r_range = (central_range - centt);
    f = fit(r_range', mean_interp_profile(central_range)', 'poly2');
    poiseuille_fit = f.p1 * ((1:numInterp) -centt) .^ 2 + f.p2 * ((1:numInterp) -centt) + f.p3;
    poiseuille_fit(poiseuille_fit < 0) = 0;
    plot(poiseuille_fit, '-r', 'LineWidth', 2);

    axis tight;
    aa = axis;
    aa(3) = -10;
    aa(4) = 30;
    axis(aa);
    hold off

    title(['interpolated time-averaged velocity profile at radius = ', num2str(rad(circleIdx)), ' pix'])
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s_interpolatedBloodVelocityProfile.png", ToolBox.main_foldername, circleIdx, name)))

end

numFrames = size(v_profiles_avg_r{circleIdx}{1}, 2);

if PW_params.exportVideos

    fig = figure("Visible", "off");
    ax = axes(fig);
    hold(ax, "on");

    fillPlot = fill(ax, NaN, NaN, Color_std, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Preallocate fill area
    curve1Plot = plot(ax, NaN, NaN, "Color", Color_std, 'LineWidth', 2);
    curve2Plot = plot(ax, NaN, NaN, "Color", Color_std, 'LineWidth', 2);
    meanPlot = plot(ax, NaN, NaN, '-k', 'LineWidth', 2);
    poiseuillePlot = plot(ax, NaN, NaN, '-r', 'LineWidth', 2);

    axis tight;
    ax.YLim = [-10; 30];
    box on

    set(gca, 'Linewidth', 2)

    parfor circleIdx = 1:numCircles
        video = zeros(420, 560, 3, numFrames, 'single'); % Preallocate video array
        title(['interpolated time-averaged velocity profile at radius = ', num2str(rad(circleIdx)), ' pix'])

        for frameIdx = 1:numFrames
            % Precompute profiles
            interp_profile = zeros([numSections(circleIdx), numInterp], 'single');
            interp_profile_std = zeros([numSections(circleIdx), numInterp], 'single');

            for sectionIdx = 1:numSections(circleIdx)
                profile_avg = v_profiles_avg_r{circleIdx}{sectionIdx}(:, frameIdx);
                profile_std = v_profiles_std_r{circleIdx}{sectionIdx}(:, frameIdx);

                if any(profile_avg < 0) % edge case when there is negative velocities
                    [~, locs] = findpeaks(-profile_avg);
                    % we find the minimums and set them as the borders of the
                    % vessel profile
                    if length(locs) > 1
                        indx = locs(1):locs(end);
                    else

                        if isempty(locs)
                            indx = find(profile_avg > 0);
                        elseif locs(1) > length(profile_avg) / 2
                            indx = 1:locs(1);
                        else
                            indx = locs(1):length(profile_avg);
                        end

                    end

                else % main case
                    indx = find(profile_avg > 0);
                end

                interp_profile(sectionIdx, :) = interp1(1:length(indx), profile_avg(indx), linspace(1, length(indx), numInterp));
                interp_profile_std(sectionIdx, :) = interp1(1:length(indx), profile_std(indx), linspace(1, length(indx), numInterp));
            end

            % Compute new plot data
            mean_interp_profile = mean(interp_profile, 1);
            std_interp_profile = mean(interp_profile_std, 1);
            curve1 = mean_interp_profile + 0.5 * std_interp_profile;
            curve2 = mean_interp_profile - 0.5 * std_interp_profile;
            inBetween = [curve1, fliplr(curve2)]';
            ft2 = [(1:numInterp), fliplr(1:numInterp)];

            % Update plot elements instead of re-creating them
            set(fillPlot, 'XData', ft2, 'YData', inBetween);
            set(curve1Plot, 'XData', 1:numInterp, 'YData', curve1);
            set(curve2Plot, 'XData', 1:numInterp, 'YData', curve2);
            set(meanPlot, 'XData', 1:numInterp, 'YData', mean_interp_profile);

            % Poiseuille fit
            warning("off")
            [~, centt] = max(mean_interp_profile);
            r_range = (1:numInterp) - centt;
            f = fit(r_range', mean_interp_profile(:), 'poly2');
            poiseuille_fit = f.p1 * r_range .^ 2 + f.p2 * r_range + f.p3;
            poiseuille_fit(poiseuille_fit < 0) = 0;
            set(poiseuillePlot, 'XData', 1:numInterp, 'YData', poiseuille_fit);
            warning("on")

            % Capture frame
            video(:, :, :, frameIdx) = rescale(frame2im(getframe(fig)));
        end

        % Write only once per circle
        writeGifOnDisc(video, sprintf("circle_%d_%s_interpolatedBloodVelocityProfile", circleIdx, name), "ToolBox", ToolBox);
    end

end

end
