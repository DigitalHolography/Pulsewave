function interpolatedBloodVelocityProfile(v_profiles_avg_r, v_profiles_std_r, numSections, rad, numInterp)

ToolBox = getGlobalToolBox;
numCircles = size(v_profiles_avg_r, 2);

for circleIdx = 1:numCircles
    figure("Visible","off");

    for j = 1:numSections(circleIdx)
        plot(mean(v_profiles_avg_r{circleIdx, j}, 2))
        hold on
    end

    colors = lines(numSections(circleIdx));

    for j = 1:numSections(circleIdx)
        profile = mean(v_profiles_avg_r{circleIdx, j}, 2);

        if any(profile < 0) % edge case when there is negative velocities
            [~, locs] = findpeaks(-profile);
            % we find the minimums and set them as the borders of the
            % vessel profile
            if length(locs) > 1
                indx = locs(1):locs(end);
            else

                if locs(1) > length(profile) / 2
                    indx = 1:locs(1);
                else
                    indx = locs(1):length(profile);
                end

            end

        else % main case
            indx = find(profile > 0);
        end

        plot(indx, ones([1 length(indx)]) * mean(mean(v_profiles_avg_r{circleIdx, j}, 2)), 'Color', colors(j, :))
        hold on
    end

    title(['Measured time-averaged velocity profiles at radius = ', num2str(rad(circleIdx)), ' pix'])
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'bloodVelocityProfiles.png')))
    
    figure("Visible","off");
    interp_profile = zeros([numSections(circleIdx), numInterp], 'single');
    interp_profile_std = zeros([numSections(circleIdx), numInterp], 'single');

    for j = 1:numSections(circleIdx)

        profile = mean(v_profiles_avg_r{circleIdx, j}, 2); % mean velocity profile
        profile_std = mean(v_profiles_std_r{circleIdx, j}, 2);

        if any(profile < 0) % edge case when there is negative velocities
            [~, locs] = findpeaks(-profile);
            % we find the minimums and set them as the borders of the
            % vessel profile
            if length(locs) > 1
                indx = locs(1):locs(end);
            else

                if locs(1) > length(profile) / 2
                    indx = 1:locs(1);
                else
                    indx = locs(1):length(profile);
                end

            end

        else % main case
            indx = find(profile > 0);
        end

        interp_profile(j, :) = interp1(1:length(indx), profile(indx), linspace(1, length(indx), numInterp));
        interp_profile_std(j, :) = interp1(1:length(indx), profile_std(indx), linspace(1, length(indx), numInterp));
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

    title(['Interpolated time-averaged velocity profile at radius = ', num2str(rad(circleIdx)), ' pix'])
    set(gca, 'Linewidth', 2)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'interpolatedBloodVelocityProfile.png')))

end
end