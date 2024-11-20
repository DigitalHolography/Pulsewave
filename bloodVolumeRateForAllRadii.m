function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, v_RMS, M0_disp_video, ToolBox, k, path, flagBloodVelocityProfile, systolesIndexes)

    PW_params = Parameters_json(path);

    mkdir(ToolBox.PW_path_png, 'volumeRate')
    mkdir(ToolBox.PW_path_eps, 'volumeRate')

    [numX, numY, numFrames] = size(v_RMS);
    [X, Y] = meshgrid(1:numY, 1:numX);

    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
    M0_disp_image = rescale(mean(M0_disp_video, 3));

    %% All circles testing

    % for the all circles output
    numCircles = PW_params.nbCircles;
    maskSectionCircles = cell(1, numCircles);
    delta_rad = (PW_params.velocityBigRadiusRatio - PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 / numCircles; %PW_params.radius_gap
    mask_allSections = createMaskSection(M0_disp_image, maskArtery, (PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2, (PW_params.velocityBigRadiusRatio) * (numY + numX) / 2, sprintf('_mask_artery_all_sections.png'), ToolBox, path);

    for i = 1:numCircles
        rad_in = (PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 + (i - 1) * delta_rad; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * delta_rad ;
        rad_out = rad_in + delta_rad;
        c1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad_in;
        c2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad_out;
        maskSectionCircles(i) = {xor(c1, c2)};

        % save mask image
        createMaskSection(M0_disp_image, maskArtery, rad_in, rad_out, sprintf('_mask_artery_section_circle_%d.png', i), ToolBox, path);
    end

    close(156);

    % for all circles output

    SubImg_locs_artery_Circles = cell(numCircles);
    SubImg_width_artery_Circles = cell(numCircles);
    nb_sections_artery = zeros(1, numCircles);

    for i = 1:numCircles
        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        [maskSectionArtery, n_] = bwlabel(maskSectionArtery);

        nb_sections_artery(i) = n_;
        masksSectionsArtery = zeros(numX, numY, nb_sections_artery(i));

        parfor section_idx = 1:nb_sections_artery(i)
            masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);
        end

        SubImg_locs_artery = zeros(nb_sections_artery(i), 2);
        SubImg_width_artery = zeros(nb_sections_artery(i));

        for section_idx = 1:nb_sections_artery(i)
            [row, col] = find(masksSectionsArtery(:, :, section_idx));
            SubImg_locs_artery(section_idx, 1) = round(mean(row));
            SubImg_locs_artery(section_idx, 2) = round(mean(col));
            SubImg_width_artery(section_idx) = 0.01 * size(maskArtery, 1);
        end

        SubImg_width_artery_Circles{i} = SubImg_width_artery;
        SubImg_locs_artery_Circles{i} = SubImg_locs_artery;
    end

    % for all circles output

    avgVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
    stdVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
    cross_section_area_artery_r = zeros(numCircles, max(nb_sections_artery), 'single');
    cross_section_mask_artery_r = zeros(numCircles, numY, numX, 'single');
    stdCrossSectionWidthR = zeros(numCircles, max(nb_sections_artery), 'single');
    velocity_profiles_r = cell([numCircles max(nb_sections_artery)]);
    std_velocity_profiles_r = cell([numCircles max(nb_sections_artery)]);
    sub_images_r = cell([numCircles max(nb_sections_artery)]);
    force_width = [];

    if ~isempty(PW_params.forcewidth)
        force_width = PW_params.forcewidth;
    end

    for i = 1:numCircles
        [avgVolumeRate_artery, stdVolumeRate_artery, cross_section_area_artery, ~, ~, cross_section_mask_artery, velocity_profiles, std_velocity_profiles, subImg_cell, ~, stdCrossSectionWidth] = crossSectionAnalysis2(SubImg_locs_artery_Circles{i}, SubImg_width_artery_Circles{i}, maskArtery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile, i, force_width, 1);

        if length(avgVolumeRate_artery) < 1
            continue
        end

        avgVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(avgVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
        stdVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(stdVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
        cross_section_area_artery_r(i, 1:nb_sections_artery(i)) = reshape(cross_section_area_artery, 1, nb_sections_artery(i));
        stdCrossSectionWidthR(i, 1:nb_sections_artery(i)) = reshape(stdCrossSectionWidth, 1, nb_sections_artery(i));
        cross_section_mask_artery_r(i, :, :) = reshape(cross_section_mask_artery, 1, numX, numY);

        for j = 1:nb_sections_artery(i)
            velocity_profiles_r{i, j} = velocity_profiles{j};
            std_velocity_profiles_r{i, j} = std_velocity_profiles{j};
            sub_images_r{i, j} = rescale(subImg_cell{j});
        end

    end

    if isempty(PW_params.forcewidth)
        index_start = systolesIndexes(1);
        index_end = systolesIndexes(end);
    else
        index_start = 1;
        index_end = numFrames;
    end

    colors = lines(numCircles);

    imgRGB = repmat(M0_disp_image, 1, 1, 3);

    for i = 1:numCircles
        indxs = find(cross_section_mask_artery_r(i, :, :) > 0);
        imgRGB(indxs) = colors(i, 1);
        imgRGB(numY * numX + indxs) = colors(i, 2);
        imgRGB(2 * numY * numX + indxs) = colors(i, 3);

        if i > 1 % intersections should be drawn in white
            indxs = find(cross_section_mask_artery_r(i, :, :) > 0 & cross_section_mask_artery_r(i - 1, :, :) > 0);
            imgRGB(indxs) = 1;
            imgRGB(numY * numX + indxs) = 1;
            imgRGB(2 * numY * numX + indxs) = 1;
        end

    end

    figure(16774)
    imshow(imgRGB)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'ateries_sections.png')))

    figure(11174)
    % fill with zero images the zeros parts
    subimage_size = size(sub_images_r{1, 1}, 1);

    for i = 1:numCircles

        for j = 1:max(nb_sections_artery)

            if isempty(sub_images_r{i, j})
                sub_images_r{i, j} = zeros(subimage_size, 'single');
            end

        end

    end

    montage(sub_images_r(1:numCircles, 1:max(nb_sections_artery)), "Size", [max(nb_sections_artery), numCircles])
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'all_sections_with_increasing_radius.png')))

    section_width_plot = figure(430);
    mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate'), 'sectionsWidth')
    mkdir(fullfile(ToolBox.PW_path_eps, 'volumeRate'), 'sectionsWidth')
    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for i = 1:numCircles
        section_width_plot.Position = [200 200 600 600];
        crossSectionWidthArtery = 2 * sqrt(cross_section_area_artery_r(i, 1:nb_sections_artery(i)) / pi) * 1000;
        etiquettes_frame_values = append(string(round(crossSectionWidthArtery, 1)), "µm");
        graphMaskTags(section_width_plot, M0_disp_image, squeeze(cross_section_mask_artery_r(i, :, :)), SubImg_locs_artery_Circles{i}, etiquettes_frame_values, x_center, y_center, Fontsize = 12);
        title(sprintf("%s", 'Cross section width in arteries (µm)'));
        set(gca, 'FontSize', 14)

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'crossSectionWidthArteryImage.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'crossSectionWidthArteryImage.eps')))
    end

    figure(16796)
    cross_section_hist = histogram(2 * sqrt(cross_section_area_artery_r(cross_section_area_artery_r ~= 0) / pi) * 1000, 50, FaceColor = 'k');
    aa = axis;
    aa(4) = aa(4) * 1.14;
    axis(aa);
    title('Histogram of sections width (µm)');
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'histogram_of_section_width.png')))
    writematrix(2 * sqrt(cross_section_area_artery_r / pi) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'section_widths.txt')));
    writematrix(stdCrossSectionWidthR * PW_params.cropSection_pixelSize / (2 ^ PW_params.k) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'standard_deviation_section_width.txt')));

    plot_Bvr_full_field = figure(1676);

    Color_std = [0.7 0.7 0.7];
    rad = ((PW_params.velocitySmallRadiusRatio * (numX + numY) / 2) + delta_rad / 2:delta_rad:(PW_params.velocityBigRadiusRatio * (numX + numY) / 2) - delta_rad / 2)'';
    BvrR = sum(avgVolumeRateArteryR, 2);
    std_BvrR = sqrt(sum(stdVolumeRateArteryR .^ 2, 2)); % sqrt of the sum of variances
    mean_BvrR = squeeze(mean(BvrR(:, :, index_start:index_end), 3))';
    mean_std_BvrR = squeeze(rms(std_BvrR(:, :, index_start:index_end), 3))'; % quadratic mean
    curve1 = mean_BvrR + 0.5 * mean_std_BvrR;
    curve2 = mean_BvrR - 0.5 * mean_std_BvrR;
    rad2 = [rad, fliplr(rad)];
    inBetween = [curve1, fliplr(curve2)]';

    fill(rad2, inBetween, Color_std);
    hold on;
    plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(rad, mean_BvrR, '-k', 'LineWidth', 2);
    yline(mean(mean_BvrR), '--k', 'LineWidth', 2);
    legend({'', '', '', '', sprintf('mean = %0.2f µL/min', mean(mean_BvrR)), '', ''});

    axis tight;
    aa = axis;
    aa(3) = -5;
    aa(4) = 95;
    axis(aa);
    hold off

    ylabel('Blood Volume Rate (µL/min)')
    xlabel('radius in pixels')
    title("Time average of Blood Volume Rate")
    set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
    set(gca, 'LineWidth', 2)

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'meanvolumeRatexradius.png')))

    plot_BvrR_variance = figure(1677);

    hold on;

    for i = 1:numCircles
        plot(fullTime, squeeze(BvrR(i, :, :)), 'LineWidth', 2);
    end

    axis tight;
    aa = axis;
    aa(3) = -10;
    aa(4) = 95;
    axis(aa);
    hold off
    box on

    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Radial variations of Blood Volume Rate")
    set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
    set(gca, 'Linewidth', 2)

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRatevariancextime.png')))

    plot_BvrT = figure(1579);

    mean_BvrT = squeeze(mean(BvrR, 1))';
    mean_BvrT_value = mean(mean_BvrT(index_start:index_end));
    max_BvrT_value = max(mean_BvrT(index_start:index_end));
    mean_std_BvrT = squeeze(rms(std_BvrR, 1))';

    hold off
    curve1 = mean_BvrT + 0.5 * mean_std_BvrT;
    curve2 = mean_BvrT - 0.5 * mean_std_BvrT;
    ft2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)]';

    fill(ft2, inBetween, Color_std);
    hold on;
    yline(0, 'k-', 'LineWidth', 2)
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, mean_BvrT, '-k', 'LineWidth', 2);
    yline(mean_BvrT_value, '--k', 'LineWidth', 2)

    plot(fullTime(index_start), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
    plot(fullTime(index_end), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
    plot(fullTime(index_start:index_end), repmat(1.07 * max_BvrT_value, index_end - index_start + 1), '-k');
    legend({'', '', '', '', '', sprintf('mean = %0.2f µL/min', mean_BvrT_value), '', ''});

    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), axP(3), axP(4) * 1.2])
    box on

    hold off

    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Radial average of Blood Volume Rate")
    set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
    set(gca, 'Linewidth', 2)

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateallradxtime.png')))

    figure(4350)
    %maskNeigbors = mat2gray(mean(imread(fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png'))), 3)) > 0; % import mask neigbors
    graphCombined(M0_disp_video, imdilate(maskArtery, strel('disk', PW_params.local_background_width)) & mask_allSections, [], [], mean_BvrT, mean_std_BvrT, ToolBox, path, 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in arteries Full Field', 'µL/min', skip = ~PW_params.exportVideos);

    if flagBloodVelocityProfile
        mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles'));

        for i = 1:numCircles
            plot_mean_velocity_profiles = figure(7579 + i);

            for j = 1:nb_sections_artery(i)
                plot(mean(velocity_profiles_r{i, j}, 2))
                hold on
            end

            colors = lines(nb_sections_artery(i));

            for j = 1:nb_sections_artery(i)
                profile = mean(velocity_profiles_r{i, j}, 2);

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

                plot(indx, ones([1 length(indx)]) * mean(mean(velocity_profiles_r{i, j}, 2)), 'Color', colors(j, :))
                hold on
            end

            title(['Measured time-averaged velocity profiles at radius = ', num2str(rad(i)), ' pix'])
            set(gca, 'Linewidth', 2)
            exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'bloodVelocityProfiles.png')))

            title(['Mean velocity profiles at radius = ', num2str(rad(i)), ' pix'])
            plot_inter_velocity_profile = figure(7503 + i);
            Ninterp = 50;
            interp_profile = zeros([nb_sections_artery(i), Ninterp], 'single');
            interp_profile_std = zeros([nb_sections_artery(i), Ninterp], 'single');

            for j = 1:nb_sections_artery(i)

                profile = mean(velocity_profiles_r{i, j}, 2); % mean velocity profile
                profile_std = mean(std_velocity_profiles_r{i, j}, 2);

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

                interp_profile(j, :) = interp1(1:length(indx), profile(indx), linspace(1, length(indx), Ninterp));
                interp_profile_std(j, :) = interp1(1:length(indx), profile_std(indx), linspace(1, length(indx), Ninterp));
            end

            mean_interp_profile = mean(interp_profile, 1);
            std_interp_profile = mean(interp_profile_std, 1);
            curve1 = mean_interp_profile + 0.5 * std_interp_profile;
            curve2 = mean_interp_profile - 0.5 * std_interp_profile;
            ft2 = [(1:Ninterp), fliplr(1:Ninterp)];
            inBetween = [curve1, fliplr(curve2)]';

            fill(ft2, inBetween, Color_std);
            hold on;
            plot(1:Ninterp, curve1, "Color", Color_std, 'LineWidth', 2);
            plot(1:Ninterp, curve2, "Color", Color_std, 'LineWidth', 2);
            plot(1:Ninterp, mean_interp_profile, '-k', 'LineWidth', 2);
            axis tight;

            % adding a poiseuille fiting (poly2)
            [~, centt] = max(mean_interp_profile);
            central_range = 1:Ninterp; %max(1,centt-round(Ninterp/6)):min(Ninterp,centt+round(Ninterp/6));
            r_range = (central_range - centt);
            f = fit(r_range', mean_interp_profile(central_range)', 'poly2');
            poiseuille_fit = f.p1 * ((1:Ninterp) -centt) .^ 2 + f.p2 * ((1:Ninterp) -centt) + f.p3;
            poiseuille_fit(poiseuille_fit < 0) = 0;
            plot(poiseuille_fit, '-r', 'LineWidth', 2);

            axis tight;
            aa = axis;
            aa(3) = -10;
            aa(4) = 30;
            axis(aa);
            hold off

            title(['Interpolated time-averaged velocity profile at radius = ', num2str(rad(i)), ' pix'])
            set(gca, 'Linewidth', 2)
            exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'interpolatedBloodVelocityProfile.png')))

        end

    end

    plot_interp_pulse = figure(7124);
    Ninterp = 1000;

    [interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(mean_BvrT, systolesIndexes, Ninterp, mean_std_BvrT);
    dt = (fullTime(2) - fullTime(1));
    pulseTime = dt * (1:Ninterp) * avgLength / Ninterp;

    [~, amin] = min(interp_BvrT);
    [~, amax] = max(interp_BvrT);
    cshiftn = Ninterp - amin;

    hold off

    % Retinal Stroke Volume
    hold on
    curve1 = circshift(interp_BvrT, cshiftn);
    curve2 = 0 * ones(size(curve1));
    ft2 = [pulseTime, fliplr(pulseTime)];
    inBetween = [curve1, fliplr(curve2)]';
    cRose = [254, 191, 210] / 255;
    fill(ft2, inBetween, cRose, 'EdgeColor', 'none');
    xline(pulseTime(end), 'k--', 'LineWidth', 2)

    % Remaining Stroke Volume
    hold on
    curve1 = circshift(interp_BvrT, cshiftn);
    curve1 = curve1(1:amax + cshiftn);
    curve2 = 0 * ones(size(curve1));
    ft2 = [pulseTime(1:amax + cshiftn), fliplr(pulseTime(1:amax + cshiftn))];
    inBetween = [curve1, fliplr(curve2)]';
    cCrimson = [222, 49, 99] / 255;
    xline(pulseTime(amax + cshiftn), 'k--', 'LineWidth', 2)
    fill(ft2, inBetween, cCrimson, 'EdgeColor', 'none');

    % Grey STD and Signal
    interp_BvrT2 = repmat(interp_BvrT, 1, 3);
    interp_std_BvrT2 = repmat(interp_std_BvrT, 1, 3);
    pulseTime2 = dt * (-Ninterp + 1:Ninterp * 2) * avgLength / Ninterp;

    hold on
    curve1 = circshift(interp_BvrT2, cshiftn) + 0.5 * circshift(interp_std_BvrT2, cshiftn);
    curve2 = circshift(interp_BvrT2, cshiftn) - 0.5 * circshift(interp_std_BvrT2, cshiftn);
    ft2 = [pulseTime2, fliplr(pulseTime2)];
    inBetween = [curve1, fliplr(curve2)]';
    cSTD = [0.7 0.7 0.7];
    fill(ft2, inBetween, cSTD, 'EdgeColor', 'none');
    plot(pulseTime2, circshift(interp_BvrT2, cshiftn), '-k', 'LineWidth', 2);
    %
    % axis padded
    % axP = axis;
    % axis([pulseTime(1)-1/2 * pulseTime(end), 3/2 * pulseTime(end), 0, axP(4)])

    yline(0, 'k--', 'LineWidth', 2)
    xline(0, 'k--', 'LineWidth', 2)
    xline(0, 'k--', 'LineWidth', 2)

    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), axP(3), axP(4)])
    xlim([pulseTime(1) - 1/2 * pulseTime(end), 3/2 * pulseTime(end)])
    box on

    ylabel('Blood Volume Rate (µL/min)')
    xlabel('Time (s)')
    ccinterpBvrT = circshift(interp_BvrT, cshiftn);
    dt2 = pulseTime2(2) - pulseTime2(1);
    stroke_volume_value = sum(ccinterpBvrT(1:amax + cshiftn)) * dt2 / 60 * 1000; % in nL
    total_volume_value = sum(ccinterpBvrT) * dt2 / 60 * 1000;
    title(sprintf("Retinal Stroke Volume : %02.0f nL and Total Volume : %02.0f nL", stroke_volume_value, total_volume_value));
    set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
    box on
    set(gca, 'LineWidth', 2)
    box ('on', Clipping = 'off')

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'strokeAndTotalVolume.png')))

    close all
    fprintf("- Blood Volume Rate for all radii took : %ds\n", round(toc))

end
