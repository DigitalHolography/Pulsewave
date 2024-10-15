function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, v_RMS, M0_disp_video, ToolBox, k, path, flagBloodVelocityProfile)

    PW_params = Parameters_json(path);
    
    mkdir(ToolBox.PW_path_png, 'bloodVolumeRate')
    mkdir(ToolBox.PW_path_eps, 'bloodVolumeRate')
    
    [numX, numY, numFrames] = size(v_RMS);
    [X, Y] = meshgrid(1:numY, 1:numX);
    
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
    M0_disp_image = rescale(mean(M0_disp_video, 3));
    
    %% All circles testing
    
    % for the all circles output
    numCircles = PW_params.nbCircles;
    maskSectionCircles = cell(1, numCircles);
    delta_rad = (PW_params.velocityBigRadiusRatio - PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 / numCircles; %PW_params.radius_gap
    
    for i = 1:numCircles
        rad_in = (PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 + (i - 1) * delta_rad; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * deltr ;
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
    
    avg_bloodVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
    std_bloodVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
    cross_section_area_artery_r = zeros(numCircles, max(nb_sections_artery), 'single');
    cross_section_mask_artery_r = zeros(numCircles, numY, numX, 'single');
    velocity_profiles_r = cell([numCircles max(nb_sections_artery)]);
    sub_images_r = cell([numCircles max(nb_sections_artery)]);
    force_width = [];
    if ~isempty(PW_params.forcewidth)
        force_width = PW_params.forcewidth;
    end
    for i = 1:numCircles
        [avg_bloodVolumeRate_artery, std_bloodVolumeRate_artery, cross_section_area_artery, ~, cross_section_mask_artery, ~, ~, velocity_profiles, subImg_cell] = crossSectionAnalysis(SubImg_locs_artery_Circles{i}, SubImg_width_artery_Circles{i}, maskArtery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile, i,force_width);
    
        if length(avg_bloodVolumeRate_artery) < 1
            continue
        end
    
        avg_bloodVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(avg_bloodVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
        std_bloodVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(std_bloodVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
        cross_section_area_artery_r(i, 1:nb_sections_artery(i)) = reshape(cross_section_area_artery, 1, nb_sections_artery(i));
        cross_section_mask_artery_r(i, :, :) = reshape(cross_section_mask_artery, 1, numX, numY);
    
        for j = 1:nb_sections_artery(i)
            velocity_profiles_r{i, j} = velocity_profiles{j};
            sub_images_r{i, j} = rescale(subImg_cell{j});
        end
    
    end
    if isempty(PW_params.forcewidth)
        index_start = systoles_indexes(1);
        index_end = systoles_indexes(end);
    else
        index_start = 1;
        index_end = N_frame;
    end
    
    colors = lines(numCircles);
    
    imgRGB = repmat(M0_disp_image, 1, 1, 3);
    
    for i = 1:numCircles
        indxs = find(cross_section_mask_artery_r(i, :, :) > 0);
        imgRGB(indxs) = colors(i, 1);
        imgRGB(numX * numY + indxs) = colors(i, 2);
        imgRGB(2 * numX * numY + indxs) = colors(i, 3);
    
        if i > 1 % intersections should be drawn in white
            indxs = find(cross_section_mask_artery_r(i, :, :) > 0 & cross_section_mask_artery_r(i - 1, :, :) > 0);
            imgRGB(indxs) = 1;
            imgRGB(numX * numY + indxs) = 1;
            imgRGB(2 * numX * numY + indxs) = 1;
        end
    
    end
    
    figure(16774)
    imshow(imgRGB)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'ateries_sections.png')))
    
    figure(11174)
    subimage_size = size(sub_images_r{1,1},1);
    for i=1:numCircles
        for j=1:max(nb_sections_artery)
            if isempty(sub_images_r{i,j})
                sub_images_r{i,j} = zeros(subimage_size,'single');
            end
        end
    end
    montage(sub_images_r(1:numCircles, 1:max(nb_sections_artery)), "Size", [numCircles, max(nb_sections_artery)])
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'all_sections_with_increasing_radius.png')))
    
    figure(16796)
    cross_section_hist = histogram(cross_section_area_artery_r(cross_section_area_artery_r ~= 0), 100);
    title('histogram of cross sections');
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'histogram_of_cross_sections.png')))
    
    %     discarded_cross_sections = ~(2e-3<cross_section_area_artery_r<8e-3);
    %     cross_section_mask_artery_r_new = cross_section_mask_artery_r;
    %     for i =1:nbCircles
    %         [labels,nbsection] = bwlabel(squeeze(cross_section_mask_artery_r_new(i,:,:)));
    %         for j =1:nbsection
    %             if discarded_cross_sections(i,j)
    %                 mask = cross_section_mask_artery_r_new(i,:,:);
    %                 mask (labels==j) = 0;
    %                 cross_section_mask_artery_r_new(i,:,:)=mask;
    %             end
    %         end
    %     end
    %
    %     imgRGB = repmat(mean_M0,1,1,3);
    %     for i =1:nbCircles
    %         indxs = find(cross_section_mask_artery_r_new(i,:,:)>0);
    %         imgRGB(indxs) = colors(i,1);
    %         imgRGB(N*M+indxs) = colors(i,2);
    %         imgRGB(2*N*M+indxs) = colors(i,3);
    %
    %         if i>1 % intersections should be drawn in white
    %             indxs = find(cross_section_mask_artery_r_new(i,:,:)>0&cross_section_mask_artery_r_new(i-1,:,:)>0);
    %             imgRGB(indxs) = 1;
    %             imgRGB(N*M+indxs) = 1;
    %             imgRGB(2*N*M+indxs) = 1;
    %         end
    %     end
    %     figure(164)
    %     imshow(imgRGB)
    %     exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'ateries_sections_after_discard.png')))
    %
    
    plot_bvr_full_field = figure(1676);
    
    Color_std = [0.7 0.7 0.7];
    rad = ((PW_params.velocitySmallRadiusRatio * (numY + numX) / 2) + delta_rad / 2:delta_rad:(PW_params.velocityBigRadiusRatio * (numY + numX) / 2) - delta_rad / 2)'';
    bvr_r = sum(avg_bloodVolumeRateArteryR, 2);
    std_bvr_r = sum(std_bloodVolumeRateArteryR, 2);
    mean_bvr_r = squeeze(mean(bvr_r(:,:,index_start:index_end), 3))';
    mean_std_bvr_r = squeeze(mean(std_bvr_r(:,:,index_start:index_end), 3))';
    curve1 = mean_bvr_r + 0.5 * mean_std_bvr_r;
    curve2 = mean_bvr_r - 0.5 * mean_std_bvr_r;
    rad_out = [rad, fliplr(rad)];
    inBetween = [curve1, fliplr(curve2)]';
    
    fill(rad_out, inBetween, Color_std);
    hold on;
    plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(rad, mean_bvr_r, '-k', 'LineWidth', 2);
    yline(mean(mean_bvr_r),'--k', 'LineWidth', 2);
    legend({'','','','',sprintf('mean = %f µL/min',mean(mean_bvr_r)),'',''});
    
    axis tight;
    aa = axis;
    aa(3) = 0;
    aa(4) = 1.3 * aa(4);
    axis(aa);
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('radius in pixels')
    title("Mean over time Total Blood Volume Rate with the radius to center CRA")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'meanbloodVolumeRatexradius.png')))
    
    plot_bvr_r_variance = figure(1677);
    
    hold on;
    
    for i = 1:numCircles
        plot(fullTime, squeeze(bvr_r(i, :, :)), 'LineWidth', 2);
    end
    
    axis tight;
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Total Blood Volume Rate over time in each artery sections")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRatevariancextime.png')))
    
    plot_bvr_t = figure(1579);
    
    mean_bvr_t = squeeze(mean(bvr_r, 1))';
    mean_bvr_t_value = mean(mean_bvr_t(index_start:index_end));
    mean_std_bvr_t = squeeze(mean(std_bvr_r, 1))';
    
    hold off
    curve1 = mean_bvr_t + 0.5 * mean_std_bvr_t;
    curve2 = mean_bvr_t - 0.5 * mean_std_bvr_t;
    ft2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)]';
    
    fill(ft2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, mean_bvr_t, '-k', 'LineWidth', 2);
    yline(mean_bvr_t_value,'--k', 'LineWidth', 2)
    
    plot(fullTime(index_start), 2.5*mean_bvr_t_value, 'k|', 'MarkerSize', 10); 
    plot(fullTime(index_end), 2.5*mean_bvr_t_value, 'k|', 'MarkerSize', 10); 
    plot(fullTime(index_start:index_end),repmat(2.5*mean_bvr_t_value,index_end-index_start+1),'-k');
    legend({'','','','',sprintf('mean = %f µL/min',mean_bvr_t_value),'',''});
    
    axis tight;
    aa = axis;
    aa(3) = 0;
    aa(4) = 1.3 * aa(4);
    axis(aa);
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Average over all radii Total Blood Volume Rate over time")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateallradxtime.png')))
    
    end
    