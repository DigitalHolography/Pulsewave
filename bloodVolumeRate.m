function [] = bloodVolumeRate(maskArtery, maskVein, vRMS, flatfieldM0, ToolBox, k, path, flagBloodVelocityProfile)

    PW_params = Parameters_json(path);

    veins_analysis = PW_params.veins_analysis;

    mkdir(ToolBox.PW_path_png, 'bloodVolumeRate')
    mkdir(ToolBox.PW_path_eps, 'bloodVolumeRate')

    [numX, numY, numFrames] = size(vRMS);
    [x, y] = meshgrid(1:numY, 1:numX);

    %maskArtery = imdilate(maskArtery, strel('disk', 5));
    %FIXME function velocity map

    v_RMS_AVG = mean(vRMS, 3);
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

    %% change mask section ?

    % radius_ratio = 0.18;
    % radius_gap = 0.05;
    % radius1 = (radius_ratio-radius_gap)* (M+N)/2;
    % radius2 = (radius_ratio+radius_gap)* (M+N)/2;

    radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (numY + numX) / 2;
    radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (numY + numX) / 2;
    cercle_mask1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius1;
    cercle_mask2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius2;

    maskSection = xor(cercle_mask1, cercle_mask2);

    

    %% Find the locations of the sections in arteries

    maskSectionArtery = maskSection .* maskArtery;

    figure(111)
    imagesc(maskSectionArtery .* v_RMS_AVG);

    maskSectionArtery = bwlabel(maskSectionArtery);

    figure(222)
    imagesc(maskSectionArtery)

    nb_sections_artery = max(maskSectionArtery, [], 'all');
    masksSectionsArtery = zeros(numX, numY, nb_sections_artery);

    parfor section_idx = 1:nb_sections_artery

        masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);

    end

    SubImg_locs_artery = zeros(nb_sections_artery, 2);
    SubImg_width_artery = zeros(nb_sections_artery, 1);

    for section_idx = 1:nb_sections_artery
        [row, col] = find(masksSectionsArtery(:, :, section_idx));
        SubImg_locs_artery(section_idx, 1) = round(mean(row));
        SubImg_locs_artery(section_idx, 2) = round(mean(col));
        SubImg_width_artery(section_idx) = 0.01 * numX;
    end

    

    %% Find the locations of the sections in veins

    if veins_analysis

        maskSectionVein = maskSection .* maskVein;

        figure(1111)
        imagesc(maskSectionVein .* v_RMS_AVG);

        maskSectionVein = bwlabel(maskSectionVein);

        figure(2222)
        imagesc(maskSectionVein)

        nb_sections_vein = max(maskSectionVein, [], 'all');
        masksSectionsVein = zeros(numX, numY, nb_sections_vein);

        parfor section_idx = 1:nb_sections_vein
            masksSectionsVein(:, :, section_idx) = (maskSectionVein == section_idx);
        end

        SubImg_locs_vein = zeros(nb_sections_vein, 2);
        SubImg_width_vein = zeros(nb_sections_vein, 1);

        for section_idx = 1:nb_sections_vein
            [row, col] = find(masksSectionsVein(:, :, section_idx));
            SubImg_locs_vein(section_idx, 1) = round(mean(row));
            SubImg_locs_vein(section_idx, 2) = round(mean(col));
            SubImg_width_vein(section_idx) = 0.01 * numX;
        end

    end

    %% Compute blood volume rate

    mask_artery = maskArtery;
    mask_vein = maskVein;

    if veins_analysis
        [avgBloodVolumeRateVein, stdBloodVolumeRateVein, crossSectionAreaVein, avgBloodVelocityVein, crossSectionMaskVein, totalAvgBloodVolumeRateVein, totalStdBloodVolumeRateVein] = cross_section_analysis(SubImg_locs_vein, SubImg_width_vein, mask_vein, vRMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'vein', flagBloodVelocityProfile,[]);
    end

    [avgBloodVolumeRateArtery, stdBloodVolumeRateArtery, crossSectionAreaArtery, avgBloodVelocityArtery, crossSectionMaskArtery, totalAvgBloodVolumeRateArtery, totalStdBloodVolumeRateArtery] = cross_section_analysis(SubImg_locs_artery, SubImg_width_artery, mask_artery, vRMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile,[]);

    labels_arteries = cell(nb_sections_artery, 1);
    strXlabel = 'Time(s)'; %createXlabelTime(1);
    strYlabel = 'Velocity (mm.s-1)';

    data_to_plot_artery = zeros(size(avgBloodVelocityArtery, 1), size(avgBloodVolumeRateArtery, 2)); % avg_blood_velocity (section, time)

    parfor section_idx = 1:nb_sections_artery
        data_to_plot_artery(section_idx, :) = smoothdata(avgBloodVelocityArtery(section_idx, :), 'lowess');
        labels_arteries{section_idx} = strcat('A', num2str(section_idx));
    end

    figure(57)
    plot(fullTime, data_to_plot_artery, 'LineWidth', 2)
    title('Blood velocity in artery sections');
    legend(labels_arteries);
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInArterySection.eps')))

    data_to_plot_artery_all = squeeze(mean(data_to_plot_artery, 1));
    plot_artery_mean = mean(data_to_plot_artery_all);

    figure(571)
    plot(fullTime, data_to_plot_artery_all, 'k', 'LineWidth', 2)
    yline(plot_artery_mean, '--k', 'LineWidth', 2)
    title('Blood velocity in artery sections');
    legend(sprintf('All Artery Sections'));
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInAllArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInAllArterySection.eps')))

    if veins_analysis

        labels_veins = cell(nb_sections_vein, 1);

        data_to_plot_vein = zeros(size(avgBloodVelocityVein, 1), size(avgBloodVolumeRateVein, 2));

        parfor section_idx = 1:nb_sections_vein
            data_to_plot_vein(section_idx, :) = smoothdata(avgBloodVelocityVein(section_idx, :), 'lowess');
            labels_veins{section_idx} = strcat('A', num2str(section_idx));
        end

        figure(58)
        plot(fullTime, data_to_plot_vein, 'LineWidth', 2)
        title('Blood velocity in vein sections');
        legend(labels_veins);
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInVeinSection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInVeinSection.eps')))

        data_to_plot_vein_all = squeeze(mean(data_to_plot_vein, 1));
        plot_vein_mean = mean(data_to_plot_vein_all);

        figure(581)
        plot(fullTime, data_to_plot_vein_all, 'k', 'LineWidth', 2)
        yline(plot_vein_mean, '--k', 'LineWidth', 2)
        title('Blood velocity in vein sections');
        legend(sprintf('All Vein Sections'));
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInAllVeinSection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityInAllVeinSection.eps')))
    end

    maskOnes = ones(numX, numY, numFrames);
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
    referenceMean = rescale(mean(flatfieldM0, 3));
    ratio_etiquette = 1.2;
    %% Saving Vein and artery data in txt

    for section_idx = 1:nb_sections_artery
        plot2txt(fullTime, avgBloodVolumeRateArtery(section_idx, :), strcat('bloodVolumeRate_artery_A', num2str(section_idx)), ToolBox)
        plot2txt(fullTime, stdBloodVolumeRateArtery(section_idx, :), strcat('bloodVolumeRate_artery_std_A', num2str(section_idx)), ToolBox)
        plot2txt(fullTime, avgBloodVelocityArtery(section_idx, :), strcat('avg_velocity_artery_A', num2str(section_idx)), ToolBox)

    end

    if veins_analysis

        for section_idx = 1:nb_sections_vein
            plot2txt(fullTime, avgBloodVolumeRateVein(section_idx, :), strcat('bloodVolumeRate_vein_V', num2str(section_idx)), ToolBox)
            plot2txt(fullTime, stdBloodVolumeRateVein(section_idx, :), strcat('bloodVolumeRate_vein_std_V', num2str(section_idx)), ToolBox)
            plot2txt(fullTime, avgBloodVelocityVein(section_idx, :), strcat('avg_velocity_vein_V', num2str(section_idx)), ToolBox)

        end

    end

    %% Vein and artery numerotation

    maskRGB_artery(:, :, 3) = referenceMean .* ~crossSectionMaskArtery;
    maskRGB_artery(:, :, 2) = referenceMean .* ~crossSectionMaskArtery;
    maskRGB_artery(:, :, 1) = referenceMean .* ~crossSectionMaskArtery + crossSectionMaskArtery;

    figure(652)
    imshow(maskRGB_artery);

    if veins_analysis
        maskRGB_combined(:, :, 3) = (referenceMean .* ~crossSectionMaskVein + crossSectionMaskVein) .* ~crossSectionMaskArtery;
        maskRGB_combined(:, :, 2) = referenceMean .* ~(crossSectionMaskArtery + crossSectionMaskVein);
        maskRGB_combined(:, :, 1) = (referenceMean .* ~crossSectionMaskArtery + crossSectionMaskArtery) .* ~crossSectionMaskVein;

        maskRGB_vein(:, :, 1) = referenceMean .* ~crossSectionMaskVein;
        maskRGB_vein(:, :, 2) = referenceMean .* ~crossSectionMaskVein;
        maskRGB_vein(:, :, 3) = referenceMean .* ~crossSectionMaskVein + crossSectionMaskVein;

        figure(653)
        imshow(maskRGB_vein);
        figure(654)
        imshow(maskRGB_combined);
    end

    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for section_idx = 1:nb_sections_artery
        new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);
        figure(652)
        text(new_x, new_y, strcat("A", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        if veins_analysis
            figure(654)
            text(new_x, new_y, strcat("A", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end
    end

    if veins_analysis

        for section_idx = 1:nb_sections_vein
            new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);
                    figure(653)
            text(new_x, new_y, strcat("V", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                    figure(654)
            text(new_x, new_y, strcat("V", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end

    end

    drawnow
    ax = gca;
    ax.Units = 'pixels';

    figure(652)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s",ToolBox.main_foldername, 'arteries_numerotation.png')))
    if veins_analysis
        figure(653)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s",ToolBox.main_foldername, 'veins_numerotation.png')))
        figure(654)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s",ToolBox.main_foldername, 'arteries_veins_numerotation.png')))
    end

    %% Volume Rate FigurePosition

    Color_std = [0.7 0.7 0.7];
    figure_width = 600;
    figure_height = figure_width;

    volume_rate_video_artery = zeros(figure_width, figure_height, 3, numFrames);
    mean_volume_rate_artery = mean(totalAvgBloodVolumeRateArtery);

    maskRGB_artery(:, :, 3) = referenceMean .* ~crossSectionMaskArtery;
    maskRGB_artery(:, :, 2) = referenceMean .* ~crossSectionMaskArtery;
    maskRGB_artery(:, :, 1) = referenceMean .* ~crossSectionMaskArtery + crossSectionMaskArtery;

    f100 = figure(100);
    colormap("gray")
    f100.Position = [200, 200, 600, 600];
    ax100 = gca;
    ax100.Units = 'pixels';

    if veins_analysis

        volume_rate_video_vein = zeros(figure_width, figure_height, 3, numFrames);
        mean_volume_rate_vein = mean(totalAvgBloodVolumeRateVein);

        maskRGB_vein(:, :, 1) = referenceMean .* ~crossSectionMaskVein;
        maskRGB_vein(:, :, 2) = referenceMean .* ~crossSectionMaskVein;
        maskRGB_vein(:, :, 3) = referenceMean .* ~crossSectionMaskVein + crossSectionMaskVein;

        f101 = figure(101);
        colormap("gray")
        f101.Position = [800, 200, 600, 600];
        ax101 = gca;
        ax101.Units = 'pixels';

    end

    videoM0_from_holowaves_norm = rescale(flatfieldM0);

    for frameIdx = 1:numFrames

        figure(100)

        hue_artery = 0 * (maskOnes(:, :, frameIdx) .* crossSectionMaskArtery);
        sat_artery = maskOnes(:, :, frameIdx) .* crossSectionMaskArtery;
        val_artery = maskOnes(:, :, frameIdx) .* crossSectionMaskArtery;

        tmp_frame_artery = hsv2rgb(hue_artery, sat_artery, val_artery) .* crossSectionMaskArtery;
        imgM0_artery = ones(size(tmp_frame_artery)) .* ~crossSectionMaskArtery .* videoM0_from_holowaves_norm(:, :, frameIdx) + tmp_frame_artery;
        imagesc(ax100, imgM0_artery);
        axis image
        axis off

        for section_idx = 1:nb_sections_artery
            new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);

            if round(avgBloodVolumeRateArtery(section_idx, frameIdx), 1) > 0
                text(new_x, new_y, string(round(avgBloodVolumeRateArtery(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            else
                text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            end

        end

        tmp_bvra = round(totalAvgBloodVolumeRateArtery(frameIdx));
        tmp_bvra = (tmp_bvra > 0) * tmp_bvra;
        title(sprintf("Blood volume rate : %02.0f µL/min in arteries", tmp_bvra));
        set(gca, 'FontSize', 18)
        drawnow

        Total_bloodVolumeRate_artery = getframe(f100);
        volume_rate_video_artery(:, :, :, frameIdx) = frame2im(Total_bloodVolumeRate_artery);

        if veins_analysis

            figure(101)

            hue_vein = 0.6 * (maskOnes(:, :, frameIdx) .* crossSectionMaskVein);
            sat_vein = maskOnes(:, :, frameIdx) .* crossSectionMaskVein;
            val_vein = maskOnes(:, :, frameIdx) .* crossSectionMaskVein;

            tmp_frame_vein = hsv2rgb(hue_vein, sat_vein, val_vein) .* crossSectionMaskVein;
            imgM0_vein = ones(size(tmp_frame_vein)) .* ~crossSectionMaskVein .* videoM0_from_holowaves_norm(:, :, frameIdx) + tmp_frame_vein;
            imagesc(ax101, imgM0_vein);
            axis image
            axis off

            for section_idx = 1:nb_sections_vein
                new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
                new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);

                if round(avgBloodVolumeRateVein(section_idx, frameIdx), 1) > 0

                    text(new_x, new_y, string(round(avgBloodVolumeRateVein(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");

                else

                    text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");

                end

            end

            tmp_bvrv = round(totalAvgBloodVolumeRateVein(frameIdx));
            tmp_bvrv = (tmp_bvrv > 0) * tmp_bvrv;
            title(sprintf('Blood volume rate : %02.0f µL/min in veins', tmp_bvrv));
            set(gca, 'FontSize', 18)
            drawnow

            Total_bloodVolumeRate_vein = getframe(f101);
            volume_rate_video_vein(:, :, :, frameIdx) = frame2im(Total_bloodVolumeRate_vein);

        end

    end
    
    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(volume_rate_video_artery),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volume_rate_video_artery.avi')));

    if veins_analysis
        parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(volume_rate_video_vein),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volume_rate_video_vein.avi')));

    end

    figure(102);
    imshow(maskRGB_artery);

    for section_idx = 1:nb_sections_artery
        new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);
        text(new_x, new_y, string(round(mean(avgBloodVolumeRateArtery(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

    title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(totalAvgBloodVolumeRateArtery))));
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    set(gca, 'FontSize', 14)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInArteries.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInArteries.eps')))
    % F_Total_bloodVolumeRate_artery = getframe(f102);

    if veins_analysis

        figure(103);
        imshow(maskRGB_vein);

        for section_idx = 1:nb_sections_vein
            new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);
            text(new_x, new_y, string(round(mean(avgBloodVolumeRateVein(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end

        title(sprintf("Total blood volume rate : %02.0f µL/min in veins", round(mean(totalAvgBloodVolumeRateVein))));
        drawnow
        ax = gca;
        ax.Units = 'pixels';
        set(gca, 'FontSize', 14)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInVeins.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInVeins.eps')))

    end

    %% Plot Volume Rate artery

    plot_vol_rate_artery = figure(104);
    plot_vol_rate_artery.Position = [200 200 600 275];

    %plot(fullTime, total_avg_bloodVolumeRate_artery, '-k', 'LineWidth', 2);
    curve1 = totalAvgBloodVolumeRateArtery + 0.5 * totalStdBloodVolumeRateArtery;
    curve2 = totalAvgBloodVolumeRateArtery - 0.5 * totalStdBloodVolumeRateArtery;
    fullTime2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)];

    fill(fullTime2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, totalAvgBloodVolumeRateArtery, '-k', 'LineWidth', 2);
    axis tight;
    hold off

    ax_vol_rate_artery = axis;
    ax_vol_rate_artery(3) = 0;

    ylabel('Blood volume rate (µL/min)')
    xlabel('Time (s)')
    title("Total blood volume rate in arteries")
    axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
    fontsize(gca, 14, "points");
    set(gca, 'Linewidth', 2)

    Plot_volume_rate_artery = getframe(gcf);
    Plot_volume_rate_video_artery = zeros(size(Plot_volume_rate_artery.cdata, 1), size(Plot_volume_rate_artery.cdata, 2), 3, numFrames);

    if veins_analysis

        plot_vol_rate_vein = figure(105);
        plot_vol_rate_vein.Position = [200 475 600 275];

        %plot(fullTime, total_avg_bloodVolumeRate_vein, '-k', 'LineWidth', 2);
        curve1 = totalAvgBloodVolumeRateVein + 0.5 * totalStdBloodVolumeRateVein;
        curve2 = totalAvgBloodVolumeRateVein - 0.5 * totalStdBloodVolumeRateVein;
        fullTime2 = [fullTime, fliplr(fullTime)];
        inBetween = [curve1, fliplr(curve2)];

        fill(fullTime2, inBetween, Color_std);
        hold on;
        plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, totalAvgBloodVolumeRateVein, '-k', 'LineWidth', 1);
        axis tight;
        hold off

        ax_vol_rate_vein = axis;
        ax_vol_rate_vein(3) = 0;

        axis tight;
        axis([ax_vol_rate_vein(1) ax_vol_rate_vein(2) ax_vol_rate_vein(3) ax_vol_rate_vein(4)]);
        Plot_volume_rate_vein = getframe(gcf);
        Plot_volume_rate_video_vein = zeros(size(Plot_volume_rate_vein.cdata, 1), size(Plot_volume_rate_vein.cdata, 2), 3, numFrames);

    end

    for frameIdx = 1:numFrames

        figure(104)

        curve1 = totalAvgBloodVolumeRateArtery(1:frameIdx) + 0.5 * totalStdBloodVolumeRateArtery(1:frameIdx);
        curve2 = totalAvgBloodVolumeRateArtery(1:frameIdx) - 0.5 * totalStdBloodVolumeRateArtery(1:frameIdx);
        tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
        inBetween = [curve1, fliplr(curve2)];

        fill(tmp_fullTime, inBetween, Color_std);
        hold on;
        plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), totalAvgBloodVolumeRateArtery(1:frameIdx), '-k', 'LineWidth', 2);
        yline(mean_volume_rate_artery, '--k', 'LineWidth', 2)
        hold off;

        ylabel('Blood volume rate (µL/min)')
        xlabel('Time (s)')
        title(sprintf("Total blood volume rate in arteries : %02.0f µL/min", round(mean_volume_rate_artery)))
        axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
        fontsize(gca, 14, "points");
        set(gca, 'Linewidth', 2)

        Plot_volume_rate_artery = getframe(gcf);
        Plot_volume_rate_video_artery(:, :, :, frameIdx) = Plot_volume_rate_artery.cdata;

        set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
        set(gca, 'Linewidth', 2)

        if veins_analysis

            figure(105);

            curve1 = totalAvgBloodVolumeRateVein(1:frameIdx) + 0.5 * totalStdBloodVolumeRateVein(1:frameIdx);
            curve2 = totalAvgBloodVolumeRateVein(1:frameIdx) - 0.5 * totalStdBloodVolumeRateVein(1:frameIdx);
            tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
            inBetween = [curve1, fliplr(curve2)];

            fill(tmp_fullTime, inBetween, Color_std);
            hold on;
            plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), totalAvgBloodVolumeRateVein(1:frameIdx), '-k', 'LineWidth', 2);
            yline(mean_volume_rate_vein, '--k', 'LineWidth', 2)
            hold off;

            ylabel('Blood volume rate (µL/min)')
            xlabel('Time (s)')
            title(sprintf("Total blood volume rate in veins : %02.0f µL/min", round(mean_volume_rate_vein)))
            axis([ax_vol_rate_vein(1) ax_vol_rate_vein(2) ax_vol_rate_vein(3) ax_vol_rate_vein(4)]);
            fontsize(gca, 14, "points");
            set(gca, 'Linewidth', 2)

            Plot_volume_rate_vein = getframe(gcf);
            Plot_volume_rate_video_vein(:, :, :, frameIdx) = Plot_volume_rate_vein.cdata;

            set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
            set(gca, 'Linewidth', 2)

        end

    end

    figure(104)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInArterySection.eps')))

    if veins_analysis

        figure(105)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInVeinSection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'bloodVolumeRateInVeinSection.eps')))

    end

    plot2txt(fullTime, totalAvgBloodVolumeRateArtery, 'TotalBloodVolumeRateArteryAVG', ToolBox)
    plot2txt(fullTime, totalStdBloodVolumeRateArtery, 'TotalBloodVolumeRateArterySTD', ToolBox)

    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(Plot_volume_rate_video_artery),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volume_rate_video_artery.avi')));


    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    
    combined_Gif_artery = cat(1, mat2gray(volume_rate_video_artery), mat2gray(Plot_volume_rate_video_artery));
    writeGifOnDisc(combined_Gif_artery,fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateArtery")),timePeriod);

    imwrite(mat2gray(Plot_volume_rate_video_artery(:,:,:,end)), fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, "bloodVolumeRateVideoArtery.png")));

    if veins_analysis

        plot2txt(fullTime, totalAvgBloodVolumeRateVein, 'TotalBloodVolumeRateVeinsAVG', ToolBox)
        plot2txt(fullTime, totalStdBloodVolumeRateVein, 'TotalBloodVolumeRateVeinsSTD', ToolBox)

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volume_rate_video_vein.avi')));

        tmp = mat2gray(Plot_volume_rate_video_vein);
        open(w)

        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateVein")), timePeriod, 0.04, numFrames);
        combined_Gif_vein = cat(1, mat2gray(volume_rate_video_vein), mat2gray(Plot_volume_rate_video_vein));

        for frameIdx = 1:numFrames
            writeVideo(w, tmp(:, :, :, frameIdx));
            gifWriter.write(combined_Gif_vein(:, :, :, frameIdx), frameIdx);
        end

        gifWriter.generate();
        gifWriter.delete();

        imwrite(mat2gray(Plot_volume_rate_video_vein(:,:,:,end)), fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername, "bloodVolumeRateVideoVein.png")));

        close(w);

    end

    %% txt file output with measured pulse wave parameters

    if veins_analysis
        fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
        fprintf(fileID, [ ...
                             'Value of total arterial blood volume rate (µL/min) :\n%d\n' ...
                         'Value of total venous blood volume rate (µL/min) :\n%d\n'], ...
            mean_volume_rate_artery, ...
            mean_volume_rate_vein);
        fclose(fileID);
    else
        fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
        fprintf(fileID, ...
            'Value of total arterial blood volume rate (µL/min) :\n%d\n', ...
            mean_volume_rate_artery);
        fclose(fileID);
    end

    for section_idx = 1:nb_sections_artery
        fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
        fprintf(fileID, [ ...
                             'Artery n°A%d : cross_section (mm^2) : \n %d \n' ...
                             'Artery n°A%d : vessel diameter (µm) : \n %d \n' ...
                             'Artery n°A%d : average velocity (mm/s) : \n %d \n' ...
                         'Artery n°A%d : blood rate (µL/min) : \n %d \n \n'], ...
            section_idx, ...
            crossSectionAreaArtery(section_idx), ...
            section_idx, ...
            2 * sqrt(crossSectionAreaArtery(section_idx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
            section_idx, ...
            avgBloodVelocityArtery(section_idx), ...
            section_idx, ...
            avgBloodVolumeRateArtery(section_idx)); % mm^3/s -> µL/min
        fclose(fileID);
    end

    if veins_analysis

        for section_idx = 1:nb_sections_vein
            fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
            fprintf(fileID, [ ...
                                 'Vein n°V%d : cross_section (mm^2) : \n %d \n ' ...
                                 'Vein n°V%d : vessel diameter (µm) : \n %d \n ' ...
                                 'Vein n°V%d : average velocity (mm/s) : \n %d \n ' ...
                             'Vein n°V%d : blood rate (µL/min) : \n %d \n \n'], ...
                section_idx, ...
                crossSectionAreaVein(section_idx), ...
                section_idx, ...
                2 * sqrt(crossSectionAreaVein(section_idx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
                section_idx, ...
                avgBloodVelocityVein(section_idx), ...
                section_idx, ...
                avgBloodVolumeRateVein(section_idx)); % mm^3/s -> µL/min
            fclose(fileID);
        end

    end

    close all
    
    %% All circles testing

    if ~ PW_params.AllCirclesFlag
        return
    end

    % for the all circles output
    nbCircles = PW_params.nbCircles;
    maskSectionCircles = cell(1,nbCircles);
    deltr = (PW_params.velocity_bigRadiusRatio - PW_params.velocity_smallRadiusRatio) * (numY + numX)/2 /nbCircles; %PW_params.radius_gap 
    for i = 1:nbCircles
        rad1 = (PW_params.velocity_smallRadiusRatio) * (numY + numX) / 2 + (i-1) * deltr ; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * deltr ;
        rad2 = rad1 + deltr ;
        c1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad1;
        c2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad2;
        maskSectionCircles(i) = {xor(c1,c2)};
        
        % save mask image
        createMaskSection(referenceMean , maskArtery,rad1,rad2,sprintf('_mask_artery_section_circle_%d.png',i), ToolBox, path);
    end
    close(156);

    % for all circles output 

    SubImg_locs_artery_Circles = zeros(nbCircles,nb_sections_artery, 2);
    SubImg_width_artery_Circles = zeros(nbCircles,nb_sections_artery, 1);

    for i = 1:nbCircles
        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        maskSectionArtery = bwlabel(maskSectionArtery);
    
        nb_sections_artery = max(maskSectionArtery, [], 'all');
        masksSectionsArtery = zeros(numX, numY, nb_sections_artery);
    
        parfor section_idx = 1:nb_sections_artery
    
            masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);
    
        end
    
    
        for section_idx = 1:nb_sections_artery
            [row, col] = find(masksSectionsArtery(:, :, section_idx));
            SubImg_locs_artery_Circles(i,section_idx, 1) = round(mean(row));
            SubImg_locs_artery_Circles(i,section_idx, 2) = round(mean(col));
            SubImg_width_artery_Circles(i,section_idx) = 0.01 * size(maskArtery, 1);
        end
    end

    % for all circles output 

    for i = 1:nbCircles
        [avgBloodVolumeRateArtery, ~, ~, avgBloodVelocityArtery, crossSectionMaskArtery, totalAvgBloodVolumeRateArtery, totalStdBloodVolumeRateArtery] = cross_section_analysis(reshape(nonzeros(SubImg_locs_artery_Circles(i,:,:)),[],2), nonzeros(SubImg_width_artery_Circles(i,:,:)), mask_artery, vRMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile,i);
        if length(avgBloodVolumeRateArtery)<1
            continue
        end
        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        maskSectionArtery = bwlabel(maskSectionArtery);
    
        nb_sections_artery = max(maskSectionArtery, [], 'all');

        labels_arteries = cell(nb_sections_artery, 1);
        strXlabel = 'Time(s)'; %createXlabelTime(1);
        strYlabel = 'Velocity (mm.s-1)';
    
        data_to_plot_artery = zeros(size(avgBloodVelocityArtery, 1), size(avgBloodVolumeRateArtery, 2));

        
    
        parfor section_idx = 1:nb_sections_artery
            data_to_plot_artery(section_idx, :) = smoothdata(avgBloodVelocityArtery(section_idx, :), 'lowess');
            labels_arteries{section_idx} = strcat('A', num2str(section_idx));
        end
    
        figure(57)
        plot(fullTime, data_to_plot_artery, 'LineWidth', 2)
        title('Blood velocity in artery sections');
        legend(labels_arteries);
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;
    
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'velocityInArterySection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'velocityInArterySection.eps')))
    
        data_to_plot_artery_all = squeeze(mean(data_to_plot_artery, 1));
        plot_artery_mean = mean(data_to_plot_artery_all);
    
        figure(571)
        plot(fullTime, data_to_plot_artery_all, 'k', 'LineWidth', 2)
        yline(plot_artery_mean, '--k', 'LineWidth', 2)
        title('Blood velocity in artery sections');
        legend(sprintf('All Artery Sections'));
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;
    
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'velocityInAllArterySection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'velocityInAllArterySection.eps')))
        

        Color_std = [0.7 0.7 0.7];
        figure_width = 600;
        figure_height = figure_width;
    
        maskRGB_artery = ones(size(maskArtery, 1), size(maskArtery, 2), 3);
        volume_rate_video_artery = zeros(figure_width, figure_height, 3, numFrames);
        mean_volume_rate_artery = mean(totalAvgBloodVolumeRateArtery);
    
        maskRGB_artery(:, :, 3) = referenceMean .* ~crossSectionMaskArtery;
        maskRGB_artery(:, :, 2) = referenceMean .* ~crossSectionMaskArtery;
        maskRGB_artery(:, :, 1) = referenceMean .* ~crossSectionMaskArtery + maskRGB_artery(:, :, 1) .* crossSectionMaskArtery;
    
        f100 = figure(100);
        colormap("gray")
        f100.Position = [200, 200, 600, 600];
        ax100 = gca;
        ax100.Units = 'pixels';
    
        
    
        videoM0_from_holowaves_norm = rescale(flatfieldM0);
    
        for frameIdx = numFrames:numFrames
    
            figure(100)
    
            hue_artery = 0 * (maskOnes(:, :, frameIdx) .* crossSectionMaskArtery);
            sat_artery = maskOnes(:, :, frameIdx) .* crossSectionMaskArtery;
            val_artery = maskOnes(:, :, frameIdx) .* crossSectionMaskArtery;
    
            tmp_frame_artery = hsv2rgb(hue_artery, sat_artery, val_artery) .* crossSectionMaskArtery;
            imgM0_artery = ones(size(tmp_frame_artery)) .* ~crossSectionMaskArtery .* videoM0_from_holowaves_norm(:, :, frameIdx) + tmp_frame_artery;
            imagesc(ax100, imgM0_artery);
            axis image
            axis off
    
            for section_idx = 1:nb_sections_artery
                new_x = x_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,2) - x_center);
                new_y = y_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,1) - y_center);
    
                if round(avgBloodVolumeRateArtery(section_idx, frameIdx), 1) > 0
                    text(new_x, new_y, string(round(avgBloodVolumeRateArtery(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                else
                    text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                end
    
            end
    
            tmp_bvra = round(totalAvgBloodVolumeRateArtery(frameIdx));
            tmp_bvra = (tmp_bvra > 0) * tmp_bvra;
            title(sprintf("Blood volume rate : %02.0f µL/min in arteries", tmp_bvra));
            set(gca, 'FontSize', 18)
            drawnow
    
            Total_bloodVolumeRate_artery = getframe(f100);
            volume_rate_video_artery(:, :, :, frameIdx) = frame2im(Total_bloodVolumeRate_artery);
    
            
    
        end
        
        parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(volume_rate_video_artery),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volume_rate_video_artery.avi')));
    
        figure(102);
        imshow(maskRGB_artery);
    
        for section_idx = 1:nb_sections_artery
            new_x = x_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,1) - y_center);
            text(new_x, new_y, string(round(mean(avgBloodVolumeRateArtery(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end
    
        title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(totalAvgBloodVolumeRateArtery))));
        drawnow
        ax = gca;
        ax.Units = 'pixels';
        set(gca, 'FontSize', 14)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'bloodVolumeRateInArteries.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'bloodVolumeRateInArteries.eps')))
        % F_Total_bloodVolumeRate_artery = getframe(f102);
    
        
        
        plot_vol_rate_artery = figure(104);
        plot_vol_rate_artery.Position = [200 200 600 275];
    
        %plot(fullTime, total_avg_bloodVolumeRate_artery, '-k', 'LineWidth', 2);
        curve1 = totalAvgBloodVolumeRateArtery + 0.5 * totalStdBloodVolumeRateArtery;
        curve2 = totalAvgBloodVolumeRateArtery - 0.5 * totalStdBloodVolumeRateArtery;
        fullTime2 = [fullTime, fliplr(fullTime)];
        inBetween = [curve1, fliplr(curve2)];
    
        fill(fullTime2, inBetween, Color_std);
        hold on;
        plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, totalAvgBloodVolumeRateArtery, '-k', 'LineWidth', 2);
        axis tight;
        hold off
    
        ax_vol_rate_artery = axis;
        ax_vol_rate_artery(3) = 0;
    
        ylabel('Blood volume rate (µL/min)')
        xlabel('Time (s)')
        title("Total blood volume rate in arteries")
        axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
        fontsize(gca, 14, "points");
        set(gca, 'Linewidth', 2)
    
        Plot_volume_rate_artery = getframe(gcf);
        Plot_volume_rate_video_artery = zeros(size(Plot_volume_rate_artery.cdata, 1), size(Plot_volume_rate_artery.cdata, 2), 3, numFrames);
    
        
    
        for frameIdx = numFrames:numFrames
    
            figure(104)
    
            curve1 = totalAvgBloodVolumeRateArtery(1:frameIdx) + 0.5 * totalStdBloodVolumeRateArtery(1:frameIdx);
            curve2 = totalAvgBloodVolumeRateArtery(1:frameIdx) - 0.5 * totalStdBloodVolumeRateArtery(1:frameIdx);
            tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
            inBetween = [curve1, fliplr(curve2)];
    
            fill(tmp_fullTime, inBetween, Color_std);
            hold on;
            plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), totalAvgBloodVolumeRateArtery(1:frameIdx), '-k', 'LineWidth', 2);
            yline(mean_volume_rate_artery, '--k', 'LineWidth', 2)
            hold off;
    
            ylabel('Blood volume rate (µL/min)')
            xlabel('Time (s)')
            title(sprintf("Total blood volume rate in arteries : %02.0f µL/min", round(mean_volume_rate_artery)))
            axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
            fontsize(gca, 14, "points");
            set(gca, 'Linewidth', 2)
    
            Plot_volume_rate_artery = getframe(gcf);
            Plot_volume_rate_video_artery(:, :, :, frameIdx) = Plot_volume_rate_artery.cdata;
    
            set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
            set(gca, 'Linewidth', 2)
    
            
    
        end
    
        figure(104)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'bloodVolumeRateInArterySection.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, 'bloodVolumeRateInArterySection.eps')))
    
        
    
        %plot2txt(fullTime, total_avg_bloodVolumeRate_artery, 'TotalBloodVolumeRateArteryAVG', ToolBox)
        %plot2txt(fullTime, total_std_bloodVolumeRate_artery, 'TotalBloodVolumeRateArterySTD', ToolBox)
    
        %parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(Plot_volume_rate_video_artery),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, 'circle_',num2str(i),'_plot_volume_rate_video_artery.avi')));
    
    
        %timePeriod = ToolBox.stride / ToolBox.fs / 1000;
        
        %combined_Gif_artery = cat(1, mat2gray(volume_rate_video_artery), mat2gray(Plot_volume_rate_video_artery));
        %writeGifOnDisc(combined_Gif_artery,fullfile(ToolBox.PW_path_gif, sprintf("%s_circle_%d_%s.gif", ToolBox.PW_folder_name,i, "volumeRateArtery")),timePeriod);
    
        imwrite(mat2gray(Plot_volume_rate_video_artery(:,:,:,end)), fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i, "bloodVolumeRateVideoArtery.png")));


    end

end
