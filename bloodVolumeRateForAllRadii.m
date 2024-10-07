function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, ~, v_RMS, dataM0, videoM0_from_holowaves, ToolBox, k, path, flagBloodVelocityProfile)

    PW_params = Parameters_json(path);

    veins_analysis = PW_params.veins_analysis;

    mkdir(ToolBox.PW_path_png, 'bloodVolumeRate')
    mkdir(ToolBox.PW_path_eps, 'bloodVolumeRate')

    [N, M] = size(maskArtery);
    N_frame = size(v_RMS, 3);
    [x, y] = meshgrid(1:M, 1:N);

    %maskArtery = imdilate(maskArtery, strel('disk', 5));
    %FIXME function velocity map

    v_RMS_AVG = mean(v_RMS, 3);
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);

    %% change mask section ?

    % radius_ratio = 0.18;
    % radius_gap = 0.05;
    % radius1 = (radius_ratio-radius_gap)* (M+N)/2;
    % radius2 = (radius_ratio+radius_gap)* (M+N)/2;

    radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (M + N) / 2;
    radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (M + N) / 2;
    cercle_mask1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius1;
    cercle_mask2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius2;

    maskSection = xor(cercle_mask1, cercle_mask2);
    dataM0 = mat2gray(dataM0);
    maskOnes = ones(size(maskArtery, 1), size(maskArtery, 2), size(dataM0, 3));
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);
    mean_M0 = rescale(mean(videoM0_from_holowaves, 3));

    maskSectionArtery = maskSection .* maskArtery;
    mask_artery = maskArtery;
    mask_vein = maskVein;
    ratio_etiquette = 1.2;

    figure(111)
    imagesc(maskSectionArtery .* v_RMS_AVG);

    maskSectionArtery = bwlabel(maskSectionArtery);

    figure(222)
    imagesc(maskSectionArtery)
    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    nb_sections_artery = max(maskSectionArtery, [], 'all');
    masksSectionsArtery = zeros(size(maskArtery, 1), size(maskArtery, 2), nb_sections_artery);

    parfor section_idx = 1:nb_sections_artery

        masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);

    end
    
    %% All circles testing

    if ~ PW_params.AllCirclesFlag
        return
    end

    % for the all circles output
    nbCircles = PW_params.nbCircles;
    maskSectionCircles = cell(1,nbCircles);
    deltr = (PW_params.velocity_bigRadiusRatio - PW_params.velocity_smallRadiusRatio) * (M + N)/2 /nbCircles; %PW_params.radius_gap 
    for i = 1:nbCircles
        rad1 = (PW_params.velocity_smallRadiusRatio) * (M + N) / 2 + (i-1) * deltr ; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * deltr ;
        rad2 = rad1 + deltr ;
        c1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad1;
        c2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad2;
        maskSectionCircles(i) = {xor(c1,c2)};
        
        % save mask image
        createMaskSection(mean_M0 , maskArtery,rad1,rad2,sprintf('_mask_artery_section_circle_%d.png',i), ToolBox, path);
    end
    close(156);

    % for all circles output 

    SubImg_locs_artery_Circles = zeros(nbCircles,nb_sections_artery, 2);
    SubImg_width_artery_Circles = zeros(nbCircles,nb_sections_artery, 1);

    for i = 1:nbCircles
        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        maskSectionArtery = bwlabel(maskSectionArtery);
    
        nb_sections_artery = max(maskSectionArtery, [], 'all');
        masksSectionsArtery = zeros(size(maskArtery, 1), size(maskArtery, 2), nb_sections_artery);
    
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
        [avg_bloodVolumeRate_artery, std_bloodVolumeRate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery, total_avg_bloodVolumeRate_artery, total_std_bloodVolumeRate_artery] = cross_section_analysis(reshape(nonzeros(SubImg_locs_artery_Circles(i,:,:)),[],2), nonzeros(SubImg_width_artery_Circles(i,:,:)), mask_artery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile,i);
        if length(avg_bloodVolumeRate_artery)<1
            continue
        end
        avg_bloodVolumeRate_artery_r(i) = avg_bloodVolumeRate_artery;
        std_bloodVolumeRate_artery_r(i) = std_bloodVolumeRate_artery;
        cross_section_area_artery_r(i) = cross_section_area_artery_r;
        cross_section_mask_artery_r(i) = cross_section_mask_artery_r;

        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        maskSectionArtery = bwlabel(maskSectionArtery);
    
        nb_sections_artery = max(maskSectionArtery, [], 'all');

        labels_arteries = cell(nb_sections_artery, 1);
        strXlabel = 'Time(s)'; %createXlabelTime(1);
        strYlabel = 'Velocity (mm.s-1)';
    
        data_to_plot_artery = zeros(size(avg_blood_velocity_artery, 1), size(avg_bloodVolumeRate_artery, 2));

        
    
        parfor section_idx = 1:nb_sections_artery
            data_to_plot_artery(section_idx, :) = smoothdata(avg_blood_velocity_artery(section_idx, :), 'lowess');
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
        volume_rate_video_artery = zeros(figure_width, figure_height, 3, N_frame);
        mean_volume_rate_artery = mean(total_avg_bloodVolumeRate_artery);
    
        maskRGB_artery(:, :, 3) = mean_M0 .* ~cross_section_mask_artery;
        maskRGB_artery(:, :, 2) = mean_M0 .* ~cross_section_mask_artery;
        maskRGB_artery(:, :, 1) = mean_M0 .* ~cross_section_mask_artery + maskRGB_artery(:, :, 1) .* cross_section_mask_artery;
    
        f100 = figure(100);
        colormap("gray")
        f100.Position = [200, 200, 600, 600];
        ax100 = gca;
        ax100.Units = 'pixels';
    
        
    
        videoM0_from_holowaves_norm = rescale(videoM0_from_holowaves);
    
        for frameIdx = N_frame:N_frame
    
            figure(100)
    
            hue_artery = 0 * (maskOnes(:, :, frameIdx) .* cross_section_mask_artery);
            sat_artery = maskOnes(:, :, frameIdx) .* cross_section_mask_artery;
            val_artery = maskOnes(:, :, frameIdx) .* cross_section_mask_artery;
    
            tmp_frame_artery = hsv2rgb(hue_artery, sat_artery, val_artery) .* cross_section_mask_artery;
            imgM0_artery = ones(size(tmp_frame_artery)) .* ~cross_section_mask_artery .* videoM0_from_holowaves_norm(:, :, frameIdx) + tmp_frame_artery;
            imagesc(ax100, imgM0_artery);
            axis image
            axis off
    
            for section_idx = 1:nb_sections_artery
                new_x = x_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,2) - x_center);
                new_y = y_center + ratio_etiquette * (SubImg_locs_artery_Circles(i,section_idx,1) - y_center);
    
                if round(avg_bloodVolumeRate_artery(section_idx, frameIdx), 1) > 0
                    text(new_x, new_y, string(round(avg_bloodVolumeRate_artery(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                else
                    text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                end
    
            end
    
            tmp_bvra = round(total_avg_bloodVolumeRate_artery(frameIdx));
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
            text(new_x, new_y, string(round(mean(avg_bloodVolumeRate_artery(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end
    
        title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(total_avg_bloodVolumeRate_artery))));
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
        curve1 = total_avg_bloodVolumeRate_artery + 0.5 * total_std_bloodVolumeRate_artery;
        curve2 = total_avg_bloodVolumeRate_artery - 0.5 * total_std_bloodVolumeRate_artery;
        fullTime2 = [fullTime, fliplr(fullTime)];
        inBetween = [curve1, fliplr(curve2)];
    
        fill(fullTime2, inBetween, Color_std);
        hold on;
        plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, total_avg_bloodVolumeRate_artery, '-k', 'LineWidth', 2);
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
        Plot_volume_rate_video_artery = zeros(size(Plot_volume_rate_artery.cdata, 1), size(Plot_volume_rate_artery.cdata, 2), 3, N_frame);
    
        
    
        for frameIdx = N_frame:N_frame
    
            figure(104)
    
            curve1 = total_avg_bloodVolumeRate_artery(1:frameIdx) + 0.5 * total_std_bloodVolumeRate_artery(1:frameIdx);
            curve2 = total_avg_bloodVolumeRate_artery(1:frameIdx) - 0.5 * total_std_bloodVolumeRate_artery(1:frameIdx);
            tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
            inBetween = [curve1, fliplr(curve2)];
    
            fill(tmp_fullTime, inBetween, Color_std);
            hold on;
            plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), total_avg_bloodVolumeRate_artery(1:frameIdx), '-k', 'LineWidth', 2);
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
