function [] = flow_rate(maskArtery, maskVein, ~, ~, v_RMS, dataM0, videoM0_from_holowaves, ToolBox, k, path)

    PW_params = Parameters_json(path);

    veins_analysis = PW_params.veins_analysis;

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

    %% Find the locations of the sections in arteries

    maskSectionArtery = maskSection .* maskArtery;

    figure(111)
    imagesc(maskSectionArtery .* v_RMS_AVG);

    maskSectionArtery = bwlabel(maskSectionArtery);

    figure(222)
    imagesc(maskSectionArtery)

    nb_sections_artery = max(maskSectionArtery, [], 'all');
    masksSectionsArtery = zeros(size(maskArtery, 1), size(maskArtery, 2), nb_sections_artery);

    for section_idx = 1:nb_sections_artery

        masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);

    end

    SubImg_locs_artery = zeros(nb_sections_artery, 2);
    SubImg_width_artery = zeros(nb_sections_artery, 1);

    for section_idx = 1:nb_sections_artery
        [row, col] = find(masksSectionsArtery(:, :, section_idx));
        SubImg_locs_artery(section_idx, 1) = round(mean(row));
        SubImg_locs_artery(section_idx, 2) = round(mean(col));
        SubImg_width_artery(section_idx) = 0.01 * size(maskArtery, 1);
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
        masksSectionsVein = zeros(size(maskVein, 1), size(maskVein, 2), nb_sections_vein);

        for section_idx = 1:nb_sections_vein
            masksSectionsVein(:, :, section_idx) = (maskSectionVein == section_idx);
        end

        SubImg_locs_vein = zeros(nb_sections_vein, 2);
        SubImg_width_vein = zeros(nb_sections_vein, 1);

        for section_idx = 1:nb_sections_vein
            [row, col] = find(masksSectionsVein(:, :, section_idx));
            SubImg_locs_vein(section_idx, 1) = round(mean(row));
            SubImg_locs_vein(section_idx, 2) = round(mean(col));
            SubImg_width_vein(section_idx) = 0.01 * size(maskVein, 1);
        end

    end

    %% Compute blood volume rate

    mask_artery = maskArtery;
    mask_vein = maskVein;

    if veins_analysis
        [avg_blood_volume_rate_vein, std_blood_volume_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein, total_avg_blood_volume_rate_vein, total_std_blood_volume_rate_vein] = cross_section_analysis(SubImg_locs_vein, SubImg_width_vein, mask_vein, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'vein');
    end

    [avg_blood_volume_rate_artery, std_blood_volume_rate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery, total_avg_blood_volume_rate_artery, total_std_blood_volume_rate_artery] = cross_section_analysis(SubImg_locs_artery, SubImg_width_artery, mask_artery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery');

    labels_arteries = cell(nb_sections_artery, 1);
    strXlabel = 'Time(s)'; %createXlabelTime(1);
    strYlabel = 'Velocity (mm.s-1)';

    data_to_plot_artery = zeros(size(avg_blood_velocity_artery, 1), size(avg_blood_volume_rate_artery, 2));

    for section_idx = 1:nb_sections_artery
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

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_artery_section.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_artery_section.eps')))

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

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_all_artery_section.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_all_artery_section.eps')))

    if veins_analysis

        labels_veins = cell(nb_sections_vein, 1);

        data_to_plot_vein = zeros(size(avg_blood_velocity_vein, 1), size(avg_blood_volume_rate_vein, 2));

        for section_idx = 1:nb_sections_vein
            data_to_plot_vein(section_idx, :) = smoothdata(avg_blood_velocity_vein(section_idx, :), 'lowess');
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

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_vein_section.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_vein_section.eps')))

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

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_all_vein_section.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_velocity_in_all_vein_section.eps')))
    end

    dataM0 = mat2gray(dataM0);
    maskOnes = ones(size(maskArtery, 1), size(maskArtery, 2), size(dataM0, 3));
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);
    mean_M0 = rescale(mean(videoM0_from_holowaves, 3));
    ratio_etiquette = 1.2;
    %% Saving Vein and artery data in txt

    for section_idx = 1:nb_sections_artery
        plot2txt(fullTime, avg_blood_volume_rate_artery(section_idx, :), strcat('blood_volume_rate_artery_A', num2str(section_idx)), ToolBox)
        plot2txt(fullTime, std_blood_volume_rate_artery(section_idx, :), strcat('blood_volume_rate_artery_std_A', num2str(section_idx)), ToolBox)
        plot2txt(fullTime, avg_blood_velocity_artery(section_idx, :), strcat('avg_velocity_artery_A', num2str(section_idx)), ToolBox)

    end

    if veins_analysis

        for section_idx = 1:nb_sections_vein
            plot2txt(fullTime, avg_blood_volume_rate_vein(section_idx, :), strcat('blood_volume_rate_vein_V', num2str(section_idx)), ToolBox)
            plot2txt(fullTime, std_blood_volume_rate_vein(section_idx, :), strcat('blood_volume_rate_vein_std_V', num2str(section_idx)), ToolBox)
            plot2txt(fullTime, avg_blood_velocity_vein(section_idx, :), strcat('avg_velocity_vein_V', num2str(section_idx)), ToolBox)

        end

    end

    %% Vein and artery numerotation

    maskRGB = ones(size(maskArtery, 1), size(maskArtery, 2), 3);

    if veins_analysis
        maskRGB(:, :, 3) = (mean_M0 .* ~cross_section_mask_vein + maskRGB(:, :, 3) .* cross_section_mask_vein) .* ~cross_section_mask_artery;
        maskRGB(:, :, 2) = mean_M0 .* ~(cross_section_mask_artery + cross_section_mask_vein) + zeros(size(maskRGB(:, :, 2))) .* (cross_section_mask_artery + cross_section_mask_vein);
        maskRGB(:, :, 1) = (mean_M0 .* ~cross_section_mask_artery + maskRGB(:, :, 1) .* cross_section_mask_artery) .* ~cross_section_mask_vein;

    else
        maskRGB(:, :, 3) = mean_M0 .* ~cross_section_mask_artery + zeros(size(maskRGB(:, :, 3))) .* cross_section_mask_artery;
        maskRGB(:, :, 2) = mean_M0 .* ~cross_section_mask_artery + zeros(size(maskRGB(:, :, 2))) .* cross_section_mask_artery;
        maskRGB(:, :, 1) = mean_M0 .* ~cross_section_mask_artery + maskRGB(:, :, 1) .* cross_section_mask_artery;
    end

    figure(652)
    imshow(maskRGB);

    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for section_idx = 1:nb_sections_artery
        new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);
        text(new_x, new_y, strcat("A", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

    if veins_analysis

        for section_idx = 1:nb_sections_vein
            new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);
            text(new_x, new_y, strcat("V", num2str(section_idx)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end

    end

    drawnow
    ax = gca;
    ax.Units = 'pixels';
    marg = 30;
    pos = ax.Position;
    rect = [-marg, -marg, pos(3) + 2 * marg, pos(4) + 2 * marg];
    F_numerotation = getframe(ax, rect);

    %% Volume Rate FigurePosition

    Color_std = [0.7 0.7 0.7];
    figure_width = 600;
    figure_height = figure_width;

    maskRGB_artery = ones(size(maskArtery, 1), size(maskArtery, 2), 3);
    volume_rate_video_artery = zeros(figure_width, figure_height, 3, N_frame);
    mean_volume_rate_artery = mean(total_avg_blood_volume_rate_artery);

    maskRGB_artery(:, :, 3) = mean_M0 .* ~cross_section_mask_artery;
    maskRGB_artery(:, :, 2) = mean_M0 .* ~cross_section_mask_artery;
    maskRGB_artery(:, :, 1) = mean_M0 .* ~cross_section_mask_artery + maskRGB_artery(:, :, 1) .* cross_section_mask_artery;

    f100 = figure(100);
    colormap("gray")
    f100.Position = [200, 200, 600, 600];
    ax100 = gca;
    ax100.Units = 'pixels';

    if veins_analysis

        maskRGB_vein = ones(size(maskVein, 1), size(maskVein, 2), 3);
        volume_rate_video_vein = zeros(figure_width, figure_height, 3, N_frame);
        mean_volume_rate_vein = mean(total_avg_blood_volume_rate_vein);

        maskRGB_vein(:, :, 1) = mean_M0 .* ~cross_section_mask_vein;
        maskRGB_vein(:, :, 2) = mean_M0 .* ~cross_section_mask_vein;
        maskRGB_vein(:, :, 3) = mean_M0 .* ~cross_section_mask_vein + maskRGB_vein(:, :, 3) .* cross_section_mask_vein;

        f101 = figure(101);
        colormap("gray")
        f101.Position = [800, 200, 600, 600];
        ax101 = gca;
        ax101.Units = 'pixels';

    end

    videoM0_from_holowaves_norm = rescale(videoM0_from_holowaves);

    for frameIdx = 1:N_frame

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
            new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);

            if round(avg_blood_volume_rate_artery(section_idx, frameIdx), 1) > 0
                text(new_x, new_y, string(round(avg_blood_volume_rate_artery(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            else
                text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            end

        end

        tmp_bvra = round(total_avg_blood_volume_rate_artery(frameIdx));
        tmp_bvra = (tmp_bvra > 0) * tmp_bvra;
        title(sprintf("Blood volume rate : %02.0f µL/min in arteries", tmp_bvra));
        set(gca, 'FontSize', 18)
        drawnow

        Total_blood_volume_rate_artery = getframe(f100);
        volume_rate_video_artery(:, :, :, frameIdx) = frame2im(Total_blood_volume_rate_artery);

        if veins_analysis

            figure(101)

            hue_vein = 0.6 * (maskOnes(:, :, frameIdx) .* cross_section_mask_vein);
            sat_vein = maskOnes(:, :, frameIdx) .* cross_section_mask_vein;
            val_vein = maskOnes(:, :, frameIdx) .* cross_section_mask_vein;

            tmp_frame_vein = hsv2rgb(hue_vein, sat_vein, val_vein) .* cross_section_mask_vein;
            imgM0_vein = ones(size(tmp_frame_vein)) .* ~cross_section_mask_vein .* videoM0_from_holowaves_norm(:, :, frameIdx) + tmp_frame_vein;
            imagesc(ax101, imgM0_vein);
            axis image
            axis off

            for section_idx = 1:nb_sections_vein
                new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
                new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);

                if round(avg_blood_volume_rate_vein(section_idx, frameIdx), 1) > 0

                    text(new_x, new_y, string(round(avg_blood_volume_rate_vein(section_idx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");

                else

                    text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");

                end

            end

            tmp_bvrv = round(total_avg_blood_volume_rate_vein(frameIdx));
            tmp_bvrv = (tmp_bvrv > 0) * tmp_bvrv;
            title(sprintf('Blood volume rate : %02.0f µL/min in veins', tmp_bvrv));
            set(gca, 'FontSize', 18)
            drawnow

            Total_blood_volume_rate_vein = getframe(f101);
            volume_rate_video_vein(:, :, :, frameIdx) = frame2im(Total_blood_volume_rate_vein);

        end

    end

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volume_rate_video_artery.avi')));
    tmp = mat2gray(volume_rate_video_artery);
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);

    if veins_analysis

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volume_rate_video_vein.avi')));
        tmp = mat2gray(volume_rate_video_vein);
        open(w)

        for frameIdx = 1:N_frame
            writeVideo(w, tmp(:, :, :, frameIdx));
        end

        close(w);

    end

    figure(102);
    imshow(maskRGB_artery);

    for section_idx = 1:nb_sections_artery
        new_x = x_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_artery(section_idx, 1) - y_center);
        text(new_x, new_y, string(round(mean(avg_blood_volume_rate_artery(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

    title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(total_avg_blood_volume_rate_artery))));
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    set(gca, 'FontSize', 14)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'blood_volume_rate_in_arteries.png')))
    % F_Total_blood_volume_rate_artery = getframe(f102);

    if veins_analysis

        figure(103);
        imshow(maskRGB_vein);

        for section_idx = 1:nb_sections_vein
            new_x = x_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_vein(section_idx, 1) - y_center);
            text(new_x, new_y, string(round(mean(avg_blood_volume_rate_vein(section_idx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        end

        title(sprintf("Total blood volume rate : %02.0f µL/min in veins", round(mean(total_avg_blood_volume_rate_vein))));
        drawnow
        ax = gca;
        ax.Units = 'pixels';
        set(gca, 'FontSize', 14)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'blood_volume_rate_in_veins.png')))

    end

    %% Plot Volume Rate artery

    plot_vol_rate_artery = figure(104);
    plot_vol_rate_artery.Position = [200 200 600 275];

    %plot(fullTime, total_avg_blood_volume_rate_artery, '-k', 'LineWidth', 2);
    curve1 = total_avg_blood_volume_rate_artery + 0.5 * total_std_blood_volume_rate_artery;
    curve2 = total_avg_blood_volume_rate_artery - 0.5 * total_std_blood_volume_rate_artery;
    fullTime2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1', fliplr(curve2')];

    fill(fullTime2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, total_avg_blood_volume_rate_artery, '-k', 'LineWidth', 2);
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

    if veins_analysis

        plot_vol_rate_vein = figure(105);
        plot_vol_rate_vein.Position = [200 475 600 275];

        %plot(fullTime, total_avg_blood_volume_rate_vein, '-k', 'LineWidth', 2);
        curve1 = total_avg_blood_volume_rate_vein + 0.5 * total_std_blood_volume_rate_vein;
        curve2 = total_avg_blood_volume_rate_vein - 0.5 * total_std_blood_volume_rate_vein;
        fullTime2 = [fullTime, fliplr(fullTime)];
        inBetween = [curve1', fliplr(curve2')];

        fill(fullTime2, inBetween, Color_std);
        hold on;
        plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime, total_avg_blood_volume_rate_vein, '-k', 'LineWidth', 1);
        axis tight;
        hold off

        ax_vol_rate_vein = axis;
        ax_vol_rate_vein(3) = 0;

        axis tight;
        axis([ax_vol_rate_vein(1) ax_vol_rate_vein(2) ax_vol_rate_vein(3) ax_vol_rate_vein(4)]);
        Plot_volume_rate_vein = getframe(gcf);
        Plot_volume_rate_video_vein = zeros(size(Plot_volume_rate_vein.cdata, 1), size(Plot_volume_rate_vein.cdata, 2), 3, N_frame);

    end

    for frameIdx = 1:N_frame

        figure(104)

        curve1 = total_avg_blood_volume_rate_artery(1:frameIdx) + 0.5 * total_std_blood_volume_rate_artery(1:frameIdx);
        curve2 = total_avg_blood_volume_rate_artery(1:frameIdx) - 0.5 * total_std_blood_volume_rate_artery(1:frameIdx);
        tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
        inBetween = [curve1', fliplr(curve2')];

        fill(tmp_fullTime, inBetween, Color_std);
        hold on;
        plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), total_avg_blood_volume_rate_artery(1:frameIdx), '-k', 'LineWidth', 2);
        yline(mean_volume_rate_artery, '--k', 'LineWidth', 2)
        hold off;

        ylabel('Blood volume rate (µL/min)')
        xlabel('Time (s)')
        title(sprintf("Total blood volume rate in arteries : %02.0f", round(mean_volume_rate_artery)))
        axis([ax_vol_rate_artery(1) ax_vol_rate_artery(2) ax_vol_rate_artery(3) ax_vol_rate_artery(4)]);
        fontsize(gca, 14, "points");
        set(gca, 'Linewidth', 2)

        Plot_volume_rate_artery = getframe(gcf);
        Plot_volume_rate_video_artery(:, :, :, frameIdx) = Plot_volume_rate_artery.cdata;

        set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
        set(gca, 'Linewidth', 2)

        if veins_analysis

            figure(105);

            curve1 = total_avg_blood_volume_rate_vein(1:frameIdx) + 0.5 * total_std_blood_volume_rate_vein(1:frameIdx);
            curve2 = total_avg_blood_volume_rate_vein(1:frameIdx) - 0.5 * total_std_blood_volume_rate_vein(1:frameIdx);
            tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
            inBetween = [curve1', fliplr(curve2')];

            fill(tmp_fullTime, inBetween, Color_std);
            hold on;
            plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), total_avg_blood_volume_rate_vein(1:frameIdx), '-k', 'LineWidth', 2);
            yline(mean_volume_rate_vein, '--k', 'LineWidth', 2)
            hold off;

            ylabel('Blood volume rate (µL/min)')
            xlabel('Time (s)')
            title(sprintf("Total blood volume rate in veins : %02.0f", round(mean_volume_rate_vein)))
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
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_blood_volume_rate_in_artery_section.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_blood_volume_rate_in_artery_section.eps')))

    if veins_analysis

        figure(105)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_blood_volume_rate_in_vein_section.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, sprintf("%s_%s", ToolBox.main_foldername, 'Plot_blood_volume_rate_in_vein_section.eps')))

    end

    plot2txt(fullTime, total_avg_blood_volume_rate_artery, 'TotalBloodVolumeRateArteryAVG', ToolBox)
    plot2txt(fullTime, total_std_blood_volume_rate_artery, 'TotalBloodVolumeRateArterySTD', ToolBox)

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volume_rate_video_artery.avi')));
    tmp = mat2gray(Plot_volume_rate_video_artery);
    open(w)

    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateArtery")), timePeriod, 0.04, N_frame);
    combined_Gif_artery = cat(1, mat2gray(volume_rate_video_artery), mat2gray(Plot_volume_rate_video_artery));

    for frameIdx = 1:N_frame

        writeVideo(w, tmp(:, :, :, frameIdx));
        gifWriter.write(combined_Gif_artery(:, :, :, frameIdx), frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();

    imwrite(mat2gray(Plot_volume_rate_video_artery(:,:,:,end)), fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, "plot_volume_rate_video_artery.png")));

    if veins_analysis

        plot2txt(fullTime, total_avg_blood_volume_rate_vein, 'TotalBloodVolumeRateVeinsAVG', ToolBox)
        plot2txt(fullTime, total_std_blood_volume_rate_vein, 'TotalBloodVolumeRateVeinsSTD', ToolBox)

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volume_rate_video_vein.avi')));

        tmp = mat2gray(Plot_volume_rate_video_vein);
        open(w)

        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateVein")), timePeriod, 0.04, N_frame);
        combined_Gif_vein = cat(1, mat2gray(volume_rate_video_vein), mat2gray(Plot_volume_rate_video_vein));

        for frameIdx = 1:N_frame
            writeVideo(w, tmp(:, :, :, frameIdx));
            gifWriter.write(combined_Gif_vein(:, :, :, frameIdx), frameIdx);
        end

        gifWriter.generate();
        gifWriter.delete();

        imwrite(mat2gray(Plot_volume_rate_video_vein(:,:,:,end)), fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, "_plot_volume_rate_video_vein.png")));

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
            cross_section_area_artery(section_idx), ...
            section_idx, ...
            2 * sqrt(cross_section_area_artery(section_idx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
            section_idx, ...
            avg_blood_velocity_artery(section_idx), ...
            section_idx, ...
            avg_blood_volume_rate_artery(section_idx)); % mm^3/s -> µL/min
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
                cross_section_area_vein(section_idx), ...
                section_idx, ...
                2 * sqrt(cross_section_area_vein(section_idx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
                section_idx, ...
                avg_blood_velocity_vein(section_idx), ...
                section_idx, ...
                avg_blood_volume_rate_vein(section_idx)); % mm^3/s -> µL/min
            fclose(fileID);
        end

    end

    %% Saving figures

    %print('-f111222', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Artery_Sections.png'))) ;
    %print('-f111223', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Vein_Sections.png'))) ;
    %print('-f121', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_mean_blood_volume_rate_in_arteries.png'))) ;
    % print('-f102', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Plot_blood_volume_rate_in_arteries.png')));
    %print('-f122', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_mean_blood_volume_rate_in_veins.png'))) ;
    %print('-f115', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_AVG_Velocity_in_arteriy_sections.png'))) ;
    %print('-f120', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Total_blood_volume_rate_in_arteries.png'))) ;

    imwrite(F_numerotation.cdata, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_arteries_veins_numerotation.png')));

    close all

end
