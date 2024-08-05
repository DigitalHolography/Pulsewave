function [] = velocity(v_RMS_all, maskArtery, maskVein, videoM0, FreqVideoRGB, ToolBox, path)

    PW_params = Parameters_json(path);
    veins_analysis = PW_params.veins_analysis;

    % TRUE MIN and MAX V_RMS but not realistic
    Im = rescale(mean(videoM0, 3));
    [Nx, Ny, N_frame] = size(v_RMS_all);

    % Ones = ones(size(v_RMS));
    % V = mean(v_RMS,3);
    % Vmax_Arteries = max(V.*maskArtery,[],'all');
    % Vmax_Veins = max(V.*maskVein,[],'all');
    % Vmin_Arteries = min(V.*maskArtery+Vmax_Arteries*Ones.*(~maskArtery),[],'all');
    % Vmin_Veins = min(V.*maskVein+Vmax_Veins*Ones.*(~maskVein),[],'all');

    v_artery = sum(v_RMS_all .* maskArtery, [1 2]) / nnz(maskArtery);
    v_vein = sum(v_RMS_all .* maskVein, [1 2]) / nnz(maskVein);
    Vmax_Arteries = max(v_artery(:));
    Vmax_Veins = max(v_vein(:));
    Vmin_Arteries = min(v_artery(:));
    Vmin_Veins = min(v_vein(:));

    %% Construct velocity map

    ImgM0 = rescale(mean(videoM0, 3));

    if veins_analysis
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(Im, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        [hue_vein, sat_vein, val_vein, cmap_vein] = createHSVmap(Im, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
        val = Im .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein;
        flowImageRGB = hsv2rgb(hue_artery + hue_vein, sat_artery + sat_vein, val);
        flowImageRGB = flowImageRGB .* (maskArtery + maskVein) + ones(Nx, Ny) .* ImgM0 .* ~(maskArtery + maskVein);
    else
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(Im, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val = Im .* (~(maskArtery)) + val_artery .* maskArtery;
        flowImageRGB = hsv2rgb(hue_artery, sat_artery, val);
        flowImageRGB = flowImageRGB .* (maskArtery) + ones(Nx, Ny) .* ImgM0 .* ~(maskArtery);
    end

    figure(321)
    imshow(flowImageRGB)
    title('Video M0')

    %% Construct Velocity video
    flowVideoRGB = zeros(Nx, Ny, 3, N_frame);
    videoM0_norm = rescale(videoM0);
    v_mean = mean(v_RMS_all, 3);
    M0_norm_mean = mean(videoM0_norm, 3);

    if veins_analysis

        [hue_artery_mean, sat_artery_mean, val_artery_mean, cmap_artery] = createHSVmap(v_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        [hue_vein_mean, sat_vein_mean, val_vein_mean, cmap_vein] = createHSVmap(v_mean, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
        val_mean = v_mean .* (~(maskArtery + maskVein)) + val_artery_mean .* maskArtery + val_vein_mean .* maskVein;
        flowVideoRGB_mean = hsv2rgb(hue_artery_mean + hue_vein_mean, sat_artery_mean + sat_vein_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery + maskVein) + ones(Nx, Ny, 3) .* ~(maskArtery + maskVein) .* M0_norm_mean;
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, "v_RMS_RGB.png")))

        for frame_idx = 1:N_frame
            v = mat2gray(v_RMS_all(:, :, frame_idx));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            [hue_vein, sat_vein, val_vein, ~] = createHSVmap(v, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
            val = v .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein;
            flowVideoRGB(:, :, :, frame_idx) = hsv2rgb(hue_artery + hue_vein, sat_artery + sat_vein, val);
            flowVideoRGB(:, :, :, frame_idx) = flowVideoRGB(:, :, :, frame_idx) .* (maskArtery + maskVein) + ones(Nx, Ny, 3) .* ~(maskArtery + maskVein) .* videoM0_norm(:, :, frame_idx);
        end

    else

        [hue_artery_mean, sat_artery_mean, val_artery_mean, cmap_artery] = createHSVmap(v_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val_mean = v_mean .* (~(maskArtery)) + val_artery_mean .* maskArtery;
        flowVideoRGB_mean = hsv2rgb(hue_artery_mean, sat_artery_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery) + ones(Nx, Ny, 3) .* ~(maskArtery) .* M0_norm_mean;
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, "v_RMS_RGB.png")))

        for frame_idx = 1:N_frame
            v = mat2gray(v_RMS_all(:, :, frame_idx));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            val = v .* (~(maskArtery)) + val_artery .* maskArtery;
            flowVideoRGB(:, :, :, frame_idx) = hsv2rgb(hue_artery, sat_artery, val);
            flowVideoRGB(:, :, :, frame_idx) = flowVideoRGB(:, :, :, frame_idx) .* (maskArtery) + ones(Nx, Ny, 3) .* ~(maskArtery) .* videoM0_norm(:, :, frame_idx);
        end

    end

    clear hue_artery sat_artery val_artery hue_vein sat_vein val_vein

    v_histo_artery = round(v_RMS_all .* maskArtery);
    v_min = min(v_histo_artery, [], 'all');
    v_max = max(v_histo_artery, [], 'all');

    X = linspace(v_min, v_max, v_max - v_min + 1);
    n = size(X, 2);
    histo = zeros(size(X, 2), N_frame);

    for frame_idx = 1:N_frame

        for xx = 1:Nx

            for yy = 1:Ny

                if (v_histo_artery(xx, yy, frame_idx) ~= 0)
                    i = find(X == v_histo_artery(xx, yy, frame_idx));
                    histo(i, frame_idx) = histo(i, frame_idx) + 1;
                end

            end

        end

    end

    figure(156)
    yAx = [v_min v_max];
    xAx = [0 n * ToolBox.stride / (1000 * ToolBox.fs)];
    imagesc(xAx, yAx, histo)
    set(gca, 'YDir', 'normal')
    colormap('hot')
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    title("Velocity distribution in arteries")

    v_histo_veins = round(v_RMS_all .* maskVein);
    v_min = min(v_histo_veins, [], 'all');
    v_max = max(v_histo_veins, [], 'all');

    X = linspace(v_min, v_max, v_max - v_min + 1);
    n = size(X, 2);
    histo = zeros(size(X, 2), N_frame);

    for frame_idx = 1:N_frame

        for xx = 1:Nx

            for yy = 1:Ny

                if (v_histo_veins(xx, yy, frame_idx) ~= 0)
                    i = find(X == v_histo_veins(xx, yy, frame_idx));
                    histo(i, frame_idx) = histo(i, frame_idx) + 1;
                end

            end

        end

    end

    if veins_analysis
        figure(157)
        yAx = [v_min v_max];
        xAx = [0 n * ToolBox.stride / (1000 * ToolBox.fs)];
        imagesc(xAx, yAx, histo)
        set(gca, 'YDir', 'normal')
        colormap('hot')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title("Velocity distribution in veins")
    end

    % h = histogram(v_RMS(:,:,ii).*maskArtery);
    % X = h.BinCounts;
    %
    % [x,y] = max(X);
    % X(y) = 0;
    % X = smoothdata(X);
    % figure
    % plot(X)

    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_flowVideo_one_cycle')));
    open(w)

    for frame_idx = 1:N_frame
        writeVideo(w, squeeze(flowVideoRGB(:, :, :, frame_idx)));
    end

    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, '_flowVideo_one_cycle')), 'MPEG-4');
    open(w)

    for frame_idx = 1:N_frame
        writeVideo(w, squeeze(flowVideoRGB(:, :, :, frame_idx)));
    end

    close(w);

    try
        % Save colorbar
        colorfig = figure(116);
        colorfig.Units = 'normalized';
        colormap(cmap_artery)
        %hCB = colorbar('north');
        hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(Vmin_Arteries, 1)), string(round(Vmax_Arteries, 1))});
        set(gca, 'Visible', false)
        set(gca, 'LineWidth', 3);
        hCB.Position = [0.10 0.3 0.81 0.35];
        colorfig.Position(4) = 0.1000;
        fontsize(gca, 15, "points");
        print('-f116', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_colorbar_velocity_arteries.png')));

        if veins_analysis
            % Save colorbar
            colorfig = figure(117);
            colorfig.Units = 'normalized';
            colormap(cmap_vein)
            %hCB = colorbar('north');
            hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(Vmin_Veins, 1)), string(round(Vmax_Veins, 1))});
            set(gca, 'Visible', false)
            set(gca, 'LineWidth', 3);
            hCB.Position = [0.10 0.3 0.81 0.35];
            colorfig.Position(4) = 0.1000;
            fontsize(gca, 15, "points");
            print('-f117', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_colorbar_velocity_veins.png')));
        end

    catch
        disp('fail saving colorbars')
    end

    f158 = figure(158);
    f158.Position = [300, 300, 600, 630];
    gifWriter = GifWriter("Animated_Velocity_map", 0.04, ToolBox);

    for frame_idx = 1:N_frame
        imagesc(flowVideoRGB(:, :, :, frame_idx));
        title('Flow Map');
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 14, "points");
        % hCB = colorbar('southoutside', 'Ticks', [0, 1], 'TickLabels', {string(round(Vmin_Arteries, 1)), string(round(Vmax_Arteries, 1))});
        % hCB.Label.String = 'Velocity (mm.s^{-1})';
        % hCB.Label.FontSize = 12;
        % colormap(cmap_artery);

        gifWriter = gifWriter.write(flowVideoRGB(:, :, :, frame_idx));

    end

    gifWriter.generate();

    print('-f156', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_velocity_distribution_arteries.png')));

    if veins_analysis
        print('-f157', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_velocity_distribution_veins.png')));
    end

    imwrite(flowImageRGB, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_flow_image.png')));

    %% Init of histogram axis

    PW_params = Parameters_json(path);
    veins_analysis = PW_params.veins_analysis;

    %%
    [maskSection] = create_mask_section(ImgM0, maskArtery, ToolBox, path);
    maskArtery_section = maskArtery & maskSection;

    %or

    % maskArtery_section = maskArtery;
    %%

    v_histo_artery = round(v_RMS_all .* maskArtery_section);
    v_min_artery = min(v_histo_artery, [], 'all');
    v_max_artery = max(v_histo_artery, [], 'all');

    if veins_analysis
        v_histo_vein = round(v_RMS_all .* maskVein);
        v_min_vein = min(v_histo_vein, [], 'all');
        v_max_vein = max(v_histo_vein, [], 'all');
        v_max_all = max(v_max_artery, v_max_vein);
        v_min_all = min(v_min_artery, v_min_vein);
    else
        v_max_all = v_max_artery;
        v_min_all = v_min_artery;
    end

    v_max_all_display = round(0.8 * v_max_all);
    v_min_all_display = round(0.8 * v_min_all);

    yAx = [v_min_all v_max_all];
    %yAx_display = [-20  80] ;
    yAx_display = yAx;
    %FIXME trouver un moyen de croper proprement sans décaler le zéro

    %% Velocity Histogram in arteries
    %FIXME prctile 10% Y = percentil(X,[5 95])

    X = linspace(v_min_all, v_max_all, v_max_all - v_min_all + 1);
    % n = size(X, 2);
    xAx = [0 N_frame * ToolBox.stride / (1000 * ToolBox.fs)];
    histo_artery = zeros(size(X, 2), N_frame);
    %histo_video_artery = zeros(size(X,2),size(dataCubeM2M0,3),3,size(dataCubeM2M0,3));

    f_distrib_artery = figure(157);
    f_distrib_artery.Position(3:4) = [600 275];
    index_min = find(X == v_min_all_display);
    index_max = find(X == v_max_all_display);
    imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
    set(gca, 'YDir', 'normal')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    colormap("hot")
    f = getframe(gcf);
    [M, N, ~] = size(f.cdata);
    histo_video_artery = zeros(M, N, 3, N_frame);

    %FIXME avoir une ligne à zéro de trois pixel
    %FIXME getframe pour la couleur

    gifWriter = GifWriter("Animated_Velocity_Histogram_Artery", 0.04, ToolBox);

    for frame_idx = 1:N_frame

        for xx = 1:Nx

            for yy = 1:Ny

                if (v_histo_artery(xx, yy, frame_idx) ~= 0)
                    i = find(X == v_histo_artery(xx, yy, frame_idx));
                    histo_artery(i, frame_idx) = histo_artery(i, frame_idx) + 1;
                end

            end

        end

        %histo_video_artery(:,:,t) = flip(histo_artery,1);
        figure(157)
        imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
        set(gca, 'YDir', 'normal')
        title("Velocity distribution in arteries")
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        f = getframe(gcf);
        histo_video_artery(:, :, :, frame_idx) = imresize(f.cdata, [M N]);
        gifWriter = gifWriter.write(f);

    end

    gifWriter.generate();
    velocity_dist_arteries = frame2im(getframe(gca));

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_velocity_histogram_artery.avi')));
    tmp = mat2gray(histo_video_artery);
    open(w)

    for frame_idx = 1:N_frame
        writeVideo(w, tmp(:, :, :, frame_idx));
    end

    close(w);

    %% Velocity Histogram in veins

    if veins_analysis

        X = linspace(v_min_all, v_max_all, v_max_all - v_min_all + 1);
        % n = size(X, 2);
        xAx = [0 N_frame * ToolBox.stride / (1000 * ToolBox.fs)];
        histo_vein = zeros(size(X, 2), N_frame);
        %histo_video_vein = zeros(size(X,2),size(dataCubeM2M0,3),size(dataCubeM2M0,3));
        f = getframe(gcf);
        [M, N, ~] = size(f.cdata);
        histo_video_vein = zeros(M, N, 3, N_frame);

        f_distrib_vein = figure(158);
        f_distrib_vein.Position(3:4) = [600 275];
        imagesc(xAx, yAx_display, histo_vein(index_min:index_max, :))
        set(gca, 'YDir', 'reverse')
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        colormap("bone")

        gifWriter = GifWriter("Animated_Velocity_Histogram_Veins", 0.04, ToolBox);

        for frame_idx = 1:N_frame

            for xx = 1:Nx

                for yy = 1:Ny

                    if (v_histo_vein(xx, yy, frame_idx) ~= 0)
                        i = find(X == v_histo_vein(xx, yy, frame_idx));
                        histo_vein(i, frame_idx) = histo_vein(i, frame_idx) + 1;
                    end

                end

            end

            %histo_video_vein(:,:,t) = flip(histo_vein,1);
            figure(158)
            imagesc(xAx, yAx_display, histo_vein(index_min:index_max, :))
            set(gca, 'YDir', 'normal')
            ylabel('Velocity (mm.s^{-1})')
            xlabel('Time (s)')
            title("Velocity distribution in veins")
            set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
            f = getframe(gcf);
            %histo_video_vein(:,:,:,t) = imresize(f.cdata,[M N]);
            histo_video_vein(:, :, :, frame_idx) = f.cdata;
            gifWriter = gifWriter.write(f);
        end

        velocity_dist_veins = frame2im(getframe(gca));
        gifWriter.generate();

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_velocity_histogram_veins.avi')));
        tmp = mat2gray(histo_video_vein);
        open(w)

        for frame_idx = 1:N_frame
            writeVideo(w, tmp(:, :, :, frame_idx));
        end

        close(w);

    end

    if veins_analysis
        freq_video(:, :, 1, :) = imresize3(squeeze(FreqVideoRGB(:, :, 1, :)), [550 550 N_frame]);
        freq_video(:, :, 2, :) = imresize3(squeeze(FreqVideoRGB(:, :, 2, :)), [550 550 N_frame]);
        freq_video(:, :, 3, :) = imresize3(squeeze(FreqVideoRGB(:, :, 3, :)), [550 550 N_frame]);
        combinedGifs = cat(2, freq_video, cat(1, mat2gray(histo_video_artery), mat2gray(histo_video_vein)));
        freq_video_mean = rescale(imresize3(imread(fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png'))), [550 550 3]));
        imwrite(cat(2, freq_video_mean, cat(1, mat2gray(histo_video_artery(:, :, :, end)), mat2gray(histo_video_vein(:, :, :, end)))), fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        freq_video(:, :, 1, :) = imresize3(squeeze(FreqVideoRGB(:, :, 1, :)), [600 600 N_frame]);
        freq_video(:, :, 2, :) = imresize3(squeeze(FreqVideoRGB(:, :, 2, :)), [600 600 N_frame]);
        freq_video(:, :, 3, :) = imresize3(squeeze(FreqVideoRGB(:, :, 3, :)), [600 600 N_frame]);
        combinedGifs = cat(1, freq_video, mat2gray(histo_video_artery));
        freq_video_mean = rescale(imresize3(imread(fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png'))), [600 600 3]));
        imwrite(cat(1, freq_video_mean, mat2gray(histo_video_artery(:, :, :, end))), fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

    gifWriter = GifWriter("Animated_Velocity_Histogram_Combined", 0.04, ToolBox);

    for frame_idx = 1:N_frame

        gifWriter = gifWriter.write(combinedGifs(:, :, :, frame_idx));

    end

    gifWriter.generate();

    print('-f157', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_velocity_distribution_arteries_fullcycle.png')));
    imwrite(velocity_dist_arteries, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Velocity_distribution_in_arteries.png')), 'png');

    if veins_analysis
        print('-f158', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_velocity_distribution_veins_fullcycle.png')));
        imwrite(velocity_dist_veins, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Velocity_distribution_in_veins.png')), 'png');
    end

    close all

end
