function [] = bloodFlowVelocity(v_RMS_all, v_one_cycle, maskArtery, maskVein, videoM0, FreqVideoRGB, ToolBox, path)

    PW_params = Parameters_json(path);
    veins_analysis = PW_params.veins_analysis;
    mkdir(ToolBox.PW_path_png, 'bloodFlowVelocity')
    mkdir(ToolBox.PW_path_eps, 'bloodFlowVelocity')

    tic

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
        val = Im .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein - (val_artery + val_vein)./2.*(maskArtery&maskVein);
        hue = (hue_artery + hue_vein)  - (hue_artery + hue_vein)./2.*(maskArtery&maskVein);
        sat = (sat_artery + sat_vein)  - (sat_artery + sat_vein)./2.*(maskArtery&maskVein);
        flowImageRGB = hsv2rgb(hue, sat, val);
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
        val_mean = v_mean .* (~(maskArtery + maskVein)) + val_artery_mean .* maskArtery + val_vein_mean .* maskVein  - (val_artery_mean + val_vein_mean)./2.*(maskArtery&maskVein);
        hue_mean = (hue_artery_mean + hue_vein_mean)  - (hue_artery_mean + hue_vein_mean)./2.*(maskArtery&maskVein);
        sat_mean = (sat_artery_mean + sat_vein_mean)  - (sat_artery_mean + sat_vein_mean)./2.*(maskArtery&maskVein);
        flowVideoRGB_mean = hsv2rgb(hue_mean, sat_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery + maskVein - (maskArtery&maskVein)) + ones(Nx, Ny, 3) .* ~(maskArtery + maskVein) .* M0_norm_mean;
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "vRMSMean.png")))

        parfor frameIdx = 1:N_frame
            v = mat2gray(v_RMS_all(:, :, frameIdx));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            [hue_vein, sat_vein, val_vein, ~] = createHSVmap(v, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
            val = v .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein - (val_artery + val_vein)./2.*(maskArtery&maskVein);
            hue = hue_artery + hue_vein - (hue_artery + hue_vein)./2.*(maskArtery&maskVein);
            sat = sat_artery + sat_vein - (sat_artery + sat_vein)./2.*(maskArtery&maskVein);
            flowVideoRGB(:, :, :, frameIdx) = hsv2rgb(hue, sat, val);
            flowVideoRGB(:, :, :, frameIdx) = flowVideoRGB(:, :, :, frameIdx) .* (maskArtery | maskVein) + ones(Nx, Ny, 3) .* ~(maskArtery | maskVein) .* videoM0_norm(:, :, frameIdx);
        end

    else

        [hue_artery_mean, sat_artery_mean, val_artery_mean, cmap_artery] = createHSVmap(v_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val_mean = v_mean .* (~(maskArtery)) + val_artery_mean .* maskArtery;
        flowVideoRGB_mean = hsv2rgb(hue_artery_mean, sat_artery_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery) + ones(Nx, Ny, 3) .* ~(maskArtery) .* M0_norm_mean;
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "vRMSMean.png")))

        parfor frameIdx = 1:N_frame
            v = mat2gray(v_RMS_all(:, :, frameIdx));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            val = v .* (~(maskArtery)) + val_artery .* maskArtery;
            val(val<0)=0;
            val(val>1)=1;
            flowVideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_artery, sat_artery, val);
            flowVideoRGB(:, :, :, frameIdx) = rescale(flowVideoRGB(:, :, :, frameIdx) .* (maskArtery) + ones(Nx, Ny, 3) .* ~(maskArtery) .* videoM0_norm(:, :, frameIdx));
        end

    end

    clear hue_artery sat_artery val_artery hue_vein sat_vein val_vein

%% Histogram of One_Cycle

    N_interp = size(v_one_cycle, 3);
    v_histo_artery = round(v_one_cycle .* maskArtery);
    v_min = min(v_histo_artery, [], 'all');
    v_max = max(v_histo_artery, [], 'all');

    X = linspace(v_min, v_max, v_max - v_min + 1);
    n = size(X, 2);
    histo = zeros(size(X, 2), N_interp);

    for frameIdx = 1:N_interp

        for xx = 1:Nx

            for yy = 1:Ny

                if (v_histo_artery(xx, yy, frameIdx) ~= 0)
                    i = find(X == v_histo_artery(xx, yy, frameIdx));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
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
    title("Velocity histogram in arteries")

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesOneCycle.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesOneCycle.eps')))

    if veins_analysis
        v_histo_veins = round(v_one_cycle .* maskVein);
        v_min = min(v_histo_veins, [], 'all');
        v_max = max(v_histo_veins, [], 'all');

        X = linspace(v_min, v_max, v_max - v_min + 1);
        n = size(X, 2);
        histo = zeros(size(X, 2), N_interp);

        for frameIdx = 1:N_interp

            for xx = 1:Nx

                for yy = 1:Ny

                    if (v_histo_veins(xx, yy, frameIdx) ~= 0)
                        i = find(X == v_histo_veins(xx, yy, frameIdx));
                        histo(i, frameIdx) = histo(i, frameIdx) + 1;
                    end

                end

            end

        end

        figure(157)
        yAx = [v_min v_max];
        xAx = [0 n * ToolBox.stride / (1000 * ToolBox.fs)];
        imagesc(xAx, yAx, histo)
        set(gca, 'YDir', 'normal')
        colormap('hot')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title("Velocity histogram in veins")

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsOneCycle.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsOneCycle.eps')))

    end

    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "vRMS.avi")));
    
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, squeeze(flowVideoRGB(:, :, :, frameIdx)));
    end

    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "vRMS.mp4")), 'MPEG-4');
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, squeeze(flowVideoRGB(:, :, :, frameIdx)));
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

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.eps')))

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

            exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityVeins.png')))
            exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityVeins.eps')))

        end

    catch
        disp('fail saving colorbars')
    end

    f158 = figure(158);
    f158.Position = [300, 300, 600, 630];
    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "vRMS")), timePeriod, 0.04, N_frame);

    for frameIdx = 1:N_frame
        imagesc(flowVideoRGB(:, :, :, frameIdx));
        title('Flow Map');
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 14, "points");
        % hCB = colorbar('southoutside', 'Ticks', [0, 1], 'TickLabels', {string(round(Vmin_Arteries, 1)), string(round(Vmax_Arteries, 1))});
        % hCB.Label.String = 'Velocity (mm.s^{-1})';
        % hCB.Label.FontSize = 12;
        % colormap(cmap_artery);

        gifWriter.write(flowVideoRGB(:, :, :, frameIdx), frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();
    
    fprintf("Velocity Map Timing :\n")
    toc

    %% Init of histogram axis

    %%
    [N, M] = size(maskArtery);
    radius1 = PW_params.velocity_bigRadiusRatio * (M + N) / 2;
    radius2 = PW_params.velocity_smallRadiusRatio * (M + N) / 2;    
    [maskSection] = createMaskSection(ImgM0, maskArtery,radius1,radius2,'_mask_artery_section_velocity_rgb.png', ToolBox, path);
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

    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityArtery")), timePeriod, 0.04, N_frame);

    for frameIdx = 1:N_frame

        for xx = 1:Nx

            for yy = 1:Ny

                if (v_histo_artery(xx, yy, frameIdx) ~= 0)
                    i = find(X == v_histo_artery(xx, yy, frameIdx)); % find the velocity range index
                    histo_artery(i, frameIdx) = histo_artery(i, frameIdx) + 1;
                end

            end

        end

        %histo_video_artery(:,:,t) = flip(histo_artery,1);
        figure(157)
        imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
        set(gca, 'YDir', 'normal')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title("Velocity distribution in arteries")
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        f = getframe(gcf);
        histo_video_artery(:, :, :, frameIdx) = imresize(f.cdata, [M N]);
        gifWriter.write(f, frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.eps')))

    % AVI

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityArtery.avi")));

    tmp = mat2gray(histo_video_artery);
    open(w)
    for frameIdx = 1:N_frame
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);

    % MP4

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityArtery.mp4")), 'MPEG-4');
    tmp = mat2gray(histo_video_artery);
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, tmp(:, :, :, frameIdx));
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

        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityVeins")), timePeriod, 0.04, N_frame);

        for frameIdx = 1:N_frame

            for xx = 1:Nx

                for yy = 1:Ny

                    if (v_histo_vein(xx, yy, frameIdx) ~= 0)
                        i = find(X == v_histo_vein(xx, yy, frameIdx));
                        histo_vein(i, frameIdx) = histo_vein(i, frameIdx) + 1;
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
            histo_video_vein(:, :, :, frameIdx) = f.cdata;
            gifWriter.write(f, frameIdx);
        end

        gifWriter.generate();
        gifWriter.delete();

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.eps')))

        % AVI

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityVein.avi")));

        tmp = mat2gray(histo_video_vein);
        open(w)
        for frameIdx = 1:N_frame
            writeVideo(w, tmp(:, :, :, frameIdx));
        end

        close(w);

        % MP4

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityVein.mp4")), 'MPEG-4');
        open(w)

        for frameIdx = 1:N_frame
            writeVideo(w, tmp(:, :, :, frameIdx));
        end

        close(w);

    end

    if veins_analysis
        freq_video(:, :, 1, :) = imresize3(squeeze(FreqVideoRGB(:, :, 1, :)), [550 550 N_frame]);
        freq_video(:, :, 2, :) = imresize3(squeeze(FreqVideoRGB(:, :, 2, :)), [550 550 N_frame]);
        freq_video(:, :, 3, :) = imresize3(squeeze(FreqVideoRGB(:, :, 3, :)), [550 550 N_frame]);
        combinedGifs = cat(2, freq_video, cat(1, mat2gray(histo_video_artery), mat2gray(histo_video_vein)));
        freq_video_mean = rescale(imresize3(imread(fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png'))), [550 550 3]));
        imwrite(cat(2, freq_video_mean, cat(1, mat2gray(histo_video_artery(:, :, :, end)), mat2gray(histo_video_vein(:, :, :, end)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        freq_video(:, :, 1, :) = imresize3(squeeze(FreqVideoRGB(:, :, 1, :)), [600 600 N_frame]);
        freq_video(:, :, 2, :) = imresize3(squeeze(FreqVideoRGB(:, :, 2, :)), [600 600 N_frame]);
        freq_video(:, :, 3, :) = imresize3(squeeze(FreqVideoRGB(:, :, 3, :)), [600 600 N_frame]);
        combinedGifs = cat(1, freq_video, mat2gray(histo_video_artery));
        freq_video_mean = rescale(imresize3(imread(fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png'))), [600 600 3]));
        imwrite(cat(1, freq_video_mean, mat2gray(histo_video_artery(:, :, :, end))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "velocityHistogramCombined")), timePeriod, 0.04, N_frame);

    for frameIdx = 1:N_frame
        gifWriter.write(combinedGifs(:, :, :, frameIdx), frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();

    close all
    toc

    
end
