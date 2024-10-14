function [] = ArterialResistivityIndex(v_RMS, flatfieldM0, maskArtery, ToolBox, path)

    PW_params = Parameters_json(path);
    exportGif = PW_params.exportGifs;

    disp('arterial resistivity and pulsatility...');

    mkdir(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex')
    mkdir(ToolBox.PW_path_eps, 'arterialResistivityPulsatilityIndex')

    [numX, numY, numFrames] = size(flatfieldM0);

    meanIm = rescale(mean(flatfieldM0, 3));
    flatfieldM0 = rescale(flatfieldM0);

    [ARI, ARImap] = construct_resistivity_index(v_RMS, maskArtery);
    [API, APImap] = construct_pulsatility_index(v_RMS, maskArtery);

    %% Arterial Resisitivity Index

    [hue_ARI, sat_ARI, val_ARI, cmap] = createARI_HSVmap(ARImap, ARI, meanIm, maskArtery, ToolBox);
    % arterial resistivity map RGB
    ARImapRGB = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
    ARImapRGB = ARImapRGB .* maskArtery + ones(numX, numY, 3) .* meanIm .* ~maskArtery;

    ARIvideoRGB = zeros(numX, numY, 3, numFrames);

    for frameIdx = 1:numFrames
        [hue_ARI, sat_ARI, val_ARI] = createARI_HSVmap(ARImap, ARI, flatfieldM0(:, :, frameIdx), maskArtery, ToolBox);
        %sat_ARI = sat_ARI.*(val_ARI.*maskArtery);
        ARIvideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
        img_M0 = flatfieldM0(:, :, frameIdx);
        ARIvideoRGB(:, :, :, frameIdx) = ARIvideoRGB(:, :, :, frameIdx) .* maskArtery + ones(numX, numY, 3, 1) .* img_M0 .* ~maskArtery;
    end

    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ARIVideo')));
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, '_ARIVideo')), 'MPEG-4');
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);

    %% ARTERIAL RESISTIVITY FIGURES
    figure(70)
    imagesc(ARImapRGB);
    title(strcat('Arterial Resistivity Index avg. : ', sprintf(" %3.2f", ARI)));
    axis image
    axis off
    set(gca, 'LineWidth', 2);
    fontsize(gca, 14, "points");
    c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
    c.Label.String = 'Arterial resistivity index';
    c.Label.FontSize = 14;
    colormap(cmap);
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.eps')))

    imwrite(ARImapRGB, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARIMapColored.png')), 'png')
    imwrite(ARImap, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARIMapRaw.png')), 'png')

    % Save colorbar
    colorfig = figure(113);
    colorfig.Units = 'normalized';
    colormap(cmap)
    hCB = colorbar('north');
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 14, "points");
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapColorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapColorbar.eps')))

    if exportGif
        f71 = figure(71);
        f71.Position = [300, 300, 570, 630];

        timePeriod = ToolBox.stride / ToolBox.fs / 1000;
        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "ArterialResistivityIndex")), timePeriod, 0.04, numFrames);

        for frameIdx = 1:numFrames
            imagesc(ARIvideoRGB(:, :, :, frameIdx));
            title(strcat('Arterial resistivity index value : ', sprintf(" %3.2f", ARI)));
            axis image
            axis off
            set(gca, 'LineWidth', 2);
            fontsize(gca, 12, "points");
            c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
            c.Label.String = 'Arterial resistivity index';
            c.Label.FontSize = 12;

            colormap(cmap);

            frame = getframe(f71, [40 10 500 600]);
            gifWriter.write(frame, frameIdx);

        end

        gifWriter.generate();
        gifWriter.delete();
    end

    %% Arterial Pulsatility Index

    [hue_API, sat_API, val_API, cmap] = createARI_HSVmap(APImap, API, meanIm, maskArtery, ToolBox);
    % arterial resistivity map RGB
    APImapRGB = hsv2rgb(hue_API, sat_API, val_API);
    APImapRGB = APImapRGB .* maskArtery + ones(numX, numY, 3) .* meanIm .* ~maskArtery;

    APIvideoRGB = zeros(numX, numY, 3, numFrames);

    for frameIdx = 1:numFrames
        [hue_API, sat_API, val_API] = createARI_HSVmap(APImap, API, flatfieldM0(:, :, frameIdx), maskArtery, ToolBox);
        APIvideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_API, sat_API, val_API);
        img_M0 = flatfieldM0(:, :, frameIdx);
        APIvideoRGB(:, :, :, frameIdx) = APIvideoRGB(:, :, :, frameIdx) .* maskArtery + ones(numX, numY, 3, 1) .* img_M0 .* ~maskArtery;
    end
    APIvideoRGB (APIvideoRGB < 0) = 0;
    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_APIVideo')));
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, squeeze(APIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, '_APIVideo')), 'MPEG-4');
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, squeeze(APIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);

    %% ARTERIAL PULSATILITY FIGURES

    figure(72)
    imagesc(APImapRGB);
    title(strcat('Arterial Pulsatility Index avg. : ', sprintf(" %3.2f", API)));
    axis image
    axis off
    set(gca, 'LineWidth', 2);
    fontsize(gca, 14, "points");
    c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
    c.Label.String = 'Arterial pulsatility index';
    c.Label.FontSize = 14;
    colormap(cmap);
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APImapFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APImapFig.eps')))

    imwrite(APImapRGB, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APIMapColored.png')), 'png')
    imwrite(APImap, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APIMapRaw.png')), 'png')

    % Save colorbar
    colorfig = figure(113);
    colorfig.Units = 'normalized';
    colormap(cmap)
    hCB = colorbar('north');
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 14, "points");
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APImapColorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityPulsatilityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'APImapColorbar.eps')))

    f73 = figure(73);
    f73.Position = [300, 300, 570, 630];

    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "ArterialPulsatilityIndex")), timePeriod, 0.04, numFrames);

    for frameIdx = 1:numFrames
        imagesc(APIvideoRGB(:, :, :, frameIdx));
        title(strcat('Arterial pulsatility index value : ', sprintf(" %3.2f", API)));
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 12, "points");
        c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
        c.Label.String = 'Arterial pulsatility index';
        c.Label.FontSize = 12;

        colormap(cmap);

        frame = getframe(f73, [40 10 500 600]);
        gifWriter.write(frame, frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();

    close all

end
