function [] = ArterialResistivityIndex(v_RMS, videoM0, maskArtery, ToolBox)

    disp('arterial resistivity...');

    mkdir(ToolBox.PW_path_png, 'arterialResistivityIndex')
    mkdir(ToolBox.PW_path_eps, 'arterialResistivityIndex')

    [Nx, Ny, N_frame] = size(videoM0);

    meanIm = rescale(mean(videoM0, 3));
    videoM0 = rescale(videoM0);

    [ARI, ARImap] = construct_resistivity_index(v_RMS, maskArtery);
    ARImap(isnan(ARImap)) = 0;

    if ARI > 1
        ARI = 1;
    end

    [hue_ARI, sat_ARI, val_ARI, cmap] = createARI_HSVmap(ARImap, ARI, meanIm, maskArtery, ToolBox);
    % arterial resistivity map RGB
    ARImapRGB = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
    ARImapRGB = ARImapRGB .* maskArtery + ones(Nx, Ny, 3) .* meanIm .* ~maskArtery;

    ARIvideoRGB = zeros(Nx, Ny, 3, N_frame);

    for frameIdx = 1:N_frame
        [hue_ARI, sat_ARI, val_ARI] = createARI_HSVmap(ARImap, ARI, videoM0(:, :, frameIdx), maskArtery, ToolBox);
        %sat_ARI = sat_ARI.*(val_ARI.*maskArtery);
        ARIvideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
        img_M0 = videoM0(:, :, frameIdx);
        ARIvideoRGB(:, :, :, frameIdx) = ARIvideoRGB(:, :, :, frameIdx) .* maskArtery + ones(Nx, Ny, 3, 1) .* img_M0 .* ~maskArtery;
    end

    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ARIVideo')));
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, '_ARIVideo')), 'MPEG-4');
    open(w)

    for frameIdx = 1:N_frame
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, frameIdx)));
    end

    close(w);

    % disp('arterial resistivity...');
    % [ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma, img_avg] = construct_resistivity_index(onePulseVideo, maskArtery,path);
    % ARImap = ARImap.*maskArtery;
    %
    % % avi
    % w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ARIvideoRGB.avi')));
    % open(w)
    % ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
    % for frame_idx = 1:N_frame % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    %     writeVideo(w,squeeze(ARIvideoRGB(:,:,:,frame_idx))) ;
    % end
    % close(w);
    %
    % % mp4
    % w = VideoWriter(fullfile(ToolBox.PW_path_mp4,strcat(ToolBox.main_foldername,'_ARIvideoRGB.mp4')),'MPEG-4');
    % open(w)
    % ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
    % for frame_idx = 1:N_frame % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    %     writeVideo(w,squeeze(ARIvideoRGB(:,:,:,frame_idx))) ;
    % end
    %
    % disp('done.');

    %% Display Figure
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
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.eps')))

    imwrite(ARImapRGB, fullfile(ToolBox.PW_path_png, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARIMapColored.png')), 'png')
    imwrite(ARImap, fullfile(ToolBox.PW_path_png, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARIMapRaw.png')), 'png')

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
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapColorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'arterialResistivityIndex', sprintf("%s_%s", ToolBox.main_foldername, 'ARImapColorbar.eps')))

    f71 = figure(71);
    f71.Position = [300, 300, 570, 630];

    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "ArterialResistivityIndex")), timePeriod, 0.04, N_frame);

    for frameIdx = 1:N_frame
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

    close all

end
