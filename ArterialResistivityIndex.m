function [] = ArterialResistivityIndex(v_RMS_one_cycle, videoM0_from_holowaves, maskArtery, ToolBox)
    PW_params = Parameters_json(path);
    disp('arterial resistivity...');
    
    name_log = strcat(ToolBox.PW_folder_name, '_log.txt');
    path_file_log = fullfile(ToolBox.PW_path_log, name_log);
    
    fileID = fopen(path_file_log, 'a+');
    fprintf(fileID, 'Arterial resistivity... \r\n');
    fclose(fileID);
    
    meanIm = rescale(mean(videoM0_from_holowaves, 3));
    videoM0_from_holowaves = rescale(videoM0_from_holowaves);
    
    [ARI, ARImap] = construct_resistivity_index(v_RMS_one_cycle, maskArtery);
    ARImap(isnan(ARImap)) = 0;
    
    if ARI > 1
        ARI = 1;
    end
    
    [hue_ARI, sat_ARI, val_ARI, cmap] = createARI_HSVmap(ARImap, ARI, meanIm, maskArtery, ToolBox);
    % arterial resistivity map RGB
    ARImapRGB = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
    ARImapRGB = ARImapRGB .* maskArtery + ones(size(ARImapRGB)) .* meanIm .* ~maskArtery;
    
    ARIvideoRGB = zeros(size(v_RMS_one_cycle, 1), size(v_RMS_one_cycle, 2), 3, size(v_RMS_one_cycle, 3));
    
    for ii = 1:size(v_RMS_one_cycle, 3)
        [hue_ARI, sat_ARI, val_ARI] = createARI_HSVmap(ARImap, ARI, videoM0_from_holowaves(:, :, ii), maskArtery, ToolBox);
        %sat_ARI = sat_ARI.*(val_ARI.*maskArtery);
        ARIvideoRGB(:, :, :, ii) = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
        img_M0 = videoM0_from_holowaves(:, :, ii);
        ARIvideoRGB(:, :, :, ii) = ARIvideoRGB(:, :, :, ii) .* maskArtery + ones(size(ARIvideoRGB(:, :, :, ii))) .* img_M0 .* ~maskArtery;
    end
    
    % save video
    % avi
    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ARIVideo')));
    open(w)
    
    for jj = 1:size(ARIvideoRGB, 4)
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, jj)));
    end
    
    close(w);
    % mp4
    w = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, '_ARIVideo')), 'MPEG-4');
    open(w)
    
    for jj = 1:size(ARIvideoRGB, 4)
        writeVideo(w, squeeze(ARIvideoRGB(:, :, :, jj)));
    end
    
    close(w);
    
    disp('done.');
    fileID = fopen(path_file_log, 'a+');
    fprintf(fileID, 'Done. \r\n');
    fclose(fileID);
    
    % disp('arterial resistivity...');
    % [ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma, img_avg] = construct_resistivity_index(onePulseVideo, maskArtery,path);
    % ARImap = ARImap.*maskArtery;
    %
    % % avi
    % w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ARIvideoRGB.avi')));
    % open(w)
    % ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
    % for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    %     writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
    % end
    % close(w);
    %
    % % mp4
    % w = VideoWriter(fullfile(ToolBox.PW_path_mp4,strcat(ToolBox.main_foldername,'_ARIvideoRGB.mp4')),'MPEG-4');
    % open(w)
    % ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
    % for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    %     writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
    % end
    %
    % disp('done.');
    
    %% Display Figure
    figure(70)
    imagesc(ARImapRGB);
    title(strcat('Arterial resistivity. avg. index value : ', sprintf(" %3.2f",ARI)));
    axis image
    axis off
    set(gca, 'LineWidth', 2);
    fontsize(gca, 12, "points");
    c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
    c.Label.String = 'Arterial resistivity index';
    c.Label.FontSize = 12;
    
    ARI_x = linspace(0, 1, 256);
    ARI_h = (ToolBox.ARI_hue_max - ToolBox.ARI_hue_min) * sigmoid(ARI_x, ToolBox.ARI_inflexion_point_hue, ToolBox.ARI_slope_hue) + ToolBox.ARI_hue_min;
    ARI_s = ones(1, 256);
    ARI_v = (ToolBox.ARI_val_max - ToolBox.ARI_val_min) * sigmoid(ARI_x, ToolBox.ARI_inflexion_point_val, ToolBox.ARI_slope_val) + ToolBox.ARI_val_min;
    % cmap = squeeze(hsv2rgb(ARI_h, ARI_s, ARI_v));
    colormap(cmap);
    
    % Save colorbar
    colorfig = figure(113);
    colorfig.Units = 'normalized';
    colormap(cmap)
    hCB = colorbar('north');
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 15, "points");
    
    f71 = figure(71);
    f71.Position = [300, 300, 570, 630];
    
    gifWriter = GifWriter("Animated_ARI",0.04,ToolBox);
    
    for tt = 1:size(ARIvideoRGB,4)
        imagesc(ARIvideoRGB(:,:,:,tt));
        title(strcat('Arterial resistivity index value : ', sprintf(" %3.2f",ARI)));
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 12, "points");
        c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
        c.Label.String = 'Arterial resistivity index';
        c.Label.FontSize = 12;
        
        colormap(cmap);
        
        frame = getframe(f71,[40 10 500 600]);
        gifWriter = gifWriter.write(frame);
        
    end
    
    gifWriter.generate();
    
    %% Save Figures
    print('-f113', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_ARI_map_colorbar.png')));
    print('-f70', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_ARI_map.png')));
    
    print('-f70', '-depsc', fullfile(ToolBox.PW_path_eps, strcat(ToolBox.main_foldername, '_resistivityMap.eps')));
    
    imwrite(ARImapRGB, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_ARI_map_img.png')), 'png');
    
    %close all
    