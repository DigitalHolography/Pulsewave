function [] = ArterialResistivityIndex(t, v_video, maskArtery, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps
cArtery = [255 22 18] / 255;

if size(v_video, 3) > 1
    arterial_signal = squeeze(sum(v_video .* maskArtery, [1 2])) / nnz(maskArtery);
else
    arterial_signal = v_video;
end

% Remove outliers
arterial_signal = filloutliers(arterial_signal, 'center', 'movmedian', 3); % Replace outliers with the median
% arterial_signal = filloutliers(arterial_signal, 'clip', 'ThresholdFactor', 3); % Alternative: Clip outliers

[vMin, idxMin] = min(arterial_signal);
[vMax, idxMax] = max(arterial_signal);
vMean = mean(arterial_signal);
vStd = std(arterial_signal);

ARI = (vMax - vMin) / vMax;
API = (vMax - vMin) / vMean;

% ARI Graph
graphSignal(sprintf('ARI_%s', name), folder, ...
    t, arterial_signal, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('ARI %s = %0.2f', name, ARI), Legend = {'Processed', 'Raw'}, ...
    yLines = [0, vMin, vMax], yLineLabels = {'', '', ''});

% API Graph
graphSignal(sprintf('API_%s', name), folder, ...
    t, arterial_signal, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('API %s = %0.2f', name, API), Legend = {'Processed', 'Raw'}, ...
    yLines = [0, vMin, vMean, vMax], yLineLabels = {'', '', '', ''});

% Save image

if size(v_video, 3) > 1 % if given a video, output the image of ARI / API

    % Remove outliers from the video data
    v_processed = filloutliers(v_video, 'center', 'movmedian', 3); % Replace outliers in the video data
    % v_processed = filloutliers(v_video, 'clip', 'ThresholdFactor', 3); % Alternative: Clip outliers

    % Compute min, max, and mean velocities
    % imgMin = min(v_processed, [], 3);
    % imgMax = max(v_processed, [], 3);
    imgMin = v_processed(:, :, idxMin);
    imgMax = v_processed(:, :, idxMax);
    imgMean = mean(v_processed, 3);

    % Compute the velocity difference within the mask
    dV = (imgMax - imgMin) .* maskArtery;

    % Avoid division by zero by masking out regions where imgMax is zero
    nonzeroMax = imgMax > 0; % Create a mask for non-zero max values
    imgARI = zeros(size(imgMax)); % Initialize ARI image
    imgARI(nonzeroMax & maskArtery) = dV(nonzeroMax & maskArtery) ./ imgMax(nonzeroMax & maskArtery); % Compute ARI only where valid
    imgARI(imgARI > 1) = 1;
    imgARI(imgARI < 0) = 0;

    % Compute the API (Arterial Pulsatility Index)
    nonzeroMean = imgMean > 0; % Create a mask for non-zero mean values
    imgAPI = zeros(size(imgMean)); % Initialize API image
    imgAPI(nonzeroMean & maskArtery) = dV(nonzeroMean & maskArtery) ./ imgMean(nonzeroMean & maskArtery); % Compute API only where valid
    imgAPI(imgAPI > 3) = 3;
    imgAPI(imgAPI < 0) = 0;

    % Generate colormap
    [cmapARI] = cmapLAB(256, [1 1 1], 0, [1 0 0], 1);

    % Create RGB images for visualization
    RGBARI = setcmap(imgARI, maskArtery, cmapARI) + rescale(mean(v_processed, 3)) .* ~maskArtery;
    RGBAPI = setcmap(imgAPI, maskArtery, cmapARI) + rescale(mean(v_processed, 3)) .* ~maskArtery;

    % Display and save the ARI image
    fig = figure(211354);
    imagesc(RGBARI), axis off, axis image;
    colorbar, colormap(cmapARI), clim([0 1]);
    title(sprintf('ARI %s = %0.2f', name, ARI));
    saveas(fig, fullfile(ToolBox.path_png, folder, strcat(ToolBox.main_foldername, '_', 'ARI', '_', name)), 'png');

    % Display and save the API image
    f = figure(211354 + 1);
    imagesc(RGBAPI), axis off, axis image;
    colorbar, colormap(cmapARI), clim([0 3]);
    title(sprintf('API %s = %0.2f', name, API));
    saveas(f, fullfile(ToolBox.path_png, folder, strcat(ToolBox.main_foldername, '_', 'API', '_', name)), 'png');

    % Close figures
    close(f), close(fig);

else

    % Save txt
    fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_main_outputs', '.txt')), 'a');

    if strcmp(name, 'velocity')
        fprintf(fileID, 'Mean Velocity artery : %f (mm/s) \r\n', vMean);
        fprintf(fileID, 'Std Velocity artery : %f (mm/s) \r\n', vStd);
        fprintf(fileID, 'Max Velocity artery : %f (mm/s) \r\n', vMax);
        fprintf(fileID, 'Min Velocity artery : %f (mm/s) \r\n', vMin);
    end

    fprintf(fileID, 'Arterial Resistivity Index (%s) : %f  \r\n', name, ARI);
    fprintf(fileID, 'Arterial Pulsatility Index (%s) : %f  \r\n', name, API);
    fclose(fileID);
end

end
