function [v_video_RGB, v_mean_RGB] = flowMap(v_video, maskSection, maskArtery, maskVein, M0_ff_video, xy_barycenter, TB)

    params = TB.getParams;
    veinsAnalysis = params.veins_analysis;
    exportVideos = params.exportVideos;

    % Rescale once
    M0_ff_video = rescale(M0_ff_video);
    M0_ff_image = rescale(mean(M0_ff_video, 3));
    [numX, numY, numFrames] = size(v_video);

    if veinsAnalysis
        % Precompute masks
        maskAV = maskArtery & maskVein;
        maskArtery = maskArtery & ~maskAV;
        maskVein = maskVein & ~maskAV;
    end

    % Precompute constants
    x_c = xy_barycenter(1) / numX;
    y_c = xy_barycenter(2) / numY;
    r1 = params.velocityBigRadiusRatio;
    r2 = params.velocitySmallRadiusRatio;
    w = 0.01;

    % Precompute circles
    circle1 = diskMask(numX, numY, r1, r1 + w, center = [x_c, y_c]);
    circle2 = diskMask(numX, numY, r2 - w, r2, center = [x_c, y_c]);
    circles = circle1 | circle2;

    % Precompute colormaps
    cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
    if veinsAnalysis
        cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
        cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [0 1 0], 2/3, [1 1 1], 1);
    end

    % Precompute v_rescaled
    v_max = max(v_video(maskSection));
    v_min = 0;
    v_rescaled = abs(v_video);
    v_rescaled = (v_rescaled - v_min) / v_max;

    % Precompute v_mean and v_mean_rescaled
    v_mean = squeeze(mean(v_video, 3));
    v_mean_rescaled = squeeze(mean(v_rescaled, 3));
    v_video_RGB = zeros(numX, numY, 3, numFrames);

    if veinsAnalysis
        % Velocity Images
        velocityIm(v_mean, maskArtery, cmapArtery, 'arteries', colorbarOn = true);
        velocityIm(v_mean, maskVein, cmapVein, 'veins', colorbarOn = true);
        velocityIm(v_mean, maskArtery | maskVein, turbo, 'vessels', colorbarOn = true);

        % velocity Colorbars of the gif
        velocityColorbar(cmapArtery, v_min, v_max, 'Arteries');
        velocityColorbar(cmapVein, v_min, v_max, 'Veins');
        
        % Average Image
        v_mean_RGB = flowMapImg(v_mean_rescaled, {maskArtery, maskVein, maskAV}, {cmapArtery, cmapVein, cmapAV} ,background = M0_ff_image, circles = circles);
        imwrite(v_mean_RGB, fullfile(TB.path_png, 'bloodFlowVelocity', sprintf("%s_%s", TB.main_foldername, 'v_mean.png')))

        % Gif frame generation
        parfor frameIdx = 1:numFrames
            v_video_RGB(:, :, :, frameIdx) = flowMapImg(v_rescaled(:, :, frameIdx), {maskArtery, maskVein, maskAV}, {cmapArtery, cmapVein, cmapAV} ,background = M0_ff_video(:, :, frameIdx), circles = circles);
        end

    else

        % Velocity Images
        velocityIm(v_mean, maskArtery, cmapArtery, 'arteries', colorbarOn = true);

        % velocity Colorbars of the gif
        velocityColorbar(cmapArtery, v_min, v_max, 'Arteries');

        % Average Image
        v_mean_RGB = flowMapImg(v_mean_rescaled, {maskArtery}, {cmapArtery} ,background = M0_ff_image, circles = circles);
        imwrite(v_mean_RGB, fullfile(TB.path_png, 'bloodFlowVelocity', sprintf("%s_%s", TB.main_foldername, 'v_mean.png')))

        % Gif frame generation
        parfor frameIdx = 1:numFrames
            v_video_RGB(:, :, :, frameIdx) = flowMapImg(v_rescaled(:, :, frameIdx), {maskArtery}, {cmapArtery} ,background = M0_ff_video(:, :, frameIdx), circles = circles);
        end

    end

    if exportVideos
        writeGifOnDisc(v_video_RGB, "flowMap");
    end

end
