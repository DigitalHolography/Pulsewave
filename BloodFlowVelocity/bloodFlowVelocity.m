function [] = bloodFlowVelocity(v_video, maskArtery, maskVein, maskSection, M0_ff_video)

close all

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);
veinsAnalysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;

mkdir(ToolBox.PW_path_png, 'bloodFlowVelocity')
mkdir(ToolBox.PW_path_eps, 'bloodFlowVelocity')

tic

% TRUE MIN and MAX V but not realistic
M0_ff_video = rescale(M0_ff_video);
M0_ff_image = rescale(mean(M0_ff_video, 3));
[numX, numY, numFrames] = size(v_video);

% AV = Artery AND Vein
maskAV = maskArtery & maskVein;
maskArterySection = maskArtery & maskSection & ~maskAV;
maskVeinSection = maskVein & maskSection & ~maskAV;

cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

%% 1) VELOCITY VIDEO
tVelocityVideo = tic;

v_video_RGB = zeros(numX, numY, 3, numFrames);

v_max = max(v_video(maskSection));
v_min = 0;

v_mean = squeeze(mean(v_video(:, :, :), 3));
v_rescaled = v_video .* (v_video > 0);
v_rescaled = (v_rescaled - v_min) / v_max;
v_mean_rescaled = squeeze(mean(v_rescaled(:, :, :), 3));

velocityIm(v_mean, maskArtery, cmapArtery, 'arteries', colorbarOn = true);
velocityColorbar(cmapArtery, v_min, v_max, 'Arteries');

v_mean_Artery = setcmap(v_mean_rescaled, maskArtery, cmapArtery);

if veinsAnalysis

    velocityIm(v_mean, (maskArtery | maskVein), turbo, 'vessels', colorbarOn = true);

    velocityIm(v_mean, maskVein, cmapVein, 'veins', colorbarOn = true);
    velocityColorbar(cmapVein, v_min, v_max, 'Veins');

    v_mean_RGB = v_mean_Artery .* maskArtery + M0_ff_image .* ~maskArtery;
    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_artery_mean.png')))

    v_mean_Vein = setcmap(v_mean_rescaled, maskVein, cmapVein);
    v_mean_AV = setcmap(v_mean_rescaled, (maskArtery&maskVein), cmapAV);
    v_mean_RGB = (v_mean_Artery + v_mean_Vein) .* ~maskAV + v_mean_AV + M0_ff_image .* ~(maskArtery | maskVein);

    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_mean.png')))

    parfor ii = 1:numFrames
        v_frame_Artery = setcmap(v_rescaled(:, :, ii), maskArtery, cmapArtery);
        v_frame_Vein = setcmap(v_rescaled(:, :, ii), maskVein, cmapVein);
        v_mean_AV = setcmap(v_rescaled(:, :, ii), (maskArtery&maskVein), cmapAV);
        v_video_RGB(:, :, :, ii) = (v_frame_Artery + v_frame_Vein) .* ~maskAV + v_mean_AV + M0_ff_video(:, :, ii) .* ~(maskArtery | maskVein);
    end

else
    v_mean_RGB = v_mean_Artery .* maskArtery + M0_ff_image .* ~maskArtery;
    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_mean.png')))

    parfor ii = 1:numFrames
        v_frame_Artery = setcmap(v_rescaled(:, :, ii), maskArtery, cmapArtery);
        v_video_RGB(:, :, :, ii) = v_frame_Artery + M0_ff_video(:, :, ii) .* ~maskArtery;
    end

end

if exportVideos
    writeGifOnDisc(v_video_RGB, "flowMap");

    % avi
    parfeval(backgroundPool, @writeVideoOnDisc, 0, v_video_RGB, fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')));

    % mp4
    parfeval(backgroundPool, @writeVideoOnDisc, 0, v_video_RGB, fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')), 'MPEG-4');
end

fprintf("- Velocity Map Timing : %ds\n", round(toc(tVelocityVideo)))

%% 2) HISTOGRAM

histoVideoArtery = VelocityHistogram(v_video, maskArterySection, 'Arteries');
if veinsAnalysis
    histoVideoVein = VelocityHistogram(v_video, maskVeinSection, 'Veins');
end


%% 3) COMBINED



if veinsAnalysis
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [550 550 3]));
    imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :)), mat2gray(histoVideoVein(:, :, :)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
else
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [600 600 3]));
    imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
end

if exportVideos

    if veinsAnalysis
        v_video_RGB4Gif(:, :, 1, :) = imresize3(squeeze(v_video_RGB(:, :, 1, :)), [550 550 numFrames]);
        v_video_RGB4Gif(:, :, 2, :) = imresize3(squeeze(v_video_RGB(:, :, 2, :)), [550 550 numFrames]);
        v_video_RGB4Gif(:, :, 3, :) = imresize3(squeeze(v_video_RGB(:, :, 3, :)), [550 550 numFrames]);
        combinedGifs = cat(2, v_video_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
    else
        v_video_RGB4Gif(:, :, 1, :) = imresize3(squeeze(v_video_RGB(:, :, 1, :)), [600 600 numFrames]);
        v_video_RGB4Gif(:, :, 2, :) = imresize3(squeeze(v_video_RGB(:, :, 2, :)), [600 600 numFrames]);
        v_video_RGB4Gif(:, :, 3, :) = imresize3(squeeze(v_video_RGB(:, :, 3, :)), [600 600 numFrames]);
        combinedGifs = cat(1, v_video_RGB4Gif, mat2gray(histoVideoArtery));
    end

    gifWriter = GifWriter("velocityHistogramCombined", numFrames);

    for frameIdx = 1:numFrames
        gifWriter.write(combinedGifs(:, :, :, frameIdx), frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();
end

close all

end
