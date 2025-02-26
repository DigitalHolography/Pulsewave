function [] = bloodFlowVelocity(v_video, maskArtery, maskVein, maskSection, M0_ff_video, xy_barycenter)

close all
ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
veinsAnalysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;

mkdir(ToolBox.PW_path_png, 'bloodFlowVelocity');
mkdir(ToolBox.PW_path_eps, 'bloodFlowVelocity');

% Precompute masks
maskAV = maskArtery & maskVein;
maskArterySection = maskArtery & maskSection & ~maskAV;
maskVeinSection = maskVein & maskSection & ~maskAV;

% Rescale once
M0_ff_video = rescale(M0_ff_video);
[~, ~, numFrames] = size(v_video);

%% 1) VELOCITY VIDEO
tVelocityVideo = tic;

[v_video_RGB, v_mean_RGB] = flowMap(v_video, maskSection, maskArtery, maskVein, M0_ff_video, xy_barycenter, ToolBox);

fprintf("- Velocity Map Timing : %ds\n", round(toc(tVelocityVideo)))

%% 2) HISTOGRAM

histoVideoArtery = VelocityHistogram(v_video, maskArterySection, 'Arteries');
if veinsAnalysis
    histoVideoVein = VelocityHistogram(v_video, maskVeinSection, 'Veins');
end

%% 3) COMBINED

if exportVideos
    if veinsAnalysis
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [550 550 3]));
        imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :, end)), mat2gray(histoVideoVein(:, :, :, end)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [600 600 3]));
        imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :, end))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

    if veinsAnalysis
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [550 550 numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [550 550 numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [550 550 numFrames]));
        combinedGifs = cat(2, v_video_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
    else
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [600 600 numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [600 600 numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [600 600 numFrames]));
        combinedGifs = cat(1, v_video_RGB4Gif, mat2gray(histoVideoArtery));
    end

    gifWriter = GifWriter("velocityHistogramCombined", numFrames);

    for frameIdx = 1:numFrames
        gifWriter.write(combinedGifs(:, :, :, frameIdx), frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();

else
    if veinsAnalysis
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [550 550 3]));
        imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :)), mat2gray(histoVideoVein(:, :, :)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [600 600 3]));
        imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

end

close all

end
