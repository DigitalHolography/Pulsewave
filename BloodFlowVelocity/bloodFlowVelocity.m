function [] = bloodFlowVelocity(v_video, maskArtery, maskVein, maskSection, M0_ff_video)

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

%% 1) VELOCITY VIDEO ~3min
tVelocityVideo = tic;

v_video_RGB = zeros(numX, numY, 3, numFrames);

v_max = max(v_video(maskSection));
v_min = min(v_video(maskSection));

v_mean = squeeze(mean(v_video(:, :, :), 3));
v_rescaled = (v_video - v_min) / v_max;
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
    v_mean_AV = setcmap(v_mean_rescaled/2, (maskArtery&maskVein), cmapAV);
    v_mean_RGB = (v_mean_Artery + v_mean_Vein) .* ~maskAV + v_mean_AV + M0_ff_image .* ~(maskArtery | maskVein);

    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_mean.png')))

    parfor ii = 1:numFrames
        v_frame_Artery = setcmap(v_rescaled(:, :, ii), maskArtery, cmapArtery);
        v_frame_Vein = setcmap(v_rescaled(:, :, ii), maskVein, cmapVein);
        v_mean_AV = setcmap(v_rescaled(:, :, ii) / 2, (maskArtery&maskVein), cmapAV);
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
%% Init of histogram axis

tic

v_histoArtery = round(v_video .* maskArterySection);
v_minArtery = min(v_histoArtery, [], 'all');
v_maxArtery = max(v_histoArtery, [], 'all');

if veinsAnalysis
    v_histoVein = round(v_video .* maskVeinSection);
    v_minVein = min(v_histoVein, [], 'all');
    v_maxVein = max(v_histoVein, [], 'all');
    v_maxAll = max(v_maxArtery, v_maxVein);
    v_minAll = min(v_minArtery, v_minVein);
else
    v_maxAll = v_maxArtery;
    v_minAll = v_minArtery;
end

v_maxAll_display = round(0.8 * v_maxAll);
v_minAll_display = round(0.8 * v_minAll);

yAx = [v_minAll v_maxAll];
%yAxDisplay = [-20  80] ;
yAxDisplay = yAx;
%FIXME trouver un moyen de croper proprement sans décaler le zéro

%% Velocity Histogram in arteries
%FIXME prctile 10% Y = percentil(X,[5 95])

X = linspace(v_minAll, v_maxAll, v_maxAll - v_minAll + 1);
n = size(X, 2);
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
histoArtery = zeros(n, numFrames);

fDistribArtery = figure(20);
fDistribArtery.Position(3:4) = [600 275];
indexMin = find(X == v_minAll_display);
indexMax = find(X == v_maxAll_display);
imagesc(xAx, yAxDisplay, histoArtery(indexMin:indexMax, :))
set(gca, 'YDir', 'normal')
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
colormap(cmapArtery)
f = getframe(gcf);
[numX_fig, numY_fig, ~] = size(f.cdata);
histoVideoArtery = zeros(numX_fig, numY_fig, 3, numFrames);

if exportVideos
    gifWriter = GifWriter("histogramVelocityArtery", numFrames);

    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if maskArterySection(xx, yy) ~= 0
                    i = find(X == v_histoArtery(xx, yy, frameIdx));
                    histoArtery(i, frameIdx) = histoArtery(i, frameIdx) + 1;
                end

            end

        end

        %histoVideoArtery(:,:,t) = flip(histoArtery,1);
        figure(20)
        imagesc(xAx, yAxDisplay, histoArtery(indexMin:indexMax, :))
        set(gca, 'YDir', 'normal')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title("Velocity distribution in arteries")
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        f = getframe(gcf);
        histoVideoArtery(:, :, :, frameIdx) = imresize(f.cdata, [numX_fig numY_fig]);
        gifWriter.write(f, frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.eps')))

    % AVI

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityArtery.avi")));

    tmp = mat2gray(histoVideoArtery);
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);

    % MP4

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityArtery.mp4")), 'MPEG-4');
    tmp = mat2gray(histoVideoArtery);
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);
else

    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if maskArterySection(xx, yy) ~= 0
                    i = find(X == v_histoArtery(xx, yy, frameIdx));
                    histoArtery(i, frameIdx) = histoArtery(i, frameIdx) + 1;
                end

            end

        end

    end

    figure(20)
    imagesc(xAx, yAxDisplay, histoArtery(indexMin:indexMax, :))
    set(gca, 'YDir', 'normal')
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    title("Velocity distribution in arteries")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesFullCycle.eps')))

end

%% Velocity Histogram in veins

if veinsAnalysis

    X = linspace(v_minAll, v_maxAll, v_maxAll - v_minAll + 1);
    n = size(X, 2);
    xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
    histoVein = zeros(n, numFrames);

    fDistribVein = figure(21);
    fDistribVein.Position(3:4) = [600 275];
    imagesc(xAx, yAxDisplay, histoVein(indexMin:indexMax, :))
    set(gca, 'YDir', 'reverse')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    colormap(cmapVein)
    f = getframe(gcf);
    [numX_fig, numY_fig, ~] = size(f.cdata);
    histoVideoVein = zeros(numX_fig, numY_fig, 3, numFrames);

    if exportVideos
        gifWriter = GifWriter("histogramVelocityVeins", numFrames);

        for frameIdx = 1:numFrames

            for xx = 1:numX

                for yy = 1:numY

                    if maskVeinSection(xx, yy) ~= 0
                        i = find(X == v_histoVein(xx, yy, frameIdx));
                        histoVein(i, frameIdx) = histoVein(i, frameIdx) + 1;
                    end

                end

            end

            %histoVideoVein(:,:,t) = flip(histoVein,1);
            figure(21)
            imagesc(xAx, yAxDisplay, histoVein(indexMin:indexMax, :))
            set(gca, 'YDir', 'normal')
            ylabel('Velocity (mm.s^{-1})')
            xlabel('Time (s)')
            title("Velocity distribution in veins")
            set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
            f = getframe(gcf);
            %histoVideoVein(:,:,:,t) = imresize(f.cdata,[M N]);
            histoVideoVein(:, :, :, frameIdx) = f.cdata;
            gifWriter.write(f, frameIdx);
        end

        gifWriter.generate();
        gifWriter.delete();

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.eps')))

        % AVI

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityVein.avi")));

        tmp = mat2gray(histoVideoVein);
        open(w)

        for frameIdx = 1:numFrames
            writeVideo(w, tmp(:, :, :, frameIdx));
        end

        close(w);

        % MP4

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "histogramVelocityVein.mp4")), 'MPEG-4');
        open(w)

        for frameIdx = 1:numFrames
            writeVideo(w, tmp(:, :, :, frameIdx));
        end

        close(w);
    else

        for frameIdx = 1:numFrames

            for xx = 1:numX

                for yy = 1:numY

                    if maskVeinSection(xx, yy) ~= 0
                        i = find(X == v_histoVein(xx, yy, frameIdx));
                        histoVein(i, frameIdx) = histoVein(i, frameIdx) + 1;
                    end

                end

            end

        end

        figure(21)
        imagesc(xAx, yAxDisplay, histoVein(indexMin:indexMax, :))
        set(gca, 'YDir', 'normal')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title("Velocity distribution in veins")
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityVeinsFullCycle.eps')))
    end

end

if exportVideos

    if veinsAnalysis
        v_video_RGB4Gif(:, :, 1, :) = imresize3(squeeze(v_video_RGB(:, :, 1, :)), [550 550 numFrames]);
        v_video_RGB4Gif(:, :, 2, :) = imresize3(squeeze(v_video_RGB(:, :, 2, :)), [550 550 numFrames]);
        v_video_RGB4Gif(:, :, 3, :) = imresize3(squeeze(v_video_RGB(:, :, 3, :)), [550 550 numFrames]);
        combinedGifs = cat(2, v_video_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [550 550 3]));
        imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :, end)), mat2gray(histoVideoVein(:, :, :, end)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        v_video_RGB4Gif(:, :, 1, :) = imresize3(squeeze(v_video_RGB(:, :, 1, :)), [600 600 numFrames]);
        v_video_RGB4Gif(:, :, 2, :) = imresize3(squeeze(v_video_RGB(:, :, 2, :)), [600 600 numFrames]);
        v_video_RGB4Gif(:, :, 3, :) = imresize3(squeeze(v_video_RGB(:, :, 3, :)), [600 600 numFrames]);
        combinedGifs = cat(1, v_video_RGB4Gif, mat2gray(histoVideoArtery));
        v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [600 600 3]));
        imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :, end))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
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
