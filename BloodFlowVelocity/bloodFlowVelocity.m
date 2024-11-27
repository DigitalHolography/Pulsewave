function [] = bloodFlowVelocity(v_video, maskArtery, maskVein, maskSection, M0_ff_video, xy_barycenter, ToolBox, path)

PW_params = Parameters_json(path);
veinsAnalysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;

mkdir(ToolBox.PW_path_png, 'bloodFlowVelocity')
mkdir(ToolBox.PW_path_eps, 'bloodFlowVelocity')

tic

% TRUE MIN and MAX V but not realistic
M0_ff_video = rescale(M0_ff_video);
M0_ff_image = mean(M0_ff_video, 3);
[numX, numY, numFrames] = size(v_video);

maskArtery = maskArtery & maskSection;
maskVein = maskVein & maskSection;

v_Artery = sum(v_video .* maskArtery, [1 2]) / nnz(maskArtery);
v_Vein = sum(v_video .* maskVein, [1 2]) / nnz(maskVein);

v_maxArteries = max(v_Artery(:));
v_maxVeins = max(v_Vein(:));
v_minArteries = min(v_Artery(:));
v_minVeins = min(v_Vein(:));

cmapArtery = cmapPerception('rocket');
cmapVein = cmapPerception('mako');

%% 1) VELOCITY VIDEO ~3min
v_video_RGB = zeros(numX, numY, 3, numFrames);
v_video_rescale = rescale(v_video);
v_mean = squeeze(mean(v_video_rescale(:, :, :), 3));

if veinsAnalysis
    % v Vessels
    figure
    colormap turbo
    imagesc(v_mean)
    c = colorbar;
    c.Label.String = 'mm/s';
    axis image; axis off;
    exportgraphics(gcf, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "v_imagesc_Mean.png")), 'BackgroundColor', 'none', 'ContentType', 'vector', 'Resolution', 300);

    % v Artery
    figure
    colormap(cmapArtery)
    imagesc(v_mean .* maskArtery)
    c = colorbar;
    c.Label.String = 'mm/s';
    axis image; axis off;
    exportgraphics(gcf, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "v_imagesc_Mean_arteries.png")), 'BackgroundColor', 'none', 'ContentType', 'vector', 'Resolution', 300);

    % v Vein
    figure
    colormap(cmapVein);
    imagesc(v_mean .* maskVein);
    c = colorbar;
    c.Label.String = 'mm/s';
    axis image; axis off;
    exportgraphics(gcf, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "v_imagesc_Mean_veins.png")), 'BackgroundColor', 'none', 'ContentType', 'vector', 'Resolution', 300);

    v_mean_Artery = setcmap(v_mean, maskArtery, cmapArtery);
    v_mean_RGB = v_mean_Artery .* maskArtery + M0_ff_image .* ~maskArtery;
    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_artery_mean.png')))

    v_mean_Artery = setcmap(v_mean, maskArtery, cmapArtery);
    v_mean_Vein = setcmap(v_mean, maskVein, cmapVein);
    v_mean_RGB = v_mean_Artery .* maskArtery + v_mean_Vein .* maskVein + M0_ff_image .* ~(maskArtery | maskVein) - (v_mean_Artery + v_mean_Vein)/2 .* (maskArtery & maskVein);
    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_mean.png')))

    parfor ii = 1:numFrames
        v_frame_Artery = setcmap(v_video_rescale(:, :, ii), maskArtery, cmapArtery);
        v_frame_Vein = setcmap(v_video_rescale(:, :, ii), maskVein, cmapVein);
        v_video_RGB(:, :, :, ii) = v_frame_Artery .* maskArtery + v_frame_Vein .* maskVein + M0_ff_video(:, :, ii) .* ~(maskArtery | maskVein) - (v_frame_Artery + v_frame_Vein)/2 .* (maskArtery & maskVein);
    end

else
    v_mean_Artery = setcmap(v_mean, maskArtery, cmapArtery);
    v_mean_RGB = v_mean_Artery .* maskArtery + M0_ff_image .* ~maskArtery;
    imwrite(v_mean_RGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_mean.png')))

    parfor ii = 1:numFrames
        v_frame_Artery = setcmap(v_video_rescale(:, :, ii), maskArtery, cmapArtery);
        v_video_RGB(:, :, :, ii) = v_frame_Artery .* maskArtery + M0_ff_video(:, :, ii) .* ~maskArtery;
    end

end

%% Saving video
% avi

parfeval(backgroundPool, @writeVideoOnDisc, 0, v_video_RGB, fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')));

if exportVideos
    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    writeGifOnDisc(v_video_RGB, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "flowMap")), timePeriod);
end

% mp4
parfeval(backgroundPool, @writeVideoOnDisc, 0, v_video_RGB, fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')), 'MPEG-4');

% save video
% avi
w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "v_video.avi")));

open(w)

for frameIdx = 1:numFrames
    writeVideo(w, squeeze(v_video_RGB(:, :, :, frameIdx)));
end

close(w);
% mp4
w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "v_video.mp4")), 'MPEG-4');
open(w)

for frameIdx = 1:numFrames
    writeVideo(w, squeeze(v_video_RGB(:, :, :, frameIdx)));
end

close(w);

% Save colorbar
colorfig = figure(11);
colorfig.Units = 'normalized';
colormap(cmapArtery)
%hCB = colorbar('north');
hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(v_minArteries, 1)), string(round(v_maxArteries, 1))});
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
hCB.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.eps')))

if veinsAnalysis
    % Save colorbar
    colorfig = figure(12);
    colorfig.Units = 'normalized';
    colormap(cmapVein)
    %hCB = colorbar('north');
    hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(v_minVeins, 1)), string(round(v_maxVeins, 1))});
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 15, "points");

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityVeins.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityVeins.eps')))

end

if exportVideos
    f3 = figure(13); % v Gif
    f3.Position = [300, 300, 600, 630];
    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "v_video")), timePeriod, 0.04, numFrames);

    for frameIdx = 1:numFrames
        imagesc(v_video_RGB(:, :, :, frameIdx));
        title('Flow Map');
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 14, "points");

        gifWriter.write(v_video_RGB(:, :, :, frameIdx), frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();
end

fprintf("- Velocity Map Timing : %ds\n", round(toc))

%% 2) HISTOGRAM
%% Init of histogram axis

tic

radius1 = PW_params.velocityBigRadiusRatio * (numY + numX) / 2;
radius2 = PW_params.velocitySmallRadiusRatio * (numY + numX) / 2;
maskArterySection = maskArtery & maskSection;
maskVeinSection = maskVein & maskSection;

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
    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityArtery")), timePeriod, 0.04, numFrames);

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
        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityVeins")), timePeriod, 0.04, numFrames);

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

    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "velocityHistogramCombined")), timePeriod, 0.04, numFrames);

    for frameIdx = 1:numFrames
        gifWriter.write(combinedGifs(:, :, :, frameIdx), frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();
end

close all

%% Velocity funnel Histogram in arteries (exactly the same but with an increasing number of points)
%FIXME prctile 10% Y = percentil(X,[5 95])

if PW_params.AllCirclesFlag
    radius0 = 0;
    radiusmid = (radius1 + radius2) / 2;
    radiusend = (numX + numY) / 2;
    deltar = radiusend / 100;
    deltarcentral = (radius1 - radius2) / 100; % two times radius1-radius2 in total

    Color_std = [0.7 0.7 0.7];

    v_frame = mean(v_video, 3);

    X = linspace(v_min_all, v_max_all, v_max_all - v_min_all + 1);
    histo_artery = zeros(size(X, 2), 100);

    for j = 1:100 % to parforize (change createMaskSection)

        r1 = radiusmid + j * deltarcentral;
        r2 = radiusmid - j * deltarcentral;

        if mod(j, 10) == 0 % save one on 10
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, xy_barycenter, sprintf('_mask_artery_section_velocity_rgb%d.png', j), ToolBox, path);
        else
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, xy_barycenter, '_mask_artery_section_velocity_rgb100.png', ToolBox, path);
        end

        maskArtery_section = maskArtery & maskSection;
        v_histo_artery = round(v_frame .* maskArtery_section);

        for xx = 1:numX

            for yy = 1:numY

                if (v_histo_artery(xx, yy) ~= 0)
                    i = find(X == v_histo_artery(xx, yy));
                    histo_artery(i, j) = histo_artery(i, j) + 1;
                end

            end

        end

        histo_artery(:, j) = histo_artery(:, j) / sum(maskArtery_section, [1, 2]);

        number_of_points(j) = sum(maskArtery_section, [1, 2]);

        r1 = radius0 + j * deltar;
        r2 = radius0 + (j - 1) * deltar;

        if mod(j, 10) == 0 % save one on 10
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, xy_barycenter, sprintf('_mask_artery_section_velocity_rgb%d.png', j), ToolBox, path);
        else
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, xy_barycenter, '_mask_artery_section_velocity_rgb100.png', ToolBox, path);
        end

        maskArtery_section_only = maskArtery & maskSection;

        non_zero_points = find(maskArtery_section_only);
        mean_velocity_in_section_only(j) = mean(v_frame(non_zero_points), [1, 2]);
        std_velocity_in_section_only(j) = std(v_frame(non_zero_points));

        rad(j) = sum(r1);

    end

    plot_velocity_funnel = figure(30);

    %
    xAx = number_of_points;

    f_distrib_artery = figure(31);
    f_distrib_artery.Position(3:4) = [500 275];
    index_min = find(X == v_min_all_display);
    index_max = find(X == v_max_all_display);
    imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
    set(gca, 'YDir', 'normal')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    colormap(cmapArtery)

    ylabel('Velocity (mm.s^{-1})')
    xlabel('Number of points')
    title("Velocity distribution in a growing artery section")

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArteriesxNumpoints.png')))

    plot_velocity_in_sections = figure(32);

    curve1 = mean_velocity_in_section_only + 0.5 * std_velocity_in_section_only;
    curve2 = mean_velocity_in_section_only - 0.5 * std_velocity_in_section_only;
    rad2 = [rad, fliplr(rad)];
    inBetween = [curve1, fliplr(curve2)];

    fill(rad2, inBetween, Color_std);
    hold on;
    plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(rad, mean_velocity_in_section_only, '-k', 'LineWidth', 2);
    axis tight;
    hold off

    ylabel('Velocity (mm.s^{-1})')
    xlabel('radius ratio')
    title("Velocity in arteries sections with the radius to center")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArteriesxradius.png')))
end

end
