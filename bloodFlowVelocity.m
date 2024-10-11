function [] = bloodFlowVelocity(v_RMS_video, vOneCycle, maskArtery, maskVein, M0_disp_video, ToolBox, path)

PW_params = Parameters_json(path);
veinsAnalysis = PW_params.veins_analysis;

mkdir(ToolBox.PW_path_png, 'bloodFlowVelocity')
mkdir(ToolBox.PW_path_eps, 'bloodFlowVelocity')

tic

% TRUE MIN and MAX V_RMS but not realistic
Im = rescale(mean(M0_disp_video, 3));
[numX, numY, numFrames] = size(v_RMS_video);

vArtery = sum(v_RMS_video .* maskArtery, [1 2]) / nnz(maskArtery);
vVein = sum(v_RMS_video .* maskVein, [1 2]) / nnz(maskVein);
VmaxArteries = max(vArtery(:));
VmaxVeins = max(vVein(:));
VminArteries = min(vArtery(:));
VminVeins = min(vVein(:));

%% Construct Velocity video ~3min
v_RMS_videoRGB = zeros(numX, numY, 3, numFrames);
M0_disp_image = mat2gray(mean(M0_disp_video, 3));
v_RMS_mean = mat2gray(squeeze(mean(v_RMS_video(:, :, :), 3)));

if veinsAnalysis
    [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(v_RMS_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
    [hue_vein_mean, sat_vein_mean, val_vein_mean, ~] = createHSVmap(v_RMS_mean, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
    
    val_mean = v_RMS_mean .* (~(maskArtery + maskVein)) ...
             + val_artery_mean .* maskArtery + val_vein_mean .* maskVein - (val_artery_mean + val_vein_mean) .* (maskArtery & maskVein) ./ 2;
    hue_mean = hue_artery_mean .* maskArtery + hue_vein_mean .* maskVein - (hue_artery_mean + hue_vein_mean) .* (maskArtery & maskVein) ./ 2;
    sat_mean = sat_artery_mean .* maskArtery + sat_vein_mean .* maskVein - (sat_artery_mean + sat_vein_mean) .* (maskArtery & maskVein) ./ 2;
    
    v_RMS_meanRGB = hsv2rgb(hue_mean, sat_mean, val_mean);
    v_RMS_meanRGB = v_RMS_meanRGB .* (maskArtery | maskVein) + M0_disp_image .* ~(maskArtery | maskVein);

    imwrite(v_RMS_meanRGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_RMS_mean.png')))

    parfor ii = 1:numFrames
        v = mat2gray(squeeze(v_RMS_video(:, :, ii)));
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        [hue_vein, sat_vein, val_vein, ~] = createHSVmap(v, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
        
        val = v .* ~(maskArtery | maskVein) ...
            + val_artery .* maskArtery + val_vein .* maskVein - (val_artery + val_vein) .* (maskArtery & maskVein) ./ 2;
        hue = hue_artery .* maskArtery + hue_vein .* maskVein - (hue_artery + hue_vein) .* (maskArtery & maskVein) ./ 2;
        sat = sat_artery .* maskArtery + sat_vein .* maskVein - (sat_artery + sat_vein) .* (maskArtery & maskVein) ./ 2;
        
        v_RMS_videoRGB(:, :, :, ii) = hsv2rgb(hue, sat, val);
        v_RMS_videoRGB(:, :, :, ii) = v_RMS_videoRGB(:, :, :, ii) .* (maskArtery | maskVein) + rescale(M0_disp_video(:, :, ii)) .* ~(maskArtery | maskVein);
    end

else
    [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(f_RMS_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
    val_mean = f_RMS_mean .* ~maskArtery + val_artery_mean .* maskArtery;
    v_RMS_meanRGB = hsv2rgb(hue_artery_mean, sat_artery_mean, val_mean);
    v_RMS_meanRGB = v_RMS_meanRGB .* maskArtery + M0_disp_image .* ~maskArtery;
    imwrite(v_RMS_meanRGB, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'v_RMS_mean.png')))

    parfor ii = 1:numFrames
        v = mat2gray(squeeze(v_RMS_video(:, :, ii)));
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val = (v .* ~maskArtery) + val_artery .* maskArtery;
        v_RMS_videoRGB(:, :, :, ii) = hsv2rgb(hue_artery, sat_artery, val);
        v_RMS_videoRGB(:, :, :, ii) = v_RMS_videoRGB(:, :, :, ii) .* maskArtery + rescale(M0_disp_video(:, :, ii)) .* ~(maskArtery);
    end

end

clear hue_artery sat_artery val_artery hue_vein sat_vein val_vein

%% Saving video
% avi

parfeval(backgroundPool, @writeVideoOnDisc, 0, v_RMS_videoRGB, fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')));

timePeriod = ToolBox.stride / ToolBox.fs / 1000;

writeGifOnDisc(v_RMS_videoRGB, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "flowMap")), timePeriod);

% mp4
parfeval(backgroundPool, @writeVideoOnDisc, 0, v_RMS_videoRGB, fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')), 'MPEG-4');

%% Histogram of One_Cycle

nInterpFrames = size(vOneCycle, 3);
vHistoArtery = round(vOneCycle .* maskArtery);
vMin = min(vHistoArtery, [], 'all');
vMax = max(vHistoArtery, [], 'all');

X = linspace(vMin, vMax, vMax - vMin + 1);
n = size(X, 2);
histo = zeros(size(X, 2), nInterpFrames);

for frameIdx = 1:nInterpFrames

    for xx = 1:numX

        for yy = 1:numY

            if (vHistoArtery(xx, yy, frameIdx) ~= 0)
                i = find(X == vHistoArtery(xx, yy, frameIdx));
                histo(i, frameIdx) = histo(i, frameIdx) + 1;
            end

        end

    end

end

figure(156)
yAx = [vMin vMax];
xAx = [0 n * ToolBox.stride / (1000 * ToolBox.fs)];
imagesc(xAx, yAx, histo)
pbaspect([1.618 1 1]);
set(gca, 'YDir', 'normal')
colormap('hot')
ylabel('Velocity (mm.s^{-1})')
xlabel('Time (s)')
title("Velocity histogram in arteries")

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesOneCycle.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'histogramVelocityArteriesOneCycle.eps')))

if veinsAnalysis
    vHistoVeins = round(vOneCycle .* maskVein);
    vMin = min(vHistoVeins, [], 'all');
    vMax = max(vHistoVeins, [], 'all');

    X = linspace(vMin, vMax, vMax - vMin + 1);
    n = size(X, 2);
    histo = zeros(size(X, 2), nInterpFrames);

    for frameIdx = 1:nInterpFrames

        for xx = 1:numX

            for yy = 1:numY

                if (vHistoVeins(xx, yy, frameIdx) ~= 0)
                    i = find(X == vHistoVeins(xx, yy, frameIdx));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
                end

            end

        end

    end

    figure(157)
    yAx = [vMin vMax];
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

for frameIdx = 1:numFrames
    writeVideo(w, squeeze(v_RMS_videoRGB(:, :, :, frameIdx)));
end

close(w);
% mp4
w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, "vRMS.mp4")), 'MPEG-4');
open(w)

for frameIdx = 1:numFrames
    writeVideo(w, squeeze(v_RMS_videoRGB(:, :, :, frameIdx)));
end

close(w);

try
    % Save colorbar
    colorfig = figure(116);
    colorfig.Units = 'normalized';
    colormap(cmapArtery)
    %hCB = colorbar('north');
    hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(VminArteries, 1)), string(round(VmaxArteries, 1))});
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 15, "points");

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarVelocityArteries.eps')))

    if veinsAnalysis
        % Save colorbar
        colorfig = figure(117);
        colorfig.Units = 'normalized';
        colormap(cmapVein)
        %hCB = colorbar('north');
        hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(VminVeins, 1)), string(round(VmaxVeins, 1))});
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
gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "vRMS")), timePeriod, 0.04, numFrames);

for frameIdx = 1:numFrames
    imagesc(v_RMS_videoRGB(:, :, :, frameIdx));
    title('Flow Map');
    axis image
    axis off
    set(gca, 'LineWidth', 2);
    fontsize(gca, 14, "points");
    % hCB = colorbar('southoutside', 'Ticks', [0, 1], 'TickLabels', {string(round(VminArteries, 1)), string(round(VmaxArteries, 1))});
    % hCB.Label.String = 'Velocity (mm.s^{-1})';
    % hCB.Label.FontSize = 12;
    % colormap(cmapArtery);

    gifWriter.write(v_RMS_videoRGB(:, :, :, frameIdx), frameIdx);

end

gifWriter.generate();
gifWriter.delete();

fprintf("Velocity Map Timing :\n")
toc

%% Init of histogram axis

radius1 = PW_params.velocityBigRadiusRatio * (numY + numX) / 2;
radius2 = PW_params.velocitySmallRadiusRatio * (numY + numX) / 2;
[maskSection] = createMaskSection(Im, maskArtery, radius1, radius2, '_maskArterySection_rgb.png', ToolBox, path);
maskArterySection = maskArtery & maskSection;

%or

% maskArterySection = maskArtery;
%%

vHistoArtery = round(v_RMS_video .* maskArterySection);
vMinArtery = min(vHistoArtery, [], 'all');
vMaxArtery = max(vHistoArtery, [], 'all');

if veinsAnalysis
    vHistoVein = round(v_RMS_video .* maskVein);
    vMinVein = min(vHistoVein, [], 'all');
    vMaxVein = max(vHistoVein, [], 'all');
    vMaxAll = max(vMaxArtery, vMaxVein);
    vMinAll = min(vMinArtery, vMinVein);
else
    vMaxAll = vMaxArtery;
    vMinAll = vMinArtery;
end

vMaxAllDisplay = round(0.8 * vMaxAll);
vMinAllDisplay = round(0.8 * vMinAll);

yAx = [vMinAll vMaxAll];
%yAxDisplay = [-20  80] ;
yAxDisplay = yAx;
%FIXME trouver un moyen de croper proprement sans décaler le zéro

%% Velocity Histogram in arteries
%FIXME prctile 10% Y = percentil(X,[5 95])

X = linspace(vMinAll, vMaxAll, vMaxAll - vMinAll + 1);
% n = size(X, 2);
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
histoArtery = zeros(size(X, 2), numFrames);
%histoVideoArtery = zeros(size(X,2),size(dataCubeM2M0,3),3,size(dataCubeM2M0,3));

fDistribArtery = figure(157);
fDistribArtery.Position(3:4) = [600 275];
indexMin = find(X == vMinAllDisplay);
indexMax = find(X == vMaxAllDisplay);
imagesc(xAx, yAxDisplay, histoArtery(indexMin:indexMax, :))
set(gca, 'YDir', 'normal')
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
colormap("hot")
f = getframe(gcf);
[numX_fig, numY_fig, ~] = size(f.cdata);
histoVideoArtery = zeros(numX_fig, numY_fig, 3, numFrames);

%FIXME avoir une ligne à zéro de trois pixel
%FIXME getframe pour la couleur

gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityArtery")), timePeriod, 0.04, numFrames);

for frameIdx = 1:numFrames

    for xx = 1:numX_fig

        for yy = 1:numY_fig

            if (vHistoArtery(xx, yy, frameIdx) ~= 0)
                i = find(X == vHistoArtery(xx, yy, frameIdx));
                histoArtery(i, frameIdx) = histoArtery(i, frameIdx) + 1;

            end

        end

    end

    %histoVideoArtery(:,:,t) = flip(histoArtery,1);
    figure(157)
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

%% Velocity Histogram in veins

if veinsAnalysis

    X = linspace(vMinAll, vMaxAll, vMaxAll - vMinAll + 1);
    % n = size(X, 2);
    xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
    histoVein = zeros(size(X, 2), numFrames);
    %histoVideoVein = zeros(size(X,2),size(dataCubeM2M0,3),size(dataCubeM2M0,3));
    f = getframe(gcf);
    [numX_fig, numY_fig, ~] = size(f.cdata);
    histoVideoVein = zeros(numX_fig, numY_fig, 3, numFrames);

    fDistribVein = figure(158);
    fDistribVein.Position(3:4) = [600 275];
    imagesc(xAx, yAxDisplay, histoVein(indexMin:indexMax, :))
    set(gca, 'YDir', 'reverse')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    colormap("bone")

    gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "histogramVelocityVeins")), timePeriod, 0.04, numFrames);

    for frameIdx = 1:numFrames

        for xx = 1:numX_fig

            for yy = 1:numY_fig

                if (vHistoVein(xx, yy, frameIdx) ~= 0)
                    i = find(X == vHistoVein(xx, yy, frameIdx));
                    histoVein(i, frameIdx) = histoVein(i, frameIdx) + 1;
                end

            end

        end

        %histoVideoVein(:,:,t) = flip(histoVein,1);
        figure(158)
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

end

if veinsAnalysis
    v_RMS_videoRGB4Gif(:, :, 1, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 1, :)), [550 550 numFrames]);
    v_RMS_videoRGB4Gif(:, :, 2, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 2, :)), [550 550 numFrames]);
    v_RMS_videoRGB4Gif(:, :, 3, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 3, :)), [550 550 numFrames]);
    combinedGifs = cat(2, v_RMS_videoRGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
    v_RMS_meanRGB4Gif = rescale(imresize3(v_RMS_meanRGB, [550 550 3]));
    imwrite(cat(2, v_RMS_meanRGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :, end)), mat2gray(histoVideoVein(:, :, :, end)))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
else
    v_RMS_videoRGB4Gif(:, :, 1, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 1, :)), [600 600 numFrames]);
    v_RMS_videoRGB4Gif(:, :, 2, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 2, :)), [600 600 numFrames]);
    v_RMS_videoRGB4Gif(:, :, 3, :) = imresize3(squeeze(v_RMS_videoRGB(:, :, 3, :)), [600 600 numFrames]);
    combinedGifs = cat(1, v_RMS_videoRGB4Gif, mat2gray(histoVideoArtery));
    v_RMS_meanRGB4Gif = rescale(imresize3(v_RMS_meanRGB, [600 600 3]));
    imwrite(cat(1, v_RMS_meanRGB4Gif, mat2gray(histoVideoArtery(:, :, :, end))), fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
end

gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "velocityHistogramCombined")), timePeriod, 0.04, numFrames);

for frameIdx = 1:numFrames
    gifWriter.write(combinedGifs(:, :, :, frameIdx), frameIdx);
end

gifWriter.generate();
gifWriter.delete();

close all
%% Velocity funnel Histogram in arteries (exactly the same but with an increasing number of points)
%FIXME prctile 10% Y = percentil(X,[5 95])

if PW_params.AllCirclesFlag
    radius0 = 0;
    radiusmid = (radius1 + radius2) / 2;
    radiusend = (M + N) / 2;
    deltar = radiusend / 100;
    deltarcentral = (radius1 - radius2) / 100; % two times radius1-radius2 in total

    Color_std = [0.7 0.7 0.7];

    v_RMS_frame = mean(v_RMS_all, 3);

    X = linspace(v_min_all, v_max_all, v_max_all - v_min_all + 1);
    histo_artery = zeros(size(X, 2), 100);

    for j = 1:100 % to parforize (change createMaskSection)

        r1 = radiusmid + j * deltarcentral;
        r2 = radiusmid - j * deltarcentral;

        if mod(j, 10) == 0 % save one on 10
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, sprintf('_mask_artery_section_velocity_rgb%d.png', j), ToolBox, path);
        else
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, '_mask_artery_section_velocity_rgb100.png', ToolBox, path);
        end

        maskArtery_section = maskArtery & maskSection;
        v_histo_artery = round(v_RMS_frame .* maskArtery_section);

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
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, sprintf('_mask_artery_section_velocity_rgb%d.png', j), ToolBox, path);
        else
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, '_mask_artery_section_velocity_rgb100.png', ToolBox, path);
        end

        maskArtery_section_only = maskArtery & maskSection;

        non_zero_points = find(maskArtery_section_only);
        mean_velocity_in_section_only(j) = mean(v_RMS_frame(non_zero_points), [1, 2]);
        std_velocity_in_section_only(j) = std(v_RMS_frame(non_zero_points));

        rad(j) = sum(r1);

    end

    plot_velocity_funnel = figure(164);

    %
    xAx = number_of_points;

    f_distrib_artery = figure(197);
    f_distrib_artery.Position(3:4) = [500 275];
    index_min = find(X == v_min_all_display);
    index_max = find(X == v_max_all_display);
    imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
    set(gca, 'YDir', 'normal')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    colormap("hot")

    ylabel('Velocity (mm.s^{-1})')
    xlabel('Number of points')
    title("Velocity distribution in a growing artery section")

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArteriesxNumpoints.png')))

    plot_velocity_in_sections = figure(164);

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