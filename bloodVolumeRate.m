function [] = bloodVolumeRate(maskArtery, maskVein, v_RMS_all, M0_disp_video, ToolBox, k, path, flagBloodVelocityProfile)

PW_params = Parameters_json(path);

veins_analysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;
force_width = [];
if ~isempty(PW_params.forcewidth)
    force_width = PW_params.forcewidth;
end

mkdir(ToolBox.PW_path_png, 'volumeRate')
mkdir(ToolBox.PW_path_eps, 'volumeRate')

[numX, numY, numFrames] = size(v_RMS_all);
[X, Y] = meshgrid(1:numY, 1:numX);

v_RMS_AVG = mean(v_RMS_all, 3);
fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

%% 0) Change mask section

radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (numY + numX) / 2;
radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (numY + numX) / 2;
cercle_mask1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius1;
cercle_mask2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius2;

maskSection = xor(cercle_mask1, cercle_mask2);

figure(1)
imagesc(maskSection), axis image

%% 1) Find the locations of the sections
%% 1) 1) Arteries

maskSectionArtery = maskSection .* maskArtery;

figure(110)
imagesc(maskSectionArtery .* v_RMS_AVG), axis image

maskSectionArtery = bwlabel(maskSectionArtery);

figure(111)
imagesc(maskSectionArtery), axis image

numSectionsArtery = max(maskSectionArtery, [], 'all');
masksSectionsArtery = zeros(numX, numY, numSectionsArtery);

parfor sectionIdx = 1:numSectionsArtery

    masksSectionsArtery(:, :, sectionIdx) = (maskSectionArtery == sectionIdx);

end

SubImg_locs_artery = zeros(numSectionsArtery, 2);
SubImg_width_artery = zeros(numSectionsArtery, 1);

for sectionIdx = 1:numSectionsArtery
    [row, col] = find(masksSectionsArtery(:, :, sectionIdx));
    SubImg_locs_artery(sectionIdx, 1) = round(mean(row));
    SubImg_locs_artery(sectionIdx, 2) = round(mean(col));
    SubImg_width_artery(sectionIdx) = 0.01 * numX;
end

%% 1) 2) Veins

if veins_analysis

    maskSectionVein = maskSection .* maskVein;

    figure(120)
    imagesc(maskSectionVein .* v_RMS_AVG), axis image

    maskSectionVein = bwlabel(maskSectionVein);

    figure(121), axis image
    imagesc(maskSectionVein)

    numSectionsVein = max(maskSectionVein, [], 'all');
    masksSectionsVein = zeros(numX, numY, numSectionsVein);

    parfor sectionIdx = 1:numSectionsVein
        masksSectionsVein(:, :, sectionIdx) = (maskSectionVein == sectionIdx);
    end

    SubImg_locs_vein = zeros(numSectionsVein, 2);
    SubImg_width_vein = zeros(numSectionsVein, 1);

    for sectionIdx = 1:numSectionsVein
        [row, col] = find(masksSectionsVein(:, :, sectionIdx));
        SubImg_locs_vein(sectionIdx, 1) = round(mean(row));
        SubImg_locs_vein(sectionIdx, 2) = round(mean(col));
        SubImg_width_vein(sectionIdx) = 0.01 * numX;
    end

end

%% 2) Compute blood volume rate: Cross_section_analysis

strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'Velocity (mm.s-1)';

%% Arteries
[avgVolumeRateArtery, stdVolumeRateArtery, crossSectionAreaArtery, avgVelocityArtery, crossSectionMaskArtery, avgVolumeRateArtery_total, stdVolumeRateArtery_total] = crossSectionAnalysis(SubImg_locs_artery, SubImg_width_artery, maskArtery, v_RMS_all, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile, [],force_width);

labelsArteries = cell(numSectionsArtery, 1);
smoothVelocityArtery = zeros(numSectionsArtery, numFrames); % avg_blood_velocity (section, time)

parfor sectionIdx = 1:numSectionsArtery
    smoothVelocityArtery(sectionIdx, :) = smoothdata(avgVelocityArtery(sectionIdx, :), 'lowess');
    labelsArteries{sectionIdx} = sprintf("A%d", sectionIdx);
end

figure(210)
plot(fullTime, smoothVelocityArtery, 'LineWidth', 2)
title('Blood velocity in artery sections');
legend(labelsArteries);
fontsize(gca, 14, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel(strYlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArtery.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArtery.eps')))

smoothVelocityArtery_every_section = squeeze(mean(smoothVelocityArtery, 1));
smoothVelocityArtery_every_section_mean = mean(smoothVelocityArtery_every_section);

figure(221)
plot(fullTime, smoothVelocityArtery_every_section, 'k', 'LineWidth', 2)
yline(smoothVelocityArtery_every_section_mean, '--k', 'LineWidth', 2)
title('Blood velocity in artery sections');
legend(sprintf('All Artery Sections'));
fontsize(gca, 14, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel(strYlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArteryEverySection.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArteryEverySection.eps')))

for sectionIdx = 1:numSectionsArtery
    plot2txt(fullTime, avgVolumeRateArtery(sectionIdx, :), strcat('volumeRate_artery_A', num2str(sectionIdx)), ToolBox)
    plot2txt(fullTime, stdVolumeRateArtery(sectionIdx, :), strcat('volumeRate_artery_std_A', num2str(sectionIdx)), ToolBox)
    plot2txt(fullTime, avgVelocityArtery(sectionIdx, :), strcat('avg_velocity_artery_A', num2str(sectionIdx)), ToolBox)
end

%% 2) 2) Veins
if veins_analysis
    [avgVolumeRateVein, stdVolumeRateVein, crossSectionAreaVein, avgVelocityVein, crossSectionMaskVein, avgVolumeRateVein_total, stdVolumeRateVein_total] = crossSectionAnalysis(SubImg_locs_vein, SubImg_width_vein, maskVein, v_RMS_all, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'vein', flagBloodVelocityProfile, [],force_width);

    labelsVeins = cell(numSectionsVein, 1);
    smoothVelocityVein = zeros(numSectionsVein, numFrames);

    parfor sectionIdx = 1:numSectionsVein
        smoothVelocityVein(sectionIdx, :) = smoothdata(avgVelocityVein(sectionIdx, :), 'lowess');
        labelsVeins{sectionIdx} = sprintf("V%d", sectionIdx);
    end

    figure(220)
    plot(fullTime, smoothVelocityVein, 'LineWidth', 2)
    title('Blood velocity in vein sections');
    legend(labelsVeins);
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVein.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVein.eps')))

    smoothVelocityVein_every_section = squeeze(mean(smoothVelocityVein, 1));
    smoothVelocityVein_mean = mean(smoothVelocityVein_every_section);

    figure(221)
    plot(fullTime, smoothVelocityVein_every_section, 'k', 'LineWidth', 2)
    yline(smoothVelocityVein_mean, '--k', 'LineWidth', 2)
    title('Blood velocity in vein sections');
    legend(sprintf('All Vein Sections'));
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinEverySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinEverySection.eps')))

    for sectionIdx = 1:numSectionsVein
        plot2txt(fullTime, avgVolumeRateVein(sectionIdx, :), strcat('volumeRate_vein_V', num2str(sectionIdx)), ToolBox)
        plot2txt(fullTime, stdVolumeRateVein(sectionIdx, :), strcat('volumeRate_vein_std_V', num2str(sectionIdx)), ToolBox)
        plot2txt(fullTime, avgVelocityVein(sectionIdx, :), strcat('avg_velocity_vein_V', num2str(sectionIdx)), ToolBox)
    end

end

maskOnes = ones(numX, numY);
M0_disp_image = rescale(mean(M0_disp_video, 3));
ratio_etiquette = 1.2;

%% 3) Vein and artery numerotation

maskArtery_RGB(:, :, 1) = M0_disp_image .* ~crossSectionMaskArtery + crossSectionMaskArtery;
maskArtery_RGB(:, :, 2) = M0_disp_image .* ~crossSectionMaskArtery;
maskArtery_RGB(:, :, 3) = M0_disp_image .* ~crossSectionMaskArtery;

figure(300)
imshow(maskArtery_RGB);

if veins_analysis
    maskVein_RGB(:, :, 1) = M0_disp_image .* ~crossSectionMaskVein;
    maskVein_RGB(:, :, 2) = M0_disp_image .* ~crossSectionMaskVein;
    maskVein_RGB(:, :, 3) = M0_disp_image .* ~crossSectionMaskVein + crossSectionMaskVein;

    maskCombined_RGB(:, :, 1) = M0_disp_image .* ~(crossSectionMaskArtery & crossSectionMaskVein) + crossSectionMaskArtery;
    maskCombined_RGB(:, :, 2) = M0_disp_image .* ~(crossSectionMaskArtery & crossSectionMaskVein);
    maskCombined_RGB(:, :, 3) = M0_disp_image .* ~(crossSectionMaskArtery & crossSectionMaskVein) + crossSectionMaskVein;

    figure(301)
    imshow(maskVein_RGB);
    figure(302)
    imshow(maskCombined_RGB);
end

x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;

for sectionIdx = 1:numSectionsArtery
    new_x = x_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 2) - x_center);
    new_y = y_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 1) - y_center);
    figure(300)
    text(new_x, new_y, sprintf("A%d", sectionIdx), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");

    if veins_analysis
        figure(302)
        text(new_x, new_y, sprintf("A%d", sectionIdx), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

end

if veins_analysis

    for sectionIdx = 1:numSectionsVein
        new_x = x_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 1) - y_center);
        figure(301)
        text(new_x, new_y, sprintf("V%d", sectionIdx), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
        figure(302)
        text(new_x, new_y, sprintf("V%d", sectionIdx), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

end

drawnow
ax = gca;
ax.Units = 'pixels';

figure(300)
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'arteries_numerotation.png')))

if veins_analysis
    figure(301)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'veins_numerotation.png')))
    figure(302)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'arteries_veins_numerotation.png')))
end

%% 4) Volume Rate FigurePosition

Color_std = [0.7 0.7 0.7];
figWidth = 600;
figHeight = figWidth;

volumeRateArtery_video = zeros(figWidth, figHeight, 3, numFrames);
volumeRateArtery_mean = mean(avgVolumeRateArtery_total);

maskArtery_RGB(:, :, 1) = M0_disp_image .* ~crossSectionMaskArtery + crossSectionMaskArtery;
maskArtery_RGB(:, :, 2) = M0_disp_image .* ~crossSectionMaskArtery;
maskArtery_RGB(:, :, 3) = M0_disp_image .* ~crossSectionMaskArtery;

if veins_analysis

    volumeRateVein_video = zeros(figWidth, figHeight, 3, numFrames);
    volumeRateVein_mean = mean(avgVolumeRateVein_total);

    maskVein_RGB(:, :, 1) = M0_disp_image .* ~crossSectionMaskVein;
    maskVein_RGB(:, :, 2) = M0_disp_image .* ~crossSectionMaskVein;
    maskVein_RGB(:, :, 3) = M0_disp_image .* ~crossSectionMaskVein + crossSectionMaskVein;

end

%% 4) 1) Video Generation

M0_disp_video_rescaled = rescale(M0_disp_video);

if exportVideos

    for frameIdx = 1:numFrames

        figure(410, 'Position',[200 200 200+figWidth 200+figHeight])

        hueArtery = 0 * (maskOnes .* crossSectionMaskArtery);
        satArtery = maskOnes .* crossSectionMaskArtery;
        valArtery = maskOnes .* crossSectionMaskArtery;

        crossSectionArtery_RGB = hsv2rgb(hueArtery, satArtery, valArtery) .* crossSectionMaskArtery + ~crossSectionMaskArtery .* M0_disp_video_rescaled(:, :, frameIdx);
        imagesc(crossSectionArtery_RGB);
        axis image
        axis off

        for sectionIdx = 1:numSectionsArtery
            new_x = x_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 1) - y_center);

            if round(avgVolumeRateArtery(sectionIdx, frameIdx), 1) > 0
                text(new_x, new_y, string(round(avgVolumeRateArtery(sectionIdx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            else
                text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "red", "BackgroundColor", "black");
            end

        end

        title(sprintf("Blood volume rate : %02.0f µL/min in arteries", round(avgVolumeRateArtery_total(frameIdx))));
        set(gca, 'FontSize', 18)
        drawnow

        volumeRateArtery_frame = getframe(gcf);
        volumeRateArtery_video(:, :, :, frameIdx) = frame2im(volumeRateArtery_frame);

        if veins_analysis

            figure(411, 'Position',[200 200 200+figWidth 200+figHeight])

            hueVein = 0.6 * (maskOnes .* crossSectionMaskVein);
            satVein = maskOnes .* crossSectionMaskVein;
            valVein = maskOnes .* crossSectionMaskVein;

            crossSectionVein_RGB = hsv2rgb(hueVein, satVein, valVein) .* crossSectionMaskVein + ~crossSectionMaskVein .* M0_disp_video_rescaled(:, :, frameIdx);
            imagesc(crossSectionVein_RGB);
            axis image
            axis off

            for sectionIdx = 1:numSectionsVein
                new_x = x_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 2) - x_center);
                new_y = y_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 1) - y_center);

                if round(avgVolumeRateVein(sectionIdx, frameIdx), 1) > 0
                    text(new_x, new_y, string(round(avgVolumeRateVein(sectionIdx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
                else
                    text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "red", "BackgroundColor", "black");
                end

            end

            title(sprintf('Blood volume rate : %02.0f µL/min in veins', round(avgVolumeRateVein_total(frameIdx))));
            set(gca, 'FontSize', 18)
            drawnow

            volumeRateVein_frame = getframe(gcf);
            volumeRateVein_video(:, :, :, frameIdx) = frame2im(volumeRateVein_frame);

        end

    end

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(volumeRateArtery_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_volumeRateArtery_video.avi')));

    if veins_analysis
        parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(volumeRateVein_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_volumeRateVein_video.avi')));

    end

end

%% 4) 2) Volume Rate Figure

figure(420);
imshow(maskArtery_RGB);

for sectionIdx = 1:numSectionsArtery
    new_x = x_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 2) - x_center);
    new_y = y_center + ratio_etiquette * (SubImg_locs_artery(sectionIdx, 1) - y_center);
    text(new_x, new_y, string(round(mean(avgVolumeRateArtery(sectionIdx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
end

title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(avgVolumeRateArtery_total))));
drawnow
ax = gca;
ax.Units = 'pixels';
set(gca, 'FontSize', 14)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArteryImage.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArteryImage.eps')))

if veins_analysis

    figure(421);
    imshow(maskVein_RGB);

    for sectionIdx = 1:numSectionsVein
        new_x = x_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_vein(sectionIdx, 1) - y_center);
        text(new_x, new_y, string(round(mean(avgVolumeRateVein(sectionIdx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

    title(sprintf("Total blood volume rate : %02.0f µL/min in veins", round(mean(avgVolumeRateVein_total))));
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    set(gca, 'FontSize', 14)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinImage.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinImage.eps')))

end

%% 5) Plot Volume Rate
%% 5) 1) Artery plot

volumeRateArtery_plot = figure(510);
volumeRateArtery_plot.Position = [200 200 600 275];

curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

fill(fullTime2, inBetween, Color_std);
hold on;
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
yline(volumeRateArtery_mean, '--k', 'LineWidth', 2)
axis tight;
hold off

volumeRateArtery_ax = axis;
volumeRateArtery_ax(3) = 0;

ylabel('Blood volume rate (µL/min)')
xlabel('Time (s)')
title("Total blood volume rate in arteries")
axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)

volumeRateArtery_plot_frame = getframe(gcf);
volumeRateArtery_plot_video = zeros(size(volumeRateArtery_plot_frame.cdata, 1), size(volumeRateArtery_plot_frame.cdata, 2), 3, numFrames);

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArterySection.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArterySection.eps')))

plot2txt(fullTime, avgVolumeRateArtery_total, 'AVGVolumeRateArteryTotal', ToolBox)
plot2txt(fullTime, stdVolumeRateArtery_total, 'STDVolumeRateArteryTotal', ToolBox)

%% 5) 2) Vein plot

if veins_analysis

    volumeRateVein_plot = figure(520);
    volumeRateVein_plot.Position = [200 475 600 275];

    curve1 = avgVolumeRateVein_total + 0.5 * stdVolumeRateVein_total;
    curve2 = avgVolumeRateVein_total - 0.5 * stdVolumeRateVein_total;
    fullTime2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)];

    fill(fullTime2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, avgVolumeRateVein_total, '-k', 'LineWidth', 2);
    yline(volumeRateVein_mean, '--k', 'LineWidth', 2)
    axis tight;
    hold off

    volumeRateVein_ax = axis;
    volumeRateVein_ax(3) = 0;

    ylabel('Blood volume rate (µL/min)')
    xlabel('Time (s)')
    title("Total blood volume rate in veins")
    axis([volumeRateVein_ax(1) volumeRateVein_ax(2) volumeRateVein_ax(3) volumeRateVein_ax(4)]);
    fontsize(gca, 14, "points");
    set(gca, 'Linewidth', 2)

    volumeRateVein_plot_frame = getframe(gcf);
    volumeRateVein_plot_video = zeros(size(volumeRateVein_plot_frame.cdata, 1), size(volumeRateVein_plot_frame.cdata, 2), 3, numFrames);

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinSection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinSection.eps')))

    plot2txt(fullTime, avgVolumeRateVein_total, 'AVGVolumeRateVeinsTotal', ToolBox)
    plot2txt(fullTime, stdVolumeRateVein_total, 'STDVolumeRateVeinsTotal', ToolBox)

end

%% 5) 3) Artery Plot & Vein Plot progression
if exportVideos

    for frameIdx = 1:numFrames

        figure(530)

        curve1 = avgVolumeRateArtery_total(1:frameIdx) + 0.5 * stdVolumeRateArtery_total(1:frameIdx);
        curve2 = avgVolumeRateArtery_total(1:frameIdx) - 0.5 * stdVolumeRateArtery_total(1:frameIdx);
        tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
        inBetween = [curve1, fliplr(curve2)];

        fill(tmp_fullTime, inBetween, Color_std);
        hold on;
        plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), avgVolumeRateArtery_total(1:frameIdx), '-k', 'LineWidth', 2);
        yline(volumeRateArtery_mean, '--k', 'LineWidth', 2)
        hold off;

        ylabel('Blood volume rate (µL/min)')
        xlabel('Time (s)')
        title(sprintf("Total blood volume rate in arteries : %02.0f µL/min", round(volumeRateArtery_mean)))
        axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
        fontsize(gca, 14, "points");
        set(gca, 'Linewidth', 2)

        volumeRateArtery_plot_frame = getframe(gcf);
        volumeRateArtery_plot_video(:, :, :, frameIdx) = volumeRateArtery_plot_frame.cdata;

        set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
        set(gca, 'Linewidth', 2)

        if veins_analysis

            figure(540);

            curve1 = avgVolumeRateVein_total(1:frameIdx) + 0.5 * stdVolumeRateVein_total(1:frameIdx);
            curve2 = avgVolumeRateVein_total(1:frameIdx) - 0.5 * stdVolumeRateVein_total(1:frameIdx);
            tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
            inBetween = [curve1, fliplr(curve2)];

            fill(tmp_fullTime, inBetween, Color_std);
            hold on;
            plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
            plot(fullTime(1:frameIdx), avgVolumeRateVein_total(1:frameIdx), '-k', 'LineWidth', 2);
            yline(volumeRateVein_mean, '--k', 'LineWidth', 2)
            hold off;

            ylabel('Blood volume rate (µL/min)')
            xlabel('Time (s)')
            title(sprintf("Total blood volume rate in veins : %02.0f µL/min", round(volumeRateVein_mean)))
            axis([volumeRateVein_ax(1) volumeRateVein_ax(2) volumeRateVein_ax(3) volumeRateVein_ax(4)]);
            fontsize(gca, 14, "points");
            set(gca, 'Linewidth', 2)

            volumeRateVein_plot_frame = getframe(gcf);
            volumeRateVein_plot_video(:, :, :, frameIdx) = volumeRateVein_plot_frame.cdata;

            set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
            set(gca, 'Linewidth', 2)

        end

    end

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(volumeRateArtery_plot_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volumeRateArtery_video.avi')));

    timePeriod = ToolBox.stride / ToolBox.fs / 1000;

    combined_Gif_artery = cat(1, mat2gray(volumeRateArtery_video), mat2gray(volumeRateArtery_plot_video));
    writeGifOnDisc(combined_Gif_artery, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateArtery")), timePeriod);

    if veins_analysis

        w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_plot_volumeRateVein_video.avi')));

        tmp = mat2gray(volumeRateVein_plot_video);
        open(w)

        gifWriter = GifWriter(fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "volumeRateVein")), timePeriod, 0.04, numFrames);
        combined_Gif_vein = cat(1, mat2gray(volumeRateVein_video), mat2gray(volumeRateVein_plot_video));

        for frameIdx = 1:numFrames
            writeVideo(w, tmp(:, :, :, frameIdx));
            gifWriter.write(combined_Gif_vein(:, :, :, frameIdx), frameIdx);
        end

        gifWriter.generate();
        gifWriter.delete();

        imwrite(mat2gray(volumeRateVein_plot_video(:, :, :, end)), fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, "volumeRateVideoVein.png")));

        close(w);

    end

end

%% txt file output with measured pulse wave parameters

if veins_analysis
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, [ ...
                            'Value of total arterial blood volume rate (µL/min) :\n%d\n' ...
                        'Value of total venous blood volume rate (µL/min) :\n%d\n'], ...
        volumeRateArtery_mean, ...
        volumeRateVein_mean);
    fclose(fileID);
else
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, ...
        'Value of total arterial blood volume rate (µL/min) :\n%d\n', ...
        volumeRateArtery_mean);
    fclose(fileID);
end

for sectionIdx = 1:numSectionsArtery
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, [ ...
                            'Artery n°A%d : cross_section (mm^2) : \n %d \n' ...
                            'Artery n°A%d : vessel diameter (µm) : \n %d \n' ...
                            'Artery n°A%d : average velocity (mm/s) : \n %d \n' ...
                        'Artery n°A%d : blood volume rate (µL/min) : \n %d \n \n'], ...
        sectionIdx, ...
        crossSectionAreaArtery(sectionIdx), ...
        sectionIdx, ...
        2 * sqrt(crossSectionAreaArtery(sectionIdx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
        sectionIdx, ...
        avgVelocityArtery(sectionIdx), ...
        sectionIdx, ...
        avgVolumeRateArtery(sectionIdx)); % mm^3/s -> µL/min
    fclose(fileID);
end

if veins_analysis

    for sectionIdx = 1:numSectionsVein
        fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
        fprintf(fileID, [ ...
                                'Vein n°V%d : cross_section (mm^2) : \n %d \n ' ...
                                'Vein n°V%d : vessel diameter (µm) : \n %d \n ' ...
                                'Vein n°V%d : average velocity (mm/s) : \n %d \n ' ...
                            'Vein n°V%d : blood volume rate (µL/min) : \n %d \n \n'], ...
            sectionIdx, ...
            crossSectionAreaVein(sectionIdx), ...
            sectionIdx, ...
            2 * sqrt(crossSectionAreaVein(sectionIdx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
            sectionIdx, ...
            avgVelocityVein(sectionIdx), ...
            sectionIdx, ...
            avgVolumeRateVein(sectionIdx)); % mm^3/s -> µL/min
        fclose(fileID);
    end

end

close all

%% 6) Arterial Resistivity with Volume Rate
maxVolumeRate = max(avgVolumeRateArtery_total(:));
minVolumeRate = min(avgVolumeRateArtery_total(:));
meanVolumeRate = mean(avgVolumeRateArtery_total(:));

%% 6) 1) Arterial Resisitivity Index

ARI = (maxVolumeRate - minVolumeRate) / maxVolumeRate;

figure(610)

curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

hold on 

fill(fullTime2, inBetween, Color_std);
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
yline(volumeRateArtery_mean, '--k', 'LineWidth', 2)
yline(maxVolumeRate, '--r', 'Linewidth', 2)
yline(minVolumeRate, '--r', 'Linewidth', 2)
axis tight;


volumeRateArtery_ax = axis;
volumeRateArtery_ax(3) = 0;
hold off

ylabel('Blood volume rate (µL/min)')
xlabel('Time (s)')
title(sprintf("Arterial Resistivity %0.2f ", ARI))
axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);

fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateARI.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateARI.eps')))


%% 6) 2) Arterial Pulsatility Index

API = (maxVolumeRate - minVolumeRate) / meanVolumeRate;

figure(620)

curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

hold on 

fill(fullTime2, inBetween, Color_std);
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
yline(maxVolumeRate, '--r', 'Linewidth', 2)
yline(meanVolumeRate, '--r', 'Linewidth', 2)
yline(minVolumeRate, '--r', 'Linewidth', 2)
axis tight;

volumeRateArtery_ax = axis;
volumeRateArtery_ax(3) = 0;
hold off

ylabel('Blood volume rate (µL/min)')
xlabel('Time (s)')
title(sprintf("Arterial Pulsatility %d ", API))
axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateAPI.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateAPI.eps')))

close all


%% All circles testing

if ~PW_params.AllCirclesFlag
    return
end

% for the all circles output
numCircles = PW_params.numCircles;
maskSectionCircles = cell(1, numCircles);
deltr = (PW_params.velocity_bigRadiusRatio - PW_params.velocity_smallRadiusRatio) * (numY + numX) / 2 / numCircles; %PW_params.radius_gap

for i = 1:numCircles
    rad1 = (PW_params.velocity_smallRadiusRatio) * (numY + numX) / 2 + (i - 1) * deltr; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * deltr ;
    rad2 = rad1 + deltr;
    c1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad1;
    c2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad2;
    maskSectionCircles(i) = {xor(c1, c2)};

    % save mask image
    createMaskSection(M0_disp_image, maskArtery, rad1, rad2, sprintf('_mask_artery_section_circle_%d.png', i), ToolBox, path);
end

close(156);

% for all circles output

SubImg_locs_artery_Circles = zeros(numCircles, numSectionsArtery, 2);
SubImg_width_artery_Circles = zeros(numCircles, numSectionsArtery, 1);

for i = 1:numCircles
    maskSectionArtery = maskSectionCircles{i} .* maskArtery;

    maskSectionArtery = bwlabel(maskSectionArtery);

    numSectionsArtery = max(maskSectionArtery, [], 'all');
    masksSectionsArtery = zeros(numX, numY, numSectionsArtery);

    parfor sectionIdx = 1:numSectionsArtery

        masksSectionsArtery(:, :, sectionIdx) = (maskSectionArtery == sectionIdx);

    end

    for sectionIdx = 1:numSectionsArtery
        [row, col] = find(masksSectionsArtery(:, :, sectionIdx));
        SubImg_locs_artery_Circles(i, sectionIdx, 1) = round(mean(row));
        SubImg_locs_artery_Circles(i, sectionIdx, 2) = round(mean(col));
        SubImg_width_artery_Circles(i, sectionIdx) = 0.01 * size(maskArtery, 1);
    end

end

% for all circles output

for i = 1:numCircles
    [avgVolumeRateArtery, ~, ~, avgVelocityArtery, crossSectionMaskArtery, avgVolumeRateArtery_total, stdVolumeRateArtery_total] = cross_section_analysis(reshape(nonzeros(SubImg_locs_artery_Circles(i, :, :)), [], 2), nonzeros(SubImg_width_artery_Circles(i, :, :)), maskArtery, v_RMS_all, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile, i);

    if length(avgVolumeRateArtery) < 1
        continue
    end

    maskSectionArtery = maskSectionCircles{i} .* maskArtery;

    maskSectionArtery = bwlabel(maskSectionArtery);

    numSectionsArtery = max(maskSectionArtery, [], 'all');

    labelsArteries = cell(numSectionsArtery, 1);
    strXlabel = 'Time(s)'; %createXlabelTime(1);
    strYlabel = 'Velocity (mm.s-1)';

    data_to_plot_artery = zeros(size(avgVelocityArtery, 1), size(avgVolumeRateArtery, 2));

    parfor sectionIdx = 1:numSectionsArtery
        data_to_plot_artery(sectionIdx, :) = smoothdata(avgVelocityArtery(sectionIdx, :), 'lowess');
        labelsArteries{sectionIdx} = strcat('A', num2str(sectionIdx));
    end

    figure(57)
    plot(fullTime, data_to_plot_artery, 'LineWidth', 2)
    title('Blood velocity in artery sections');
    legend(labelsArteries);
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'velocityInArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'velocityInArterySection.eps')))

    data_to_plot_artery_all = squeeze(mean(data_to_plot_artery, 1));
    plot_artery_mean = mean(data_to_plot_artery_all);

    figure(571)
    plot(fullTime, data_to_plot_artery_all, 'k', 'LineWidth', 2)
    yline(plot_artery_mean, '--k', 'LineWidth', 2)
    title('Blood velocity in artery sections');
    legend(sprintf('All Artery Sections'));
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'velocityInAllArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'velocityInAllArterySection.eps')))

    Color_std = [0.7 0.7 0.7];
    figWidth = 600;
    figHeight = figWidth;

    maskArtery_RGB = ones(size(maskArtery, 1), size(maskArtery, 2), 3);
    volumeRateArtery_video = zeros(figWidth, figHeight, 3, numFrames);
    volumeRateArtery_mean = mean(avgVolumeRateArtery_total);

    maskArtery_RGB(:, :, 3) = M0_disp_image .* ~crossSectionMaskArtery;
    maskArtery_RGB(:, :, 2) = M0_disp_image .* ~crossSectionMaskArtery;
    maskArtery_RGB(:, :, 1) = M0_disp_image .* ~crossSectionMaskArtery + maskArtery_RGB(:, :, 1) .* crossSectionMaskArtery;

    f400 = figure(100);
    colormap("gray")
    f400.Position = [200, 200, 600, 600];
    ax400 = gca;
    ax400.Units = 'pixels';

    for frameIdx = 1:numFrames

        figure(100)

        hueArtery = 0 * (maskOnes .* crossSectionMaskArtery);
        satArtery = maskOnes .* crossSectionMaskArtery;
        valArtery = maskOnes .* crossSectionMaskArtery;

        crossSectionArtery_RGB = hsv2rgb(hueArtery, satArtery, valArtery) .* crossSectionMaskArtery + ~crossSectionMaskArtery .* M0_disp_video_rescaled(:, :, frameIdx);
        imagesc(ax400, crossSectionArtery_RGB);
        axis image
        axis off

        for sectionIdx = 1:numSectionsArtery
            new_x = x_center + ratio_etiquette * (SubImg_locs_artery_Circles(i, sectionIdx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (SubImg_locs_artery_Circles(i, sectionIdx, 1) - y_center);

            if round(avgVolumeRateArtery(sectionIdx, frameIdx), 1) > 0
                text(new_x, new_y, string(round(avgVolumeRateArtery(sectionIdx, frameIdx), 1)), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            else
                text(new_x, new_y, string(0), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
            end

        end

        tmp_bvra = round(avgVolumeRateArtery_total(frameIdx));
        tmp_bvra = (tmp_bvra > 0) * tmp_bvra;
        title(sprintf("Blood volume rate : %02.0f µL/min in arteries", tmp_bvra));
        set(gca, 'FontSize', 18)
        drawnow

        volumeRateArtery_frame = getframe(f400);
        volumeRateArtery_video(:, :, :, frameIdx) = frame2im(volumeRateArtery_frame);

    end

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(volumeRateArtery_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_ volumeRateArtery_video.avi')));

    figure(102);
    imshow(maskArtery_RGB);

    for sectionIdx = 1:numSectionsArtery
        new_x = x_center + ratio_etiquette * (SubImg_locs_artery_Circles(i, sectionIdx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (SubImg_locs_artery_Circles(i, sectionIdx, 1) - y_center);
        text(new_x, new_y, string(round(mean(avgVolumeRateArtery(sectionIdx, :), 2))), "FontWeight", "bold", "FontSize", 14, "Color", "white", "BackgroundColor", "black");
    end

    title(sprintf("Total blood volume rate : %02.0f µL/min in arteries", round(mean(avgVolumeRateArtery_total))));
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    set(gca, 'FontSize', 14)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'volumeRateInArteries.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'volumeRateInArteries.eps')))
    % F_Total_volumeRate_artery = getframe(f102);

    volumeRateArtery_plot = figure(104);
    volumeRateArtery_plot.Position = [200 200 600 275];

    %plot(fullTime, total_avg_volumeRate_artery, '-k', 'LineWidth', 2);
    curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
    curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
    fullTime2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)];

    fill(fullTime2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
    axis tight;
    hold off

    volumeRateArtery_ax = axis;
    volumeRateArtery_ax(3) = 0;

    ylabel('Blood volume rate (µL/min)')
    xlabel('Time (s)')
    title("Total blood volume rate in arteries")
    axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
    fontsize(gca, 14, "points");
    set(gca, 'Linewidth', 2)

    volumeRateArtery_plot_frame = getframe(gcf);
    volumeRateArtery_plot_video = zeros(size(volumeRateArtery_plot_frame.cdata, 1), size(volumeRateArtery_plot_frame.cdata, 2), 3, numFrames);

    for frameIdx = numFrames:numFrames

        figure(104)

        curve1 = avgVolumeRateArtery_total(1:frameIdx) + 0.5 * stdVolumeRateArtery_total(1:frameIdx);
        curve2 = avgVolumeRateArtery_total(1:frameIdx) - 0.5 * stdVolumeRateArtery_total(1:frameIdx);
        tmp_fullTime = [fullTime(1:frameIdx), fliplr(fullTime(1:frameIdx))];
        inBetween = [curve1, fliplr(curve2)];

        fill(tmp_fullTime, inBetween, Color_std);
        hold on;
        plot(fullTime(1:frameIdx), curve1, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), curve2, "Color", Color_std, 'LineWidth', 2);
        plot(fullTime(1:frameIdx), avgVolumeRateArtery_total(1:frameIdx), '-k', 'LineWidth', 2);
        yline(volumeRateArtery_mean, '--k', 'LineWidth', 2)
        hold off;

        ylabel('Blood volume rate (µL/min)')
        xlabel('Time (s)')
        title(sprintf("Total blood volume rate in arteries : %02.0f µL/min", round(volumeRateArtery_mean)))
        axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
        fontsize(gca, 14, "points");
        set(gca, 'Linewidth', 2)

        volumeRateArtery_plot_frame = getframe(gcf);
        volumeRateArtery_plot_video(:, :, :, frameIdx) = volumeRateArtery_plot_frame.cdata;

        set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
        set(gca, 'Linewidth', 2)

    end

    figure(104)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'volumeRateInArterySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, 'volumeRateInArterySection.eps')))

    %plot2txt(fullTime, total_avg_volumeRate_artery, 'TotalVolumeRateArteryAVG', ToolBox)
    %plot2txt(fullTime, total_std_volumeRate_artery, 'TotalVolumeRateArterySTD', ToolBox)

    %parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(volumeRateArtery_plot_video),fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, 'circle_',num2str(i),'_plot_volumeRateArtery_video.avi')));

    %timePeriod = ToolBox.stride / ToolBox.fs / 1000;

    %combined_Gif_artery = cat(1, mat2gray(volumeRateArtery_video), mat2gray(volumeRateArtery_plot_video));
    %writeGifOnDisc(combined_Gif_artery,fullfile(ToolBox.PW_path_gif, sprintf("%s_circle_%d_%s.gif", ToolBox.PW_folder_name,i, "volumeRateArtery")),timePeriod);

    imwrite(mat2gray(volumeRateArtery_plot_video(:, :, :, end)), fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, i, "volumeRateVideoArtery.png")));

end

end
