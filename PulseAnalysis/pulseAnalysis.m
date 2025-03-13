function [v_RMS_video] = pulseAnalysis(f_RMS_video, maskArtery, maskVein, maskSection, maskNeighbors)
% pulseAnalysis.m computes the velocities
% Inputs:
%       VIDEOS:
%   f_RMS_video     Size: numX x numY x numFrames double
%   M0_disp_video   Size: numX x numY x numFrames double
%       IMAGES:
%   f_AVG_image     Size: numX x numY double
%   maskArtery      Size: numX x numY logical
%   maskBackground  Size: numX x numY logical
%   maskSection     Size: numX x numY logical
%   maskVein        Size: numX x numY logical
%       TRIVIA:
%   sysIdxList:     Size: numSystoles
%
% Output:
%   v_RMS_video     Size: numX x numY x numFrames double

TB = getGlobalToolBox;
params = TB.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;

maskArterySection = maskArtery & maskSection;
maskVeinSection = maskVein & maskSection;

mkdir(TB.path_png, 'pulseAnalysis')
mkdir(TB.path_eps, 'pulseAnalysis')
folder = 'pulseAnalysis';

[numX, numY, numFrames] = size(f_RMS_video);
strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * TB.stride / TB.fs / 1000, numFrames);
cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

%% 1) Local BKG Artery and Veins %~1min

tic

f_RMS_background = zeros(numX, numY, numFrames, 'single');

w =  params.local_background_width;
k =  params.k;

if veinsAnalysis
    maskVessel = maskArtery | maskVein;
else
    maskVessel = maskArtery;
end

parfor frameIdx = 1:numFrames
    f_RMS_background(:, :, frameIdx) = single(maskedAverage(f_RMS_video(:, :, frameIdx), 10 * w * 2^k, maskNeighbors, maskVessel));
end

imwrite(rescale(squeeze(mean(f_RMS_background, 3))), fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_frequency_RMS_bkg.png')), 'png');

graphSignal('1_Arteries_fRMS', folder, ...
    t, squeeze(sum(f_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '-', cArtery, ...
    t, squeeze(sum(f_RMS_background .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '--', cBlack, ...
    Title = 'Average f_{RMS} in Arteries', xlabel = strXlabel, ylabel = strYlabel, ...
    Legend = {'Arteries', 'Local Background'});

fileID = fopen(fullfile(TB.path_txt, strcat(TB.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
fprintf(fileID, 'Mean fRMS difference artery : %f (kHz) \r\n',mean(squeeze(sum(f_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection)))-mean(squeeze(sum(f_RMS_background .* maskArterySection, [1, 2]) / nnz(maskArterySection))));
fclose(fileID);

if veinsAnalysis
    graphSignal('1_Veins_fRMS', folder, ...
        t, squeeze(sum(f_RMS_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)), '-', cVein, ...
        t, squeeze(sum(f_RMS_background .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)), '--', cBlack, ...
        Title = 'Average f_{RMS} in Veins', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Veins', 'Local Background'});
    fileID = fopen(fullfile(TB.path_txt, strcat(TB.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
    fprintf(fileID, 'Mean fRMS difference vein : %f (kHz) \r\n',mean(squeeze(sum(f_RMS_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)))-mean(squeeze(sum(maskVeinSection .* maskVeinSection, [1, 2]) / nnz(maskVeinSection))));
    fclose(fileID);
end

fprintf("    1. Local BKG Artery and Veins calculation took %ds\n", round(toc))

%% 2) Difference calculation

tic

if params.DiffFirstCalculationsFlag == 0 %SIGNED DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    delta_f_RMS = sign(tmp) .* sqrt(abs(tmp));
    clear tmp

elseif params.DiffFirstCalculationsFlag == 1 % DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    tmp = tmp .* (tmp > 0);
    delta_f_RMS = sqrt(tmp);
    clear tmp

else % DIFFERENCE LAST

    delta_f_RMS = f_RMS_video - f_RMS_background;

end

v_RMS_video = TB.ScalingFactorVelocityInPlane * delta_f_RMS;

if veinsAnalysis
    graphSignal('2_Vessels_velocity', folder, ...
        t, squeeze(sum(v_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '-', cArtery, ...
        t, squeeze(sum(v_RMS_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)), '-', cVein, ...
        Title = 'Average estimated velocity in Arteries and Veins', xlabel = strXlabel, ylabel = 'mm/s');

else
    graphSignal('2_Arteries_velocity', folder, ...
        t, squeeze(sum(v_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '-', cArtery, ...
        Title = 'Average estimated velocity in Arteries', xlabel = strXlabel, ylabel = 'mm/s');

end
fprintf("    2. Difference calculation took %ds\n", round(toc))

ArterialResistivityIndex(t, v_RMS_video, maskArtery, 'velocity', folder);

%% 3) Plots of f_RMS mean Local Background in vessels and Delta frequency in vessels and their colorbars
tic

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];

LocalBackground_in_vessels = mean(f_RMS_background, 3) .* maskVessel + ones(numX, numY) * mean(sum(f_RMS_background.*maskVessel, [1,2])/nnz(maskVessel),3) .* ~maskVessel;
imagesc(LocalBackground_in_vessels);
colormap gray
title('Local Background in vessels');
fontsize(gca, 14, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;

imwrite(rescale(LocalBackground_in_vessels), fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_LocalBackground_in_vessels.png')))

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
LocalBackground_colorbar = colorbar('north');
clim(range)
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
LocalBackground_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(LocalBackground_colorbar, 'Title');
titleString = 'Local Background RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorBarLocalBackground_in_vessels.png')))
exportgraphics(gca, fullfile(TB.path_eps, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorBarLocalBackground_in_vessels.eps')))



f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];
in_vessels = mean(delta_f_RMS, 3) .* maskVessel;
imagesc(in_vessels);
colormap gray
title('Delta f in vessels');
fontsize(gca, 14, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'Delta Doppler RMS frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(in_vessels), fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_Df_in_vessels.png')))

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
Df_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
Df_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(Df_colorbar, 'Title');
titleString = 'Delta Doppler RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorBarDf_in_vessels.png')))
exportgraphics(gca, fullfile(TB.path_eps, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorBarDf_in_vessels.eps')))


figure("Visible", "off")
imagesc(squeeze(mean(f_RMS_video, 3)));
colormap gray
title('RMS frequency map RAW');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(squeeze(mean(f_RMS_video, 3))), fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_frequency_RMS.png')), 'png');

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
f_RMS_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
f_RMS_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_RMS_colorbar, 'Title');
titleString = 'RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(TB.path_png, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorbarRMSFrequency.png')))
exportgraphics(gca, fullfile(TB.path_eps, 'pulseAnalysis', sprintf("%s_%s", TB.main_foldername, '3_ColorbarRMSFrequency.eps')))

fprintf("    3. Plotting heatmaps took %ds\n", round(toc))


if exportVideos
    f_RMS_video_rescale = rescale(f_RMS_video);
    f_RMS_background_rescale = rescale(f_RMS_background);

    writeGifOnDisc(f_RMS_background_rescale, "f_RMS_bkg")
    writeGifOnDisc(f_RMS_video_rescale, "f_RMS")
end

clear LocalBackground_in_vessels f_RMS_background

return;

end
