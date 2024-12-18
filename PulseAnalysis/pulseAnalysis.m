function [v_RMS_video] = pulseAnalysis(f_RMS_video, maskArtery, maskVein, maskSection)
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
%   flagExtended    Size: 1
%
% Output:
%   v_RMS_video     Size: numX x numY x numFrames double

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
veinsAnalysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;

maskArterySection = maskArtery & maskSection;
maskVeinSection = maskVein & maskSection;

mkdir(ToolBox.PW_path_png, 'pulseAnalysis')
mkdir(ToolBox.PW_path_eps, 'pulseAnalysis')
folder = 'pulseAnalysis';

[numX, numY, numFrames] = size(f_RMS_video);
strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

%% 1) Local BKG Artery and Veins %~1min

tic

if veinsAnalysis
    maskVesselDilated = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width));
    imwrite(maskVesselDilated, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');
    imwrite(maskVesselDilated, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');

else
    maskVesselDilated = imdilate(maskArtery, strel('disk', PW_params.local_background_width));
    imwrite(maskVesselDilated, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');
    imwrite(maskVesselDilated, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');
end

f_RMS_background = zeros(numX, numY, numFrames, 'single');

w =  PW_params.local_background_width;
k =  PW_params.k;

parfor frameIdx = 1:numFrames
    f_RMS_background(:, :, frameIdx) = single(maskedAverage(f_RMS_video(:, :, frameIdx), 10 * w * 2^k, ~maskVesselDilated));
end

graphSignal('1_Arteries_fRMS', folder, ...
    t, squeeze(sum(f_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '-', cArtery, ...
    t, squeeze(sum(f_RMS_background .* maskArterySection, [1, 2]) / nnz(maskArterySection)), '--', cBlack, ...
    Title = 'Average f_{RMS} in Arteries', xlabel = strXlabel, ylabel = strYlabel, ...
    Legend = {'Arteries', 'Background'});

if veinsAnalysis
    graphSignal('1_Veins_fRMS', folder, ...
        t, squeeze(sum(f_RMS_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)), '-', cVein, ...
        t, squeeze(sum(f_RMS_background .* maskVeinSection, [1, 2]) / nnz(maskVeinSection)), '--', cBlack, ...
        Title = 'Average f_{RMS} in Veins', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Veins', 'Background'});
end

fprintf("    1. Local BKG Artery and Veins calculation took %ds\n", round(toc))

%% 2) Difference calculation

tic

if PW_params.DiffFirstCalculationsFlag == 0 %SIGNED DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    delta_f_RMS = sign(tmp) .* sqrt(abs(tmp));
    clear tmp

elseif PW_params.DiffFirstCalculationsFlag == 1 % DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    tmp = tmp .* (tmp > 0);
    delta_f_RMS = sqrt(tmp);
    clear tmp

else % DIFFERENCE LAST

    delta_f_RMS = f_RMS_video - f_RMS_background;

end

v_RMS_video = ToolBox.ScalingFactorVelocityInPlane * delta_f_RMS;

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

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];
LocalBackground_in_vessels = mean(f_RMS_background, 3) .* maskVesselDilated + ones(numX, numY) * mean(f_RMS_background, 'all') .* ~maskVesselDilated;
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
imwrite(rescale(LocalBackground_in_vessels), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '6_LocalBackground_in_vessels.png')))

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];
in_vessels = mean(delta_f_RMS, 3) .* maskVesselDilated;
imagesc(in_vessels);
colormap gray
title('Delta f in vessels');
fontsize(gca, 14, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
imwrite(rescale(in_vessels), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '6_Df_in_vessels.png')))

if exportVideos
    f_RMS_video_rescale = rescale(f_RMS_video);
    f_RMS_background_rescale = rescale(f_RMS_background);

    writeGifOnDisc(f_RMS_background_rescale, "f_RMS_bkg")
    writeGifOnDisc(f_RMS_video_rescale, "f_RMS")

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.avi')));
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.mp4')), 'MPEG-4');

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_video), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'f_AVG_vessels.avi')));
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_video), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'f_AVG_vessels.mp4')), 'MPEG-4');
end

fprintf("    2. Difference calculation took %ds\n", round(toc))

clear LocalBackground_in_vessels f_RMS_background

return;

end
