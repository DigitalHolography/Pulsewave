function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(M0_ff_video, ~, f_AVG_mean, path, ToolBox)

PW_params = Parameters_json(path);
exportVideos = PW_params.exportVideos;

mkdir(ToolBox.PW_path_png, 'mask')
mkdir(fullfile(ToolBox.PW_path_png, 'mask'), 'steps')
mkdir(fullfile(ToolBox.PW_path_eps, 'mask'), 'steps')
close all

%% 1) First Masks and Correlation

[numX, numY, numFrames] = size(M0_ff_video);
[X, Y] = meshgrid(1:numX, 1:numY);
maskDiaphragm = sqrt((X - numX / 2) .^ 2 + (Y - numY / 2) .^ 2) <= PW_params.masks_diaphragmRadius * (numY + numX) / 2;
imwrite(rescale(maskDiaphragm), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_0_maskDiaphragm.png')))

M0_ff_img = squeeze(mean(M0_ff_video, 3));
M0_ff_video_centered = M0_ff_video - M0_ff_img;
imwrite(rescale(M0_ff_img), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_0_M0.png')))

if exportVideos
    timePeriod = ToolBox.stride / ToolBox.fs / 1000;
    M0_ff_rescaled_video = rescale(M0_ff_video);
    writeGifOnDisc(M0_ff_rescaled_video, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "M0")), timePeriod)
end

%% 1) 1) Compute vesselness response
vesselnessM0 = vesselness_filter(M0_ff_img, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);
maskVesselness = logical(imbinarize(vesselnessM0 .* maskDiaphragm));

imwrite(rescale(vesselnessM0), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_1_Vesselness.png')))
imwrite(maskVesselness, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_2_vesselMask.png')))

%% 1) 2) Compute the barycentres and the circle mask

if ~isempty(PW_params.forcebarycenter)

    ToolBox.y_barycentre = PW_params.forcebarycenter(1);
    ToolBox.x_barycentre = PW_params.forcebarycenter(2);
    cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

else

    vascularImage = double(M0_ff_img .* f_AVG_mean);
    blurred_mask = imgaussfilt(vascularImage, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0);
    [ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));
    [y_CRV, x_CRV] = find(blurred_mask == min(blurred_mask, [], 'all'));
    cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;
    cercleMask = cercleMask | sqrt((X - x_CRV) .^ 2 + (Y - y_CRV) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

end

maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | cercleMask, 1, 4);
imwrite(rescale(maskVesselnessClean | cercleMask), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_3_choroidClean.png')))

%%  1) 3) Compute first correlation
% compute pulse in 3 dimentions for correlation in all vessels
vascularPulse = mean(M0_ff_video .* maskVesselnessClean, [1 2]);
vascularPulse = vascularPulse ./ nnz(maskVesselnessClean);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
figure, plot(t, squeeze(vascularPulse), 'k-', 'LineWidth', 2)
title('Vascular Pulse')
fontsize(gca, 14, "points");
xlabel('Time(s)', 'FontSize', 14);
ylabel('Power Doppler Pulse (s.u.)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vascularPulse.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vascularPulse.eps')))

vascularPulse_centered = vascularPulse - mean(vascularPulse, 3);

% compute local-to-average pulse wave zero-lag correlation
R_VascularPulse = mean(M0_ff_video_centered .* vascularPulse_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularPulse_centered, [], 3));
imwrite(rescale(R_VascularPulse), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_4_Correlation.png')))

vascularImage = double(M0_ff_img .* f_AVG_mean .* R_VascularPulse);
blurred_mask = imgaussfilt(vascularImage, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0) .* maskDiaphragm;
[ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));
cercleMask = cercleMask | sqrt((X - x_CRV) .^ 2 + (Y - y_CRV) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

%% 1) 4) Segment Vessels
cmapArtery = [0 0 0; 255 22 18]/255;
cmapVein = [0 0 0; 18 237 255]/255;

% Number of classes for Vessels: 5
% 1 = Background, 2 & 3 = Veins & CoroidalVessels, 4 = CoroidalVessel, 5 = Arteries
numClassesVessels = 5;
firstThresholds = multithresh(R_VascularPulse(maskVesselnessClean), numClassesVessels - 2);
firstThresholds = [-1 firstThresholds];
quantizedVesselCorrelation = imquantize(R_VascularPulse - ~maskVesselnessClean * 2, firstThresholds);

imwrite(rescale(quantizedVesselCorrelation), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_4_quantizedCorrelation.png')))

% Create first artery mask & first vein mask
firstMaskArtery = quantizedVesselCorrelation == 4 | quantizedVesselCorrelation == 5;
firstMaskVein = quantizedVesselCorrelation == 3 | quantizedVesselCorrelation == 2;

imwrite(firstMaskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_1_5_FirstMask.png')))
imwrite(firstMaskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_1_5_FirstMask.png')))

% Remove small blobs
% minSize = PW_params.masks_minSize * (numX * numY);
firstMaskArteryClean = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVeinClean = bwareaopen(firstMaskVein, PW_params.masks_minSize);

imwrite(firstMaskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_1_6_FirstMaskClean.png')))
imwrite(firstMaskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_1_6_FirstMaskClean.png')))


%% 1) 5) WIP: Choroid Segmentations
cmapChoroid = [0 0 0; 23 255 18] / 255;
numClassesChoroidalVessels = 5;
firstThresholdsChoroid = multithresh(R_VascularPulse(maskVesselness), numClassesChoroidalVessels - 2);
firstThresholdsChoroid = [-1 firstThresholdsChoroid];
quantizedChoroidalVesselCorrelation = imquantize(R_VascularPulse - ~maskVesselness * 2, firstThresholdsChoroid);

firstMaskChoroid = quantizedChoroidalVesselCorrelation == 4 | quantizedChoroidalVesselCorrelation == 3 | quantizedChoroidalVesselCorrelation == 2;
firstMaskChoroid = firstMaskChoroid & ~firstMaskArtery & ~firstMaskVein;
imwrite(firstMaskChoroid, cmapChoroid, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroid_1_6_FirstMask.png')))

RGBM0(:, :, 1) = rescale(M0_ff_img) + firstMaskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img) + firstMaskChoroid;
RGBM0(:, :, 3) = rescale(M0_ff_img) + firstMaskVeinClean;

imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_6_RGB.png')))

clear quantizedVesselCorrelation firstThresholds;
%% 2)  Compute the second correlation to fine tune the segmentation

%% 2) 1) Compute the new pulse

% compute pulse in 3 dimentions for correlation in main arteries
arterialPulse = mean(M0_ff_video .* firstMaskArteryClean .* maskDiaphragm, [1 2]);
arterialPulse = arterialPulse ./ nnz(firstMaskArteryClean .* maskDiaphragm);
arterialPulse_centered = arterialPulse - mean(arterialPulse, 3);

figure, plot(t, squeeze(arterialPulse), '-', 'LineWidth', 2, 'Color', [255/255 22/255 18/255])
title('Arterial Pulse')
fontsize(gca, 14, "points");
xlabel('Time(s)', 'FontSize', 14);
ylabel('Power Doppler Pulse (s.u.)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'arterialPulse.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'arterialPulse.eps')))

% compute pulse in 3 dimentions for correlation in main veins
venousPulse = mean(M0_ff_video .* firstMaskVeinClean, [1 2]);
venousPulse = venousPulse ./ nnz(firstMaskVeinClean);
venousPulse_centered = venousPulse - mean(venousPulse, 3);

figure, plot(t, squeeze(venousPulse), '-', 'LineWidth', 2, 'Color', [18/255 237/255 255/255])
title('Venous Pulse')
fontsize(gca, 14, "points");
xlabel('Time(s)', 'FontSize', 14);
ylabel('Power Doppler Pulse (s.u.)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'venousPulse.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'venousPulse.eps')))

% compute pulse in 3 dimentions for correlation in choroids

choroidalPulse = mean(M0_ff_video .* firstMaskChoroid, [1 2]);
choroidalPulse = choroidalPulse ./ nnz(firstMaskChoroid);
choroidalPulse_centered = choroidalPulse - mean(choroidalPulse, 3);

figure, plot(t, squeeze(choroidalPulse), '-', 'LineWidth', 2, 'Color', [18/255 255/255 23/255])
title('Choroidal Pulse')
fontsize(gca, 14, "points");
xlabel('Time(s)', 'FontSize', 14);
ylabel('Power Doppler Pulse (s.u.)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroidalPulse.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroidalPulse.eps')))

% compute local-to-average pulse wave zero-lag correlation
R_ArterialPulse = mean((M0_ff_video_centered .* arterialPulse_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(arterialPulse_centered, [], 3));
R_VenousPulse = mean((M0_ff_video_centered .* venousPulse_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(arterialPulse_centered, [], 3));
R_ChoroidalPulse = mean((M0_ff_video_centered .* choroidalPulse_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(choroidalPulse_centered, [], 3));

imwrite(rescale(R_ArterialPulse), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_1_CorrelMatrix.png')))
imwrite(rescale(R_VenousPulse), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_1_CorrelMatrix.png')))
imwrite(rescale(R_ChoroidalPulse), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroid_2_1_CorrelMatrix.png')))

figure
hold on
plot(t, squeeze(vascularPulse), 'k--', 'LineWidth', 2)
plot(t, squeeze(arterialPulse), '-', 'LineWidth', 2, 'Color', [255/255 22/255 18/255])
plot(t, squeeze(venousPulse), '-', 'LineWidth', 2, 'Color', [18/255 237/255 255/255])
title('Pulses used for correlation maps')
legend('Vascular Pulse', 'Arterial Pulse', 'Venous Pulse')
fontsize(gca, 14, "points");
xlabel('Time(s)', 'FontSize', 14);
ylabel('Power Doppler Pulse (s.u.)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
box on
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'everyPulse.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'everyPulse.eps')))

clear arterialPulse arterialPulse_centered venousPulse venousPulse_centered vascularPulse

%% 2) 2) New Vesselness

se = strel('disk', 5);
% maskVesselnessDilated = imdilate(maskVesselness, se); % As Arteries tend to be more dilated than what they appear to be

R_ArteryVessel = R_ArterialPulse .* maskVesselnessClean;
R_ArteryVessel = R_ArteryVessel .* (R_ArteryVessel > 0);

R_VeinVessel = R_VenousPulse .* maskVesselnessClean;
R_VeinVessel = R_VeinVessel .* (R_VeinVessel > 0);

imwrite(rescale(R_ArteryVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_2_R_ArteryVessel.png')))
imwrite(rescale(R_VeinVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_2_R_VeinVessel.png')))

% Segmentation based on Otsu's criteria for arteries and veins
numClassesArtery = 4;

levelArtery = multithresh(R_ArteryVessel(maskVesselnessClean), numClassesArtery - 2);
levelArtery = [-1 levelArtery];
maskArteryQ = imquantize(R_ArteryVessel - 2 * ~maskVesselnessClean, levelArtery);
maskArtery = maskArteryQ == 4 | firstMaskArteryClean;

imwrite(rescale(maskArteryQ), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_3_Quantize.png')))
imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_4_Thresh.png')))

numClassesVein = 4;

levelVein = multithresh(R_VeinVessel(maskVesselnessClean), numClassesVein - 2);
levelVein = [-1 levelVein];
maskVeinQ = imquantize(R_VeinVessel - 2 * ~maskVesselnessClean, levelVein);
maskVein = maskVeinQ == 4 | firstMaskVeinClean;

imwrite(rescale(maskVeinQ), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_3_Quantize.png')))
imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_4_Thresh.png')))

%% Clearing 3)
%% 3) 1) Morphological Operations

% Remove the small bits
maskArteryClean = maskArtery & bwareaopen(maskArtery | maskVein, PW_params.masks_minSize, 4);
maskVeinClean = maskVein & bwareaopen(maskArtery | maskVein, PW_params.masks_minSize, 4);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_1_AreaOpen.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_1_AreaOpen.png')))

maskArteryClean = imclose(maskArteryClean, se);
maskVeinClean = imclose(maskVeinClean, se);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_2_ImClose.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_2_ImClose.png')))

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVeinClean;

imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_2_RGB.png')))

clear vesselnessVein vesselnessArtery ;
%% 3) 2) Cleaning Choroid

maskArtery = maskArteryClean;
maskVein = maskVeinClean;
maskVessel = maskArtery | maskVein;

imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_3_Final.png')))
imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_3_Final.png')))
imwrite(maskVessel, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_3_Final.png')))
imwrite(rescale(uint8(cat(3, uint8(M0_ff_img) + uint8(maskArteryClean) * 255, uint8(M0_ff_img) + uint8(cercleMask), uint8(M0_ff_img) + uint8(maskVeinClean) * 255))), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_4_Final.png')))

%% SENSITIVITY & SPECIFICITY
% Add targetMaskArtery.png and targetMaskVein.png in a mask folder inside
% the pulsewave folder to have the specificity and the senisivity of the
% algorithm used with your parameters.
try
    targetMaskArtery = imread(fullfile(ToolBox.PW_path_main, 'mask', 'targetMaskArtery.png'));
    targetMaskVein = imread(fullfile(ToolBox.PW_path_main, 'mask', 'targetMaskVein.png'));

    TPArtery = nnz(targetMaskArtery & maskArtery);
    FPArtery = nnz(~targetMaskArtery & maskArtery);
    FNArtery = nnz(targetMaskArtery & ~maskArtery);
    TNArtery = nnz(~targetMaskArtery & ~maskArtery);

    TPVein = nnz(targetMaskVein & maskVein);
    FPVein = nnz(~targetMaskVein & maskVein);
    FNVein = nnz(targetMaskVein & ~maskVein);
    TNVein = nnz(~targetMaskVein & ~maskVein);

    fprintf("Artery Sensitivity: %0.2f %%\n", 100 * TPArtery / (TPArtery + FNArtery))
    fprintf("Artery Specificity: %0.2f %%\n", 100 * TNArtery / (TNArtery + FPArtery))
    fprintf("Vein Sensitivity: %0.2f %%\n", 100 * TPVein / (TPVein + FNVein))
    fprintf("Vein Specificity: %0.2f %%\n", 100 * TNVein / (TNVein + FPVein))
catch
    fprintf("No target masks to perform sensitivity & specificity evaluation\n")
end

%% Force Create Masks
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')) && isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png'))
    maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
    maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;
end

%% Create Background Mask
maskBackground = not(maskVessel);

%% Create CRA and CRV Mask
f_AVG_std = std2(f_AVG_mean);
maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);
clear f_AVG_std f_AVG_mean

%% Create Mask Section
radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (numY + numX) / 2;
radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (numY + numX) / 2;

circleMask1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius1;
circleMask2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius2;

maskSection = xor(circleMask1, circleMask2);

%% Create Colormap Artery/Vein

M0_ff_img = mat2gray(M0_ff_img);
[hueArtery, satArtery, val] = createHSVmap(M0_ff_img, maskArtery - maskArtery .* maskSection, 0, 0);
[hueVein, satVein, ~] = createHSVmap(M0_ff_img, maskVein - maskVein .* maskSection - maskVein .* maskArtery, 0.7, 0.7);
[hueSectionA, satSectionA, ~] = createHSVmap(M0_ff_img, maskSection .* maskArtery, 0.15, 0.15);
[hueSectionV, satSectionV, ~] = createHSVmap(M0_ff_img, maskSection .* maskVein, 0.5, 0.5);
satSectionArtery = satSectionA;
satSectionVein = satSectionV; %+~maskSectionArtery.*(~maskArtery);
val = val .* (~maskSection) + val .* maskSection + maskSection .* (~(maskArtery + maskVein));
vesselImageRGB = hsv2rgb(hueArtery + hueVein + hueSectionA + hueSectionV, satArtery + satVein + satSectionArtery + satSectionVein, val);

figure(101)
imshow(vesselImageRGB)

%% Create Colormap Artery Only
[hueArtery, satArtery, val] = createHSVmap(M0_ff_img, maskArtery - maskArtery .* maskSection, 0, 0);
[hueSectionA, satSectionA, ~] = createHSVmap(M0_ff_img, maskSection .* maskArtery, 0.15, 0.15);
satSectionArtery = satSectionA;
val = val .* (~maskSection) + val .* maskSection + maskSection .* (~maskArtery);
VesselImageRGB_Artery = hsv2rgb(hueArtery + hueSectionA, satArtery + satSectionArtery, val);

figure(102)
imshow(VesselImageRGB_Artery)

%% Create Segmentation Map

segmentationMap = zeros(numX, numY, 3);
segmentationMap(:, :, 1) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img + maskArtery;
segmentationMap(:, :, 2) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img;
segmentationMap(:, :, 3) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img + maskVein;
segmentationMapArtery(:, :, 1) = M0_ff_img - (maskArtery) .* M0_ff_img + maskArtery;
segmentationMapArtery(:, :, 2) = M0_ff_img - maskArtery .* M0_ff_img;
segmentationMapArtery(:, :, 3) = M0_ff_img - maskArtery .* M0_ff_img;
imwrite(segmentationMap, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arteryVeinSegmentation.png')), 'png');
imwrite(segmentationMapArtery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arterySegmentation.png')), 'png');

%% Saving masks as PNG
foldername = ToolBox.main_foldername;

imwrite(mat2gray(double(vesselnessM0)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselness.png')), 'png');
imwrite(mat2gray(single(maskArtery)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskArtery.png')), 'png');
imwrite(mat2gray(single(maskVein)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVein.png')), 'png');
imwrite(mat2gray(single(maskVessel)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVessel.png')), 'png');
imwrite(mat2gray(single(maskBackground)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskBackground.png')), 'png');
%vesselMap = uint8( cat( 3, uint8(M0_disp_img)+ uint8(maskArtery)*255, uint8(M0_disp_img) , uint8(M0_disp_img) + uint8(maskVein)*255 ));
imwrite(vesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMap.png')), 'png');
imwrite(VesselImageRGB_Artery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMapArtery.png')), 'png');
imwrite(mat2gray(maskCRA), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRA.png')), 'png');
imwrite(mat2gray(maskCRV), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRV.png')), 'png');

%% new masks

labeled = bwlabel(and(maskArtery, not(circleMask1)));
imwrite(rescale(labeled), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskLabeled.png')), 'png');
imwrite(bwskel(maskArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSkeletonArtery.png')), 'png');
imwrite(bwskel(maskVein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSkeletonVein.png')), 'png');

%% Saving AVI

% w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_RGVideoArtery.avi')));
% tmp = mat2gray(RGVideoArtery);
% open(w)
% for j = 1:size(RGVideoArtery,3)
%     writeVideo(w,tmp(:,:,j)) ;
% end
% close(w);
%
% w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_RGVideoVein.avi')));
% tmp = mat2gray(RGVideoVein);
% open(w)
% for j = 1:size(RGVideoVein,3)
%     writeVideo(w,tmp(:,:,j)) ;
% end
% close(w);
%
% w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_RGVideoVessel.avi')));
% tmp = mat2gray(rgVideoVessel);
% open(w)
%
% for j = 1:size(rgVideoVessel, 3)
%     writeVideo(w, tmp(:, :, j));
% end
%
% close(w);

end
