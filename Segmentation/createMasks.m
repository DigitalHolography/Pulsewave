 function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection, maskNeighbors] = createMasks(M0_ff_video, f_AVG_video, path, ToolBox)

PW_params = Parameters_json(path);
exportVideos = PW_params.exportVideos;

mkdir(ToolBox.PW_path_png, 'mask')
mkdir(fullfile(ToolBox.PW_path_png, 'mask'), 'steps')
mkdir(fullfile(ToolBox.PW_path_eps, 'mask'), 'steps')
folder = fullfile('mask', 'steps');

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
f_AVG_mean = mean(f_AVG_video, 3);

vascularImage = single(squeeze(M0_ff_img .* f_AVG_mean));
blurred_mask = imgaussfilt(vascularImage, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0);
[ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));
[y_CRV, x_CRV] = find(blurred_mask == min(blurred_mask, [], 'all'));
cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;
cercleMask = cercleMask | sqrt((X - x_CRV) .^ 2 + (Y - y_CRV) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | cercleMask, 1, 4);
imwrite(maskVesselnessClean | cercleMask, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_3_choroidClean.png')))

%%  1) 3) Compute first correlation
cVascular = [0 0 0];
% compute signal in 3 dimentions for correlation in all vessels
vascularSignal = mean(M0_ff_video .* maskVesselnessClean, [1 2]);
vascularSignal = vascularSignal ./ nnz(maskVesselnessClean);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
tLabel = 'Time(s)';
yLabel = 'Power Doppler (a.u.)';

graphSignal(ToolBox, 'vascularSignal', folder, t, squeeze(vascularSignal), '-', cVascular, Title='Venous Signal', xlabel=tLabel, ylabel=yLabel)

% compute local-to-average signal wave zero-lag correlation
vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);
R_VascularSignal = mean(M0_ff_video_centered .* vascularSignal_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularSignal_centered, [], 3));
imwrite(rescale(R_VascularSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_4_Correlation.png')))

if ~isempty(PW_params.forcebarycenter)
    ToolBox.y_barycentre = PW_params.forcebarycenter(1);
    ToolBox.x_barycentre = PW_params.forcebarycenter(2);
    y_CRV = numY / 2;
    x_CRV = numX / 2;
else
    vascularImage = double(M0_ff_img .* f_AVG_mean .* R_VascularSignal);
    blurred_mask = imgaussfilt(vascularImage, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0) .* maskDiaphragm;
    [ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));
end

cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;
cercleMask = cercleMask | sqrt((X - x_CRV) .^ 2 + (Y - y_CRV) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

%% 1) 4) Segment Vessels
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

cmapArtery = [0 0 0; cArtery];
cmapVein = [0 0 0; cVein];

if PW_params.masks_vascular_threshold >= -1 && PW_params.masks_vascular_threshold <= 1
    % IF Manual Thresholds have been set between -1 and 1 then they are used

    firstMaskArtery = (R_VascularSignal > PW_params.masks_vascular_threshold) .* maskVesselnessClean;
    firstMaskVein = (R_VascularSignal < PW_params.masks_vascular_threshold) .* maskVesselnessClean;

else
    % ELSE automatic Otsu segmentation is performed
    % Number of classes for Vessels: 4
    % 1 & 2 = Veins & CoroidalVessels, 3 = CoroidalVessel, 4 = Arteries
    vascularClasses = PW_params.masks_vascular_classes;
    [firstMaskArtery, firstMaskVein] = autoOtsuThresholding(R_VascularSignal, maskVesselnessClean, vascularClasses, 'all_1_5', ToolBox);
end

imwrite(logical(firstMaskArtery), cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_1_5_FirstMask.png')))
imwrite(logical(firstMaskVein), cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_1_5_FirstMask.png')))

% Remove small blobs
% minSize = PW_params.masks_minSize * (numX * numY);
firstMaskArteryClean = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVeinClean = bwareaopen(firstMaskVein, PW_params.masks_minSize);

imwrite(firstMaskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_1_6_FirstMaskClean.png')))
imwrite(firstMaskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_1_6_FirstMaskClean.png')))

clear quantizedVesselCorrelation firstThresholds;

%% 1) 5) Choroid Segmentations
cChoroid = [0 179 0] / 255;
cmapChoroid = [0 0 0; cChoroid];

firstMaskChoroid = (~maskVesselnessClean & maskVesselness) & ~firstMaskArtery & ~firstMaskVein;
imwrite(firstMaskChoroid, cmapChoroid, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroid_1_6_FirstMask.png')))

RGBM0(:, :, 1) = rescale(M0_ff_img) + firstMaskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img) + firstMaskChoroid;
RGBM0(:, :, 3) = rescale(M0_ff_img) + firstMaskVeinClean;
imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_1_6_RGB.png')))

%% 2)  Improvements of the first mask

%% 2) 1) Compute a new signal

% compute signal in 3 dimentions for correlation in main arteries
arterialSignal = mean(M0_ff_video .* firstMaskArteryClean .* maskDiaphragm, [1 2]);
arterialSignal = arterialSignal ./ nnz(firstMaskArteryClean .* maskDiaphragm);
arterialSignal_centered = arterialSignal - mean(arterialSignal, 3);

graphSignal(ToolBox, 'arterialSignal', folder, t, squeeze(arterialSignal), '-', cArtery, Title='Arterial Signal', xlabel=tLabel, ylabel=yLabel)

% compute signal in 3 dimentions for correlation in main veins
venousSignal = mean(M0_ff_video .* firstMaskVeinClean, [1 2]);
venousSignal = venousSignal ./ nnz(firstMaskVeinClean);
venousSignal_centered = venousSignal - mean(venousSignal, 3);

graphSignal(ToolBox, 'venousSignal', folder, t, squeeze(venousSignal), '-', cVein, Title='Venous Signal', xlabel=tLabel, ylabel=yLabel)

% compute signal in 3 dimentions for correlation in choroids

choroidalSignal = mean(M0_ff_video .* firstMaskChoroid, [1 2]);
choroidalSignal = choroidalSignal ./ nnz(firstMaskChoroid);
choroidalSignal_centered = choroidalSignal - mean(choroidalSignal, 3);

graphSignal(ToolBox, 'choroidalSignal', folder, t, squeeze(choroidalSignal), '-', cChoroid, Title='Choroidal Signal', xlabel=tLabel, ylabel=yLabel)


graphSignal(ToolBox, 'everySignal', folder, ...
    t, squeeze(vascularSignal), '--', cVascular, ...
    t, squeeze(arterialSignal), '-', cArtery, ...
    t, squeeze(venousSignal), '-', cVein, ...
    Title='Signals used for correlation maps', xlabel = tLabel, ylabel = yLabel, Legends={'Vessel','Artery','Vein'})

clear arterialSignal venousSignal choroidalSignal vascularSignal

%% 2) 2) New Correlations

% compute local-to-average signal wave zero-lag correlation
R_ArterialSignal = mean((M0_ff_video_centered .* arterialSignal_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(arterialSignal_centered, [], 3));
R_VenousSignal = mean((M0_ff_video_centered .* venousSignal_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(arterialSignal_centered, [], 3));
R_ChoroidalSignal = mean((M0_ff_video_centered .* choroidalSignal_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(choroidalSignal_centered, [], 3));

imwrite(rescale(R_ArterialSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_1_CorrelMatrix.png')))
imwrite(rescale(R_VenousSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_1_CorrelMatrix.png')))
imwrite(rescale(R_ChoroidalSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroid_2_1_CorrelMatrix.png')))

clear arterialSignal_centered venousSignal_centered choroidalSignal_centered vascularSignal_centered
correlation_Vesselness_flag = PW_params.masks_correlation_vesselness_flag;

if correlation_Vesselness_flag
    maskSegmentationClean = maskDiaphragm;

    R_ArteryVessel = imgaussfilt(R_ArterialSignal) .* maskSegmentationClean;
    graphThreshHistogram(R_ArteryVessel, 0, maskDiaphragm, [0 0 0; cArtery], 'artery_2_2_R_Histo', ToolBox)
    R_ArteryVessel = R_ArteryVessel .* (R_ArteryVessel > 0);

    R_ChoroidVessel = R_ChoroidalSignal .* maskDiaphragm;
    graphThreshHistogram(R_ChoroidVessel, 0, maskDiaphragm, [0 0 0; cChoroid], 'choroid_2_2_R_Histo', ToolBox)

    R_VeinVessel = imgaussfilt(R_ChoroidalSignal - R_ArterialSignal) .* maskSegmentationClean;
    graphThreshHistogram(R_VeinVessel, 0, maskDiaphragm, [0 0 0; cVein], 'vein_2_2_R_Histo', ToolBox)
    R_VeinVessel = R_VeinVessel .* (R_VeinVessel > 0);

else
    maskSegmentationClean = maskVesselnessClean;

    R_ArteryVessel = R_ArterialSignal .* maskSegmentationClean;
    graphThreshHistogram(R_ArteryVessel, 0, maskDiaphragm, [0 0 0; cArtery], 'artery_2_2_R_Histo', ToolBox)
    R_ArteryVessel = R_ArteryVessel .* (R_ArteryVessel > 0);

    R_VeinVessel = R_VenousSignal .* maskSegmentationClean;
    graphThreshHistogram(R_VeinVessel, 0, maskDiaphragm, [0 0 0; cVein], 'vein_2_2_R_Histo', ToolBox)
    R_VeinVessel = R_VeinVessel .* (R_VeinVessel > 0);

end

imwrite(rescale(R_ArteryVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_2_R_ArteryVessel.png')))
imwrite(rescale(R_VeinVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_2_R_VeinVessel.png')))
if correlation_Vesselness_flag
    imwrite(rescale(R_ChoroidVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'choroid_2_2_R_ChoroidalVessel.png')))
end

%% 2) 3) Automatic or Manual thresholds for arteries and veins

if PW_params.masks_arterial_threshold >= -1 && PW_params.masks_arterial_threshold <= 1
    % Manual Threshold
    maskArtery = R_ArteryVessel >= PW_params.masks_arterial_threshold;

else
    % Automatic Otsu Threshold
    arterialClasses = PW_params.masks_arterial_classes;
    maskArtery = autoOtsuThresholding(R_ArteryVessel, maskSegmentationClean, arterialClasses, 'artery_2_3', ToolBox);
    maskArtery = maskArtery | firstMaskArteryClean;

end

imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_2_4_Thresh.png')))

if PW_params.masks_venous_threshold >= -1 && PW_params.masks_venous_threshold <= 1
    % Manual Threshold
    maskVein = R_VeinVessel >= PW_params.masks_venous_threshold;

else
    % Automatic Otsu Threshold
    venousClasses = PW_params.masks_venous_classes;
    [~, maskVein] = autoOtsuThresholding(R_VeinVessel, maskSegmentationClean, venousClasses, 'vein_2_3', ToolBox);
    maskVein = maskVein | firstMaskVeinClean;

end

imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_2_4_Thresh.png')))

%% Clearing 3)
%% 3) 1) Morphological Operations
% Remove the small bits
maskArteryClean = bwareaopen(maskArtery, PW_params.masks_minSize, 4);
maskVeinClean = bwareaopen(maskVein, PW_params.masks_minSize, 4);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_1_AreaOpen.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_1_AreaOpen.png')))

% Remove the small gaps
imcloseSE = strel('disk', PW_params.masks_imclose_radius);
maskArteryClean = imclose(maskArteryClean, imcloseSE);
maskVeinClean = imclose(maskVeinClean, imcloseSE);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_2_ImClose.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_2_ImClose.png')))

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVeinClean;

imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_2_RGB.png')))
%% 3) 2) Minimum Mask Width
minWidthSE = strel('disk', PW_params.masks_min_width);

skelArtery = bwskel(maskArteryClean);
skelVein = bwskel(maskVeinClean);

maskArteryClean = maskArteryClean | imdilate(skelArtery, minWidthSE);
maskVeinClean = maskVeinClean | imdilate(skelVein, minWidthSE);

clear skelArtery skelVein
%% 3) Bonus) Force Create Masks
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')) && isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png'))
    fprintf("FORCED MASK\n")
    maskArteryClean = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
    maskVeinClean = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;
end

%% 3) 3) Create Vessel and Background Mask
maskArtery = maskArteryClean;
maskVein = maskVeinClean;

maskVessel = maskArtery | maskVein;
maskBackground = not(maskVessel);

if PW_params.veins_analysis
    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width)) - (maskArtery | maskVein);
else
    maskNeighbors = imdilate(maskArtery , strel('disk', PW_params.local_background_width)) - maskArtery;
end

imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'artery_3_3_Final.png')))
imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'vein_3_3_Final.png')))
imwrite(maskVessel, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_3_Final.png')))
imwrite(maskNeighbors, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskNeighbors.png')));
imwrite(rescale(uint8(cat(3, uint8(M0_ff_img) + uint8(maskArteryClean) * 255, uint8(M0_ff_img) + uint8(cercleMask), uint8(M0_ff_img) + uint8(maskVeinClean) * 255))), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", ToolBox.main_foldername, 'all_3_4_Final.png')))

%% SENSITIVITY & SPECIFICITY
% Add targetMaskArtery.png and targetMaskVein.png in a mask folder inside
% the pulsewave folder to have the specificity and the senisivity of the
% algorithm used with your parameters.
try
    targetMaskArtery = imread(fullfile(ToolBox.PW_path_main, 'mask', 'targetMaskArtery.png'));
    targetMaskVein = imread(fullfile(ToolBox.PW_path_main, 'mask', 'targetMaskVein.png'));

    TPArtery = nnz(targetMaskArtery & maskArteryClean);
    FPArtery = nnz(~targetMaskArtery & maskArteryClean);
    FNArtery = nnz(targetMaskArtery & ~maskArteryClean);
    TNArtery = nnz(~targetMaskArtery & ~maskArteryClean);

    TPVein = nnz(targetMaskVein & maskVeinClean);
    FPVein = nnz(~targetMaskVein & maskVeinClean);
    FNVein = nnz(targetMaskVein & ~maskVeinClean);
    TNVein = nnz(~targetMaskVein & ~maskVeinClean);

    fileID = fopen(fullfile(ToolBox.PW_path_log, sprintf("%s_confusionMatrix.txt", ToolBox.main_foldername)), 'w');
    fprintf(fileID, "Artery Sensitivity: %0.2f %%\n", 100 * TPArtery / (TPArtery + FNArtery));
    fprintf(fileID, "Artery Specificity: %0.2f %%\n", 100 * TNArtery / (TNArtery + FPArtery));
    fprintf(fileID, "Vein Sensitivity: %0.2f %%\n", 100 * TPVein / (TPVein + FNVein));
    fprintf(fileID, "Vein Specificity: %0.2f %%\n", 100 * TNVein / (TNVein + FPVein));
    fclose(fileID);
catch
    fprintf("No target masks to perform sensitivity & specificity evaluation\n")
end

%% Force Create Masks
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')) && isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png'))
    maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
    maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;
    maskVessel = maskArtery | maskVein;
    maskBackground = not(maskVessel);
end

%% Force Barycenter
if ~isempty(PW_params.forcebarycenter)
    ToolBox.y_barycentre = PW_params.forcebarycenter(1);
    ToolBox.x_barycentre = PW_params.forcebarycenter(2);
    cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;
end

%% Create CRA and CRV Mask
f_AVG_std = std2(f_AVG_mean);
maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);
clear f_AVG_std f_AVG_mean

%% Create Mask Section
radius1 = PW_params.velocitySmallRadiusRatio * (numY + numX) / 2;
radius2 = PW_params.velocityBigRadiusRatio * (numY + numX) / 2;

circleMask1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius1;
circleMask2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius2;

maskSection = xor(circleMask1, circleMask2);

%% Create Colormap Artery/Vein

M0_ff_img = mat2gray(M0_ff_img);
M0_ff_img_RGB = cat(3, M0_ff_img, M0_ff_img, M0_ff_img);
maskSectionRGB = cat(3, maskSection, maskSection, maskSection);

[hueA, satA, valA] = createHSVmap(M0_ff_img, ~maskSection & maskArtery, 0, 0);
[hueV, satV, valV] = createHSVmap(M0_ff_img, ~maskSection & maskVein, 0.7, 0.7);
[hueSectionA, satSectionA, valSectionA] = createHSVmap(M0_ff_img, maskSection & maskArtery, 0.15, 0.15);
[hueSectionV, satSectionV, valSectionV] = createHSVmap(M0_ff_img, maskSection & maskVein, 0.5, 0.5);

hueVessel = hueA .* (~maskSection & maskArtery) + hueV .* ( ~maskSection & maskVein) + hueSectionA  .* (maskSection & maskArtery) + hueSectionV .* (maskSection & maskVein);
satVessel = satA .* (~maskSection & maskArtery) + satV .* ( ~maskSection & maskVein) + satSectionA  .* (maskSection & maskArtery) + satSectionV .* (maskSection & maskVein);
valVessel = valA .* (~maskSection & maskArtery) + valV .* ( ~maskSection & maskVein) + valSectionA  .* (maskSection & maskArtery) + valSectionV .* (maskSection & maskVein);

vesselImageRGB = hsv2rgb(hueVessel, satVessel, valVessel) + maskSectionRGB .* ~maskVessel + (M0_ff_img_RGB .* ~maskVessel);

%% Create Colormap Artery Only

hueArtery = hueA .* (~maskSection & maskArtery) + hueSectionA .* (maskSection & maskArtery);
satArtery = satA .* (~maskSection & maskArtery) + satSectionA .* (maskSection & maskArtery);
valArtery = valA .* (~maskSection & maskArtery) + valSectionA .* (maskSection & maskArtery);

VesselImageRGB_Artery = hsv2rgb(hueArtery, satArtery, valArtery) + maskSectionRGB .* ~maskArtery + (M0_ff_img_RGB .* ~maskArtery);

%% Create Segmentation Map

segmentationMap = zeros(numX, numY, 3);
segmentationMap(:, :, 1) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img + maskArtery;
segmentationMap(:, :, 2) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img;
segmentationMap(:, :, 3) = M0_ff_img - (maskArtery + maskVein) .* M0_ff_img + maskVein;
segmentationMapArtery(:, :, 1) = M0_ff_img - (maskArtery) .* M0_ff_img + maskArtery;
segmentationMapArtery(:, :, 2) = M0_ff_img - maskArtery .* M0_ff_img;
segmentationMapArtery(:, :, 3) = M0_ff_img - maskArtery .* M0_ff_img;
figure, imagesc(segmentationMapArtery), axis image
imwrite(segmentationMap, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arteryVeinSegmentation.png')), 'png');
imwrite(segmentationMapArtery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arterySegmentation.png')), 'png');

%% Saving masks as PNG
foldername = ToolBox.main_foldername;

imwrite(mat2gray(double(vesselnessM0)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselness.png')), 'png');
imwrite(mat2gray(single(maskArtery)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskArtery.png')), 'png');
imwrite(mat2gray(single(maskVein)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVein.png')), 'png');
imwrite(mat2gray(single(maskVessel)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVessel.png')), 'png');
imwrite(mat2gray(single(maskBackground)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskBackground.png')), 'png');
imwrite(vesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMap.png')), 'png');
imwrite(VesselImageRGB_Artery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMapArtery.png')), 'png');
imwrite(mat2gray(maskCRA), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRA.png')), 'png');
imwrite(mat2gray(maskCRV), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRV.png')), 'png');

%% new masks

labeled = bwlabel(and(maskArtery, not(cercleMask)));
imwrite(rescale(labeled), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskLabeled.png')), 'png');
imwrite(bwskel(maskArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSkeletonArtery.png')), 'png');
imwrite(bwskel(maskVein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSkeletonVein.png')), 'png');

end
