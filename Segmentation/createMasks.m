function [maskArtery, maskVein, maskSection, maskNeighbors , xy_barycenter] = createMasks(M0_ff_video, f_AVG_mean)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);

if ~isfolder(fullfile(ToolBox.PW_path_png, 'mask'))
    mkdir(ToolBox.PW_path_png, 'mask')
    mkdir(fullfile(ToolBox.PW_path_png, 'mask'), 'steps')
    mkdir(fullfile(ToolBox.PW_path_eps, 'mask'), 'steps')
end
folder_steps = fullfile('mask', 'steps');
main_folder = ToolBox.main_foldername;

%% 0) Test the input for usual cases of wrong aquisition data

if max(f_AVG_mean, [], 'all') <=0
    figure(3); imshow(rescale(f_AVG_mean));
    error("Measurement error from the Moment 1 input. Check the input or re-do the measurement.")
end

%% 1) First Masks and Correlation

[numX, numY, numFrames] = size(M0_ff_video);
M0_ff_video = rescale(M0_ff_video);
maskDiaphragm = diskMask(numX, numY, PW_params.masks_diaphragmRadius);
imwrite(rescale(maskDiaphragm), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_0_maskDiaphragm.png')))

M0_ff_img = squeeze(mean(M0_ff_video, 3));
M0_ff_video_centered = M0_ff_video - M0_ff_img;
imwrite(rescale(M0_ff_img), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_0_M0.png')))
if ~isfile(fullfile(ToolBox.PW_path_gif, sprintf("%s_M0.gif", ToolBox.PW_folder_name)))
    writeGifOnDisc(M0_ff_video, "M0")
end

%% 1) 1) Compute vesselness response
vesselnessM0 = vesselness_filter(M0_ff_img, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);
maskVesselness = logical(imbinarize(vesselnessM0 .* maskDiaphragm));

imwrite(rescale(vesselnessM0), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_1_Vesselness.png')))
imwrite(maskVesselness, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_2_vesselMask.png')))

%% 1) 2) Compute the barycenters and the circle mask

if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
    x_CRV = numX / 2;
    y_CRV = numY / 2;
else
    imwrite(rescale(f_AVG_mean), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_2_fAVG.png')))
    if PW_params.gauss_filt_size_for_barycenter ~= 0
        averaged_fAVG = imgaussfilt(f_AVG_mean, PW_params.gauss_filt_size_for_barycenter * numX, 'Padding', 0) .* maskDiaphragm;
    else
        averaged_fAVG = f_AVG_mean;
    end
    [y_CRA, x_CRA] = find(averaged_fAVG == max(averaged_fAVG, [], 'all'));
    [y_CRV, x_CRV] = find(averaged_fAVG == min(averaged_fAVG, [], 'all'));
end

xy_barycenter = [x_CRA, y_CRA];
maskCircle = diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRA/numX, y_CRA/numY]);
maskCircle = maskCircle | diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRV/numX, y_CRV/numY]);

maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 4);
imwrite(maskVesselnessClean | maskCircle, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_3_choroidClean.png')))

%%  1) 3) Compute first correlation
cVascular = [0 0 0];
% compute signal in 3 dimentions for correlation in all vessels
vascularSignal = sum(M0_ff_video .* maskVesselnessClean, [1 2]);
vascularSignal = vascularSignal ./ nnz(maskVesselnessClean);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
tLabel = 'Time(s)';
yLabel = 'Power Doppler (a.u.)';

graphSignal('all_1_3_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, Title = 'Venous Signal', xlabel = tLabel, ylabel = yLabel);

% compute local-to-average signal wave zero-lag correlation
vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);
R_VascularSignal = mean(M0_ff_video_centered .* vascularSignal_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularSignal_centered, [], 3));
imwrite(rescale(R_VascularSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_4_Correlation.png')))

%% 1) 4) Segment Vessels
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

cmapArtery = [0 0 0; cArtery];
cmapVein = [0 0 0; cVein];
cmapVessels = [cVein; cArtery];

if PW_params.masks_vascular_threshold >= -1 && PW_params.masks_vascular_threshold <= 1
    % IF Manual Thresholds have been set between -1 and 1 then they are used

    firstMaskArtery = (R_VascularSignal > PW_params.masks_vascular_threshold) .* maskVesselnessClean;
    firstMaskVein = (R_VascularSignal < PW_params.masks_vascular_threshold) .* maskVesselnessClean;
    graphThreshHistogram(R_VascularSignal, PW_params.masks_vascular_threshold, maskVesselnessClean, cmapVessels, 'all_1_5')

else
    % ELSE automatic Otsu segmentation is performed
    % Number of classes for Vessels: 4
    % 1 & 2 = Veins & CoroidalVessels, 3 = CoroidalVessel, 4 = Arteries
    vascularClasses = PW_params.masks_vascular_classes;
    [firstMaskArtery, firstMaskVein] = autoOtsuThresholding(R_VascularSignal, maskVesselnessClean, vascularClasses, 'all_1_5');
end

imwrite(logical(firstMaskArtery), cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_1_5_FirstMask.png')))
imwrite(logical(firstMaskVein), cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_1_5_FirstMask.png')))

% Remove small blobs
% minSize = PW_params.masks_minSize * (numX * numY);
firstMaskArteryClean = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVeinClean = bwareaopen(firstMaskVein, PW_params.masks_minSize);

imwrite(firstMaskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_1_6_FirstMaskClean.png')))
imwrite(firstMaskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_1_6_FirstMaskClean.png')))

clear quantizedVesselCorrelation firstThresholds;

RGBM0(:, :, 1) = rescale(M0_ff_img) + firstMaskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + firstMaskVeinClean;
imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_1_6_RGB.png')))

%% 2)  Improvements of the first mask

%% 2) 1) Compute a new signal

% compute signal in 3 dimentions for correlation in main arteries
arterialSignal = sum(M0_ff_video .* firstMaskArteryClean, [1 2]);
arterialSignal = arterialSignal ./ nnz(firstMaskArteryClean);
arterialSignal_centered = arterialSignal - mean(arterialSignal, 3);

graphSignal('arterialSignal', folder_steps, t, squeeze(arterialSignal), '-', cArtery, Title = 'Arterial Signal', xlabel = tLabel, ylabel = yLabel);

% compute signal in 3 dimentions for correlation in main veins
venousSignal = sum(M0_ff_video .* firstMaskVeinClean, [1 2]);
venousSignal = venousSignal ./ nnz(firstMaskVeinClean);
venousSignal_centered = venousSignal - mean(venousSignal, 3);

graphSignal('venousSignal', folder_steps, t, squeeze(venousSignal), '-', cVein, Title = 'Venous Signal', xlabel = tLabel, ylabel = yLabel);

graphSignal('everySignal', folder_steps, ...
    t, squeeze(vascularSignal), '--', cVascular, ...
    t, squeeze(arterialSignal), '-', cArtery, ...
    t, squeeze(venousSignal), '-', cVein, ...
    Title = 'Signals used for correlation maps', xlabel = tLabel, ylabel = yLabel, Legends = {'Vessel', 'Artery', 'Vein'});

clear arterialSignal venousSignal vascularSignal

%% 2) 2) New Correlations

% compute local-to-average signal wave zero-lag correlation
R_ArterialSignal = mean((M0_ff_video_centered .* arterialSignal_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(arterialSignal_centered, [], 3));
R_VenousSignal = mean((M0_ff_video_centered .* venousSignal_centered), 3) ./ (std((M0_ff_video_centered), [], 3) * std(venousSignal_centered, [], 3));

imwrite(rescale(R_ArterialSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_2_1_CorrelMatrix.png')))
imwrite(rescale(R_VenousSignal), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_2_1_CorrelMatrix.png')))

clear arterialSignal_centered venousSignal_centered vascularSignal_centered

R_ArteryVessel = R_ArterialSignal .* maskVesselnessClean;
graphThreshHistogram(R_ArteryVessel, 0, maskDiaphragm, [0 0 0; cArtery], 'artery_2_2_R_Histo')
R_ArteryVessel = R_ArteryVessel .* (R_ArteryVessel > 0);

R_VeinVessel = R_VenousSignal .* maskVesselnessClean;
graphThreshHistogram(R_VeinVessel, 0, maskDiaphragm, [0 0 0; cVein], 'vein_2_2_R_Histo')
R_VeinVessel = R_VeinVessel .* (R_VeinVessel > 0);

imwrite(rescale(R_ArteryVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_2_2_R_ArteryVessel.png')))
imwrite(rescale(R_VeinVessel), fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_2_2_R_VeinVessel.png')))

imwrite(immultiply(R_ArterialSignal, M0_ff_img), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'artery_Corr.png')));
imwrite(immultiply(R_VenousSignal, M0_ff_img), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'vein_Corr.png')));

%% 2) 3) Automatic or Manual thresholds for arteries and veins

if PW_params.masks_arterial_threshold >= -1 && PW_params.masks_arterial_threshold <= 1
    % Manual Threshold
    maskArtery = R_ArteryVessel >= PW_params.masks_arterial_threshold;
    graphThreshHistogram(R_ArteryVessel, PW_params.masks_arterial_threshold, maskArtery, cmapArtery, 'artery_2_3')
else
    % Automatic Otsu Threshold
    arterialClasses = PW_params.masks_arterial_classes;
    maskArtery = autoOtsuThresholding(R_ArteryVessel, maskVesselnessClean, arterialClasses, 'artery_2_3');
    maskArtery = maskArtery | firstMaskArteryClean;

end

imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_2_4_Thresh.png')))

if PW_params.masks_venous_threshold >= -1 && PW_params.masks_venous_threshold <= 1
    % Manual Threshold
    maskVein = R_VeinVessel >= PW_params.masks_venous_threshold;
    graphThreshHistogram(R_VeinVessel, PW_params.masks_venous_threshold, maskVein, cmapVein, 'vein_2_3')

else
    % Automatic Otsu Threshold
    venousClasses = PW_params.masks_venous_classes;
    [~, maskVein] = autoOtsuThresholding(R_VeinVessel, maskVesselnessClean, venousClasses, 'vein_2_3');
    maskVein = maskVein | firstMaskVeinClean;

end

imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_2_4_Thresh.png')))

%% Clearing 3)
%% 3) 1) Morphological Operations
% Remove the small bits
maskArteryClean = bwareaopen(maskArtery, PW_params.masks_minSize, 4);
maskVeinClean = bwareaopen(maskVein, PW_params.masks_minSize, 4);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_3_1_AreaOpen.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_3_1_AreaOpen.png')))

% Remove the small gaps
imcloseSE = strel('disk', PW_params.masks_imclose_radius);
maskArteryClean = imclose(maskArteryClean, imcloseSE);
maskVeinClean = imclose(maskVeinClean, imcloseSE);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_3_2_ImClose.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_3_2_ImClose.png')))

%% 3) 2) Minimum Mask Width
minWidthSE = strel('disk', PW_params.masks_min_width);

skelArtery = bwskel(maskArteryClean);
skelVein = bwskel(maskVeinClean);

maskArteryClean = maskArteryClean | imdilate(skelArtery, minWidthSE);
maskVeinClean = maskVeinClean | imdilate(skelVein, minWidthSE);

clear skelArtery skelVein

%% 3) 3) Final Dilation
dilationSE = strel('disk', PW_params.masks_imdilateFinal);

maskArteryClean = imdilate(maskArteryClean, dilationSE);
maskVeinClean = imdilate(maskVeinClean, dilationSE);

imwrite(maskArteryClean, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_3_3_ImDilate.png')))
imwrite(maskVeinClean, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_3_3_ImDilate.png')))

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVeinClean;

imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'all_3_2_RGB.png')))

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

    fileID = fopen(fullfile(ToolBox.PW_path_log, sprintf("%s_confusionMatrix.txt", main_folder)), 'w');
    fprintf(fileID, "    Artery Sensitivity: %0.2f %%\n", 100 * TPArtery / (TPArtery + FNArtery));
    fprintf(fileID, "    Artery Specificity: %0.2f %%\n", 100 * TNArtery / (TNArtery + FPArtery));
    fprintf(fileID, "    Vein Sensitivity: %0.2f %%\n", 100 * TPVein / (TPVein + FNVein));
    fprintf(fileID, "    Vein Specificity: %0.2f %%\n", 100 * TNVein / (TNVein + FPVein));
    fclose(fileID);
catch
    fprintf("    No target masks to perform sensitivity & specificity evaluation\n")
end

%% 3) 3) Create Vessel and Background Mask

%% Force Create Masks
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png'))
    maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
else
    maskArtery = maskArteryClean;
end
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png'))
    maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;
else
    maskVein = maskVeinClean;
end

maskVessel = maskArtery | maskVein;
maskBackground = not(maskVessel);

%% Create Neighbours Mask
veinsAnalysis = PW_params.veins_analysis;
if veinsAnalysis
    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width)) - (maskArtery | maskVein);
else
    maskNeighbors = imdilate(maskArtery, strel('disk', PW_params.local_background_width)) - maskArtery;
end

%% Create CRA and CRV Mask
f_AVG_std = std2(f_AVG_mean);
maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);
clear f_AVG_std f_AVG_mean

%% Save Images

imwrite(maskArtery, cmapArtery, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'artery_3_3_Final.png')))
imwrite(maskVein, cmapVein, fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, 'vein_3_3_Final.png')))

imwrite(maskArtery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskArtery.png')))
imwrite(maskVein, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskVein.png')))
imwrite(maskVessel, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskVessel.png')))
imwrite(maskNeighbors, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskNeighbors.png')));
imwrite(maskBackground, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskBackground.png')));

if veinsAnalysis
    neighborsMaskSeg(:, :, 1) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskVein;
    imwrite(neighborsMaskSeg, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'NeighborsOnVessels.png')));
else
    neighborsMaskSeg(:, :, 1) = M0_ff_img .* ~(maskArtery  | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = M0_ff_img .* ~(maskArtery  | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = M0_ff_img .* ~(maskArtery | maskNeighbors);
    imwrite(neighborsMaskSeg, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'NeighborsOnVessels.png')));
end

imwrite(maskCRA, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskCRA.png')));
imwrite(maskCRV, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskCRV.png')));

%% Force Barycenter
if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
    maskCircle = diskMask(numX, numY, PW_params.masks_crop_radius, "options", [x_CRA, y_CRA]);
end

%% Create Mask Section
L = (numY + numX) / 2;
r1 = PW_params.velocitySmallRadiusRatio * L;
r2 = PW_params.velocityBigRadiusRatio * L;
maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMapArtery', maskArtery);

createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMap', maskArtery, maskVein);

%% Create Segmentation Map
if veinsAnalysis
    segmentationMap(:, :, 1) = M0_ff_img .* ~(maskArtery | maskVein) + maskArtery;
    segmentationMap(:, :, 2) = M0_ff_img .* ~(maskArtery | maskVein);
    segmentationMap(:, :, 3) = M0_ff_img .* ~(maskArtery | maskVein) + maskVein;
    imwrite(segmentationMap, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'Segmentation.png')), 'png');
else
    segmentationMapArtery(:, :, 1) = M0_ff_img .* ~maskArtery + maskArtery;
    segmentationMapArtery(:, :, 2) = M0_ff_img .* ~maskArtery;
    segmentationMapArtery(:, :, 3) = M0_ff_img .* ~maskArtery;
    imwrite(segmentationMapArtery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'Segmentation.png')), 'png');
end

%% new masks

labeled = bwlabel(and(maskArtery, not(maskCircle)));
imwrite(rescale(labeled), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'maskLabeled.png')), 'png');
imwrite(bwskel(maskArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'skeletonArtery.png')), 'png');
imwrite(bwskel(maskVein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, 'skeletonVein.png')), 'png');

end
