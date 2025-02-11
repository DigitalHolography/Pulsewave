function [maskArtery, maskVein, maskSection, maskNeighbors , xy_barycenter] = createMasks(M0_ff_video, f_AVG_mean)

ToolBox = getGlobalToolBox;
PW_params = ToolBox.getParams;

if ~isfolder(fullfile(ToolBox.PW_path_png, 'mask'))
    mkdir(ToolBox.PW_path_png, 'mask')
    mkdir(fullfile(ToolBox.PW_path_png, 'mask'), 'steps')
    mkdir(fullfile(ToolBox.PW_path_eps, 'mask'), 'steps')
end
folder_steps = fullfile('mask', 'steps');

%% 0) Test the input for usual cases of wrong aquisition data

if max(f_AVG_mean, [], 'all') <=0
    figure(3); imshow(rescale(f_AVG_mean));
    error("Measurement error from the Moment 1 input. Check the input or re-do the measurement.")
end

%% 1) First Masks and Correlation

[numX, numY, numFrames] = size(M0_ff_video);

M0_ff_img = squeeze(mean(M0_ff_video, 3));
M0_ff_video_centered = M0_ff_video - M0_ff_img;
saveImage(M0_ff_img, ToolBox, 'all_10_M0.png', isStep = true)
if ~isfile(fullfile(ToolBox.PW_path_gif, sprintf("%s_M0.gif", ToolBox.PW_folder_name)))
    writeGifOnDisc(rescale(M0_ff_video), "M0")
end

maskDiaphragm = diskMask(numX, numY, PW_params.masks_diaphragmRadius);
saveImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, ToolBox, 'all_11_maskDiaphragm.png', isStep = true)

% 1) 1) Compute vesselness response

[maskVesselness] = compute_vesselness(M0_ff_img, maskDiaphragm, 'all_12', ToolBox);

% 1) 2) Compute the barycenters and the circle mask

if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
    x_CRV = numX / 2;
    y_CRV = numY / 2;
else
    saveImage(f_AVG_mean, ToolBox, 'all_13_fAVG.png', isStep = true)
    if PW_params.gauss_filt_size_for_barycenter ~= 0
        averaged_fAVG = imgaussfilt(f_AVG_mean, PW_params.gauss_filt_size_for_barycenter, 'Padding', 0) .* maskDiaphragm;
    else
        averaged_fAVG = f_AVG_mean;
    end
    [y_CRA, x_CRA] = find(averaged_fAVG == max(averaged_fAVG, [], 'all'));
    [y_CRV, x_CRV] = find(averaged_fAVG == min(averaged_fAVG, [], 'all'));
end

maskCircle = diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRA/numX, y_CRA/numY]);
maskCircle = maskCircle | diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRV/numX, y_CRV/numY]);

maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness|maskCircle, 1, 8);
saveImage(maskVesselnessClean, ToolBox,  'all_14_connected_mask.png', isStep = true)

%  1) 3) Compute first correlation

cVascular = [0 0 0];
% compute signal in 3 dimentions for correlation in all vessels
vascularSignal = sum(M0_ff_video .* maskVesselnessClean, [1 2]);
vascularSignal = vascularSignal ./ nnz(maskVesselnessClean);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
tLabel = 'Time(s)';
yLabel = 'Power Doppler (a.u.)';

graphSignal('all_15_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, Title = 'Vascular Signal', xlabel = tLabel, ylabel = yLabel);

% compute local-to-average signal wave zero-lag correlation
vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);
R_VascularSignal = mean(M0_ff_video_centered .* vascularSignal_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularSignal_centered, [], 3));
saveImage(R_VascularSignal, ToolBox, 'all_15_Correlation.png', isStep = true)

% 1) 4) Segment Vessels

cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

cmapArtery = [0 0 0; cArtery];
cmapVein = [0 0 0; cVein];
cmapVessels = [cVein; cArtery];

if PW_params.masks_vascular_threshold >= -1 && PW_params.masks_vascular_threshold <= 1
    % IF Manual Thresholds have been set between -1 and 1 then they are used

    firstMaskArtery = (R_VascularSignal > PW_params.masks_vascular_threshold) .* maskVesselnessClean;
    firstMaskVein = (R_VascularSignal < PW_params.masks_vascular_threshold) .* maskVesselnessClean;
    graphThreshHistogram(R_VascularSignal, PW_params.masks_vascular_threshold, maskVesselnessClean, cmapVessels, 'all_16')

else
    % ELSE automatic Otsu segmentation is performed
    % Number of classes for Vessels: 4
    % 1 & 2 = Veins & CoroidalVessels, 3 = CoroidalVessel, 4 = Arteries
    vascularClasses = PW_params.masks_vascular_classes;
    [firstMaskArtery, firstMaskVein] = autoOtsuThresholding(R_VascularSignal, maskVesselnessClean, vascularClasses, 'all_16');
end

saveImage(firstMaskArtery, ToolBox, 'artery_17_FirstMask.png', isStep = true, cmap = cmapArtery)
saveImage(firstMaskVein, ToolBox, 'vein_17_FirstMask.png', isStep = true, cmap = cmapVein)

% Remove small blobs
firstMaskArteryClean = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVeinClean = bwareaopen(firstMaskVein, PW_params.masks_minSize);

saveImage(firstMaskArteryClean, ToolBox,  'artery_18_FirstMaskClean.png', isStep = true, cmap = cmapArtery)
saveImage(firstMaskVeinClean, ToolBox,  'vein_18_FirstMaskClean.png', isStep = true, cmap = cmapVein)

RGBM0(:, :, 1) = rescale(M0_ff_img) + firstMaskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + firstMaskVeinClean;
saveImage(RGBM0, ToolBox,  'all_19_RGB.png', isStep = true)

%% 2)  Improvements of the first mask

diasysAnalysis = PW_params.params.CreationOfMasks.DiaSysAnalysis;

% Parameters for arteries
arteryParams.threshold = PW_params.masks_arterial_threshold;
arteryParams.classes = PW_params.masks_arterial_classes;

% Parameters for veins
veinParams.threshold = PW_params.masks_venous_threshold;
veinParams.classes = PW_params.masks_venous_classes;

results = cell(2, 1);

if diasysAnalysis % Systole/Diastole Analysis

    % 2) 0) Computation of the M0 in Diastole and in Systole

    [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_ff_video, firstMaskArteryClean);

    % 2) 1) New Vesselness Mask

    Systole_Vesselness = compute_vesselness(M0_Systole_img, maskDiaphragm, 'artery_20', ToolBox);
    Diastole_Vesselness = compute_vesselness(M0_Diastole_img, maskDiaphragm, 'artery_20', ToolBox);
    saveImage(rescale(M0_Systole_img), ToolBox,  'artery_20_systole_img.png', isStep = true)
    saveImage(rescale(M0_Diastole_img), ToolBox,  'vein_20_diastole_img.png', isStep = true)

    maskVesselness = Systole_Vesselness | Diastole_Vesselness;
    maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness|maskCircle, 1, 8);
    saveImage(maskVesselnessClean, ToolBox,  'all_20_new_vesselMask.png', isStep = true)

    % 2) 2) Diastole-Systole Image

    diasysArtery = M0_Systole_img - M0_Diastole_img;
    diasysVein = M0_ff_img - diasysArtery;
    diasys = diasysArtery - diasysVein;
    saveImage(diasys, ToolBox,  'all_20_diasys_img.png', isStep = true)
    saveImage(diasysArtery, ToolBox,  'artery_20_diasys_img.png', isStep = true)
    saveImage(diasysVein, ToolBox,  'vein_20_diasys_img.png', isStep = true)

    L = 100 .* rescale(M0_ff_img);
    a = 256 .* rescale(diasys) - 128;
    b = 256 .* rescale(diasys) - 128;

    Labdiasys(:, :, 1) = L;
    Labdiasys(:, :, 2) = a;
    Labdiasys(:, :, 3) = b;
    RGBdiasys = lab2rgb(Labdiasys);
    saveImage(RGBdiasys, ToolBox, 'vessel_40_diasys_rgb.png', isStep = true)
    saveImage(RGBdiasys, ToolBox, 'DiaSysRGB.png')

    % 2) 3) Diastole-Systole based Segmentation

    for  i = 1:2
        if i == 1
            maskArtery = processDiaSysSignal(diasysArtery, maskVesselnessClean, arteryParams, cmapArtery, 'artery_23');
        else
            [~, maskVein] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23');
        end
    end

else % Second Correlation Analysis

    % 2) 4) Artery-Vein correlation based Segmentation

    for  i = 1:2
        if i == 1
            % Process arteries
            results{i} = processVascularSignal(M0_ff_video, firstMaskArteryClean, maskVesselnessClean, arteryParams, cmapArtery, 'artery_22', ToolBox);
        else
            % Process veins
            results{i} = processVascularSignal(M0_ff_video, firstMaskVeinClean, maskVesselnessClean, veinParams, cmapVein, 'vein_22', ToolBox);
        end
    end

    maskArtery = results{1};
    maskVein = results{2};

end

%% 3) Mask Clearing

% 3) 0) Morphological Operations
results = cell(2, 1);

parfor i = 1:2
    if i == 1
        % Process artery mask
        results{i} = clearMasks(maskArtery, 'artery_30', cmapArtery, ToolBox);
    else
        % Process vein mask
        results{i} = clearMasks(maskVein, 'vein_30', cmapVein, ToolBox);
    end
end

maskArteryCleared = results{1};
maskVeinCleared = results{2};

% 3) 1) Final Blob removal

maskArteryNoBlob = maskArteryCleared & bwareafilt(maskArteryCleared|maskCircle, 1, 8);
maskVeinNoBlob = maskVeinCleared & bwareafilt(maskVeinCleared|maskCircle, 1, 8);
saveImage(RGBdiasys, ToolBox, 'vessel_40_diasys_rgb.png', isStep = true)
saveImage(RGBdiasys, ToolBox, 'vessel_40_diasys_rgb.png', isStep = true)

% 3) 2) Force Create Masks in case they exist

if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png'))
    maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
else
    maskArtery = maskArteryNoBlob;
end
if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png'))
    maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;
else
    maskVein = maskVeinNoBlob;
end

% 3) 3) Segmention force width

if PW_params.masks_force_width > 0
    dilationSE = strel('disk', PW_params.masks_force_width);
    maskArteryNoBlob = imdilate(bwskel(maskArteryNoBlob), dilationSE);
    maskVeinNoBlob = imdilate(bwskel(maskVeinNoBlob), dilationSE);

    maskArtery = imdilate(bwskel(maskArtery), dilationSE);
    maskVein = imdilate(bwskel(maskVein), dilationSE);
end

% 3) 4) Segmentation Scores Calculation

segmentationScores(maskArteryNoBlob, maskVeinNoBlob);

% 3) 5) Create Vessel and Background Mask

maskVessel = maskArtery | maskVein;
maskBackground = not(maskVessel);

%% 4) FINAL FIGURES

% 4) 1) RGB Figures
saveImage(maskArtery, ToolBox,'artery_40_Final.png', isStep = true, cmap = cmapArtery)
saveImage(maskVein, ToolBox, 'vein_40_Final.png', isStep = true, cmap = cmapVein)

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArtery;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVein;

saveImage(RGBM0, ToolBox, 'vessel_40_RGB.png', isStep = true)

% 4) 2) Neighbours Mask

veinsAnalysis = PW_params.veins_analysis;
if veinsAnalysis
    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width)) - (maskArtery | maskVein);
else
    maskNeighbors = imdilate(maskArtery, strel('disk', PW_params.local_background_width)) - maskArtery;
end

% 4) 3) CRA and CRV Masks

f_AVG_std = std2(f_AVG_mean);
maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);

saveImage(maskCRA, ToolBox, 'maskCRA.png')
saveImage(maskCRV, ToolBox, 'maskCRV.png')

% 4) 4) Save all images

saveImage(maskArtery, ToolBox, 'maskArtery.png')
saveImage(maskVein, ToolBox, 'maskVein.png')
saveImage(maskVessel, ToolBox, 'maskVessel.png')
saveImage(maskNeighbors, ToolBox, 'maskNeighbors.png')
saveImage(maskBackground, ToolBox, 'maskBackground.png')
saveImage(bwskel(maskArtery), ToolBox, 'skeletonArtery.png')
saveImage(bwskel(maskVein), ToolBox, 'skeletonVein.png')

% 4) 5) Mask Section & Force Barycenter

L = (numY + numX) / 2;
r1 = PW_params.velocitySmallRadiusRatio * L;
r2 = PW_params.velocityBigRadiusRatio * L;

if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
end

xy_barycenter = [x_CRA, y_CRA];
maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMapArtery', maskArtery,thin=10);
createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMap', maskArtery, maskVein,thin=10);

% 4) 6) Create Segmentation Map

if veinsAnalysis
    segmentationMap(:, :, 1) = rescale(M0_ff_img) .* ~(maskArtery | maskVein) + maskArtery;
    segmentationMap(:, :, 2) = rescale(M0_ff_img) .* ~(maskArtery | maskVein);
    segmentationMap(:, :, 3) = rescale(M0_ff_img) .* ~(maskArtery | maskVein) + maskVein;
    saveImage(segmentationMap, ToolBox, 'Segmentation.png')
else
    segmentationMapArtery(:, :, 1) = rescale(M0_ff_img) .* ~maskArtery + maskArtery;
    segmentationMapArtery(:, :, 2) = rescale(M0_ff_img) .* ~maskArtery;
    segmentationMapArtery(:, :, 3) = rescale(M0_ff_img) .* ~maskArtery;
    saveImage(segmentationMapArtery, ToolBox, 'Segmentation.png')
end

if veinsAnalysis
    neighborsMaskSeg(:, :, 1) = rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors) + maskVein;
    saveImage(neighborsMaskSeg, ToolBox, 'NeighborsOnVessels.png')
else
    neighborsMaskSeg(:, :, 1) = rescale(M0_ff_img) .* ~(maskArtery  | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = rescale(M0_ff_img) .* ~(maskArtery  | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = rescale(M0_ff_img) .* ~(maskArtery | maskNeighbors);
    saveImage(neighborsMaskSeg, ToolBox, 'NeighborsOnVessels.png')
end

end
