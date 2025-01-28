function [maskArtery, maskVein, maskSection, maskNeighbors , xy_barycenter] = createMasks(M0_ff_video, f_AVG_mean)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);

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
M0_ff_video = rescale(M0_ff_video);
maskDiaphragm = diskMask(numX, numY, PW_params.masks_diaphragmRadius);
saveImage(rescale(maskDiaphragm), ToolBox, 'all_1_0_maskDiaphragm.png', isStep = true)

M0_ff_img = squeeze(mean(M0_ff_video, 3));
M0_ff_video_centered = M0_ff_video - M0_ff_img;
saveImage(rescale(M0_ff_img), ToolBox, 'all_1_0_M0.png', isStep = true)
if ~isfile(fullfile(ToolBox.PW_path_gif, sprintf("%s_M0.gif", ToolBox.PW_folder_name)))
    writeGifOnDisc(M0_ff_video, "M0")
end

%% 1) 1) Compute vesselness response

vesselnessM0 = vesselness_filter(M0_ff_img, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);
maskVesselness = logical(imbinarize(vesselnessM0 .* maskDiaphragm)); % Trivial Thresholding

saveImage(rescale(vesselnessM0), ToolBox, 'all_1_1_Vesselness.png', isStep = true)
saveImage(maskVesselness, ToolBox, 'all_1_2_vesselMask.png', isStep = true)

%% 1) 2) Compute the barycenters and the circle mask

if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
    x_CRV = numX / 2;
    y_CRV = numY / 2;
else
    saveImage(rescale(f_AVG_mean), ToolBox, 'all_1_2_fAVG.png', isStep = true)
    if PW_params.gauss_filt_size_for_barycenter ~= 0
        averaged_fAVG = imgaussfilt(f_AVG_mean, PW_params.gauss_filt_size_for_barycenter * numX, 'Padding', 0) .* maskDiaphragm;
    else
        averaged_fAVG = f_AVG_mean;
    end
    [y_CRA, x_CRA] = find(averaged_fAVG == max(averaged_fAVG, [], 'all'));
    [y_CRV, x_CRV] = find(averaged_fAVG == min(averaged_fAVG, [], 'all'));
end

maskCircle = diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRA/numX, y_CRA/numY]);
maskCircle = maskCircle | diskMask(numX, numY, PW_params.masks_crop_radius, 'center', [x_CRV/numX, y_CRV/numY]);

maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 4);
saveImage(maskVesselnessClean + maskCircle .* 0.5, ToolBox, 'all_1_3_choroidClean.png', isStep = true)

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
saveImage(rescale(R_VascularSignal), ToolBox, 'all_1_4_Correlation.png', isStep = true)

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

saveImage(logical(firstMaskArtery), ToolBox, 'artery_1_5_FirstMask.png', isStep = true, cmap = cmapArtery)
saveImage(logical(firstMaskVein), ToolBox, 'vein_1_5_FirstMask.png', isStep = true, cmap = cmapVein)

% Remove small blobs
% minSize = PW_params.masks_minSize * (numX * numY);
firstMaskArteryClean = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVeinClean = bwareaopen(firstMaskVein, PW_params.masks_minSize);

saveImage(firstMaskArteryClean, ToolBox,  'artery_1_6_FirstMaskClean.png', isStep = true, cmap = cmapArtery)
saveImage(firstMaskVeinClean, ToolBox,  'vein_1_6_FirstMaskClean.png', isStep = true, cmap = cmapVein)

clear quantizedVesselCorrelation firstThresholds;

RGBM0(:, :, 1) = rescale(M0_ff_img) + firstMaskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + firstMaskVeinClean;
saveImage(RGBM0, ToolBox,  'all_1_6_RGB.png', isStep = true)

%% 2)  Improvements of the first mask

%% 2) 1) Compute a new signal
diasysAnalysis = PW_params.params.CreationOfMasks.DiaSysAnalysis;

if diasysAnalysis
    %% 2) 0) Systole/Diastole

    [sys_index_list, ~] = find_systole_index(M0_ff_video, firstMaskArteryClean);

    numSys = numel(sys_index_list); % number of systoles
    fpCycle = round(numFrames / numSys); % Frames per cycle

    sysindexes = [];
    diaindexes = [];

    for idx = 1:numSys
        try
            sysindexes = [sysindexes, (sys_index_list(idx):sys_index_list(idx) + round(fpCycle * 0.1))];
        catch
        end
        try
            diaindexes = [diaindexes, (sys_index_list(idx) - round(fpCycle * 0.15) :sys_index_list(idx) - round(fpCycle * 0.05))];
        catch
        end
    end

    SystoleM0 = mean(M0_ff_video(:, :, sysindexes), 3);
    DiastoleM0 = mean(M0_ff_video(:, :, diaindexes), 3);
    saveImage(SystoleM0, ToolBox,  'all_2_0_M0_Systole.png', isStep = true)
    saveImage(DiastoleM0, ToolBox,  'all_2_0_M0_Diastole.png', isStep = true)

    diasys = SystoleM0 - DiastoleM0;
    saveImage(rescale(diasys), ToolBox,  'all_2_0_DIASYSGRAY.png', isStep = true)

    Labdiasys(:, :, 1) = 100 .* rescale(M0_ff_img);
    Labdiasys(:, :, 2) = 256 .* rescale(diasys) - 128;
    Labdiasys(:, :, 3) = 256 .* rescale(diasys) - 128;
    RGBdiasys = lab2rgb(Labdiasys);
    HSVdiasys = rgb2hsv(RGBdiasys);
    saveImage(RGBdiasys, ToolBox, 'all_2_0_DIASYSRGB.png', isStep = true)

    %% NEW ARTERY MASK

    vesselnessM0 = vesselness_filter(diasys, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);
    maskVesselness = logical(imbinarize(vesselnessM0 .* maskDiaphragm));
    saveImage(maskVesselness, ToolBox, 'artery_2_1_choroid.png', isStep = true)

    maskArtery = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 4);
    saveImage(maskArtery + maskCircle .* 0.5, ToolBox, 'artery_2_1_choroidClean.png', isStep = true)

    %% NEW VEIN MASK

    vesselnessM0 = vesselness_filter(DiastoleM0, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);
    maskVesselness = logical(imbinarize(vesselnessM0 .* maskDiaphragm));
    saveImage(maskVesselness, ToolBox, 'vein_2_1_choroid.png', isStep = true)

    maskVein = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 4);
    maskVein = maskVein & ~maskArtery;
    saveImage(maskVein + maskCircle .* 0.5, ToolBox, 'vein_2_1_choroidClean.png', isStep = true)

else
    results = cell(2, 1);

    % Parameters for arteries
    arteryParams.threshold = PW_params.masks_arterial_threshold;
    arteryParams.classes = PW_params.masks_arterial_classes;

    % Parameters for veins
    veinParams.threshold = PW_params.masks_venous_threshold;
    veinParams.classes = PW_params.masks_venous_classes;

    for  i = 1:2
        if i == 1
            % Process arteries
            results{i} = processVascularSignal(M0_ff_video, firstMaskArteryClean, maskVesselnessClean, arteryParams, cmapArtery, 'artery', main_folder, ToolBox);
        else
            % Process veins
            results{i} = processVascularSignal(M0_ff_video, firstMaskVeinClean, maskVesselnessClean, veinParams, cmapVein, 'vein', main_folder, ToolBox);
        end
    end

    % Assign results to their respective variables
    maskArtery = results{1};
    maskVein = results{2};
end

%% Clearing 3)
%% 3) 1) Morphological Operations
% Preallocate a cell array to store the results
results = cell(2, 1);

parfor i = 1:2
    if i == 1
        % Process artery mask
        results{i} = clearMasks(ToolBox, maskArtery, 'artery', cmapArtery);
    else
        % Process vein mask
        results{i} = clearMasks(ToolBox, maskVein, 'vein', cmapVein);
    end
end

% Assign results to their respective variables
maskArteryClean = results{1};
maskVeinClean = results{2};

%% RGB

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArteryClean;
RGBM0(:, :, 2) = rescale(M0_ff_img);
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVeinClean;

saveImage(RGBM0, ToolBox, 'all_3_2_RGB.png', isStep = true)

%% Segmentation Scores Calculation

segmentationScores(maskArtery, maskVein);

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

saveImage(maskCRA, ToolBox, 'maskCRA.png')
saveImage(maskCRV, ToolBox, 'maskCRV.png')

clear f_AVG_std f_AVG_mean
%% Save Images

saveImage(maskArtery, ToolBox,'artery_3_3_Final.png', isStep = true, cmap = cmapArtery)
saveImage(maskVein, ToolBox, 'vein_3_3_Final.png', isStep = true, cmap = cmapVein)

saveImage(maskArtery, ToolBox, 'maskArtery.png')
saveImage(maskVein, ToolBox, 'maskVein.png')
saveImage(maskVessel, ToolBox, 'maskVessel.png')
saveImage(maskNeighbors, ToolBox, 'maskNeighbors.png')
saveImage(maskBackground, ToolBox, 'maskBackground.png')

if veinsAnalysis
    neighborsMaskSeg(:, :, 1) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors) + maskVein;
    saveImage(neighborsMaskSeg, ToolBox, 'NeighborsOnVessels.png')
else
    neighborsMaskSeg(:, :, 1) = M0_ff_img .* ~(maskArtery  | maskNeighbors) + maskArtery;
    neighborsMaskSeg(:, :, 2) = M0_ff_img .* ~(maskArtery  | maskNeighbors) + maskNeighbors;
    neighborsMaskSeg(:, :, 3) = M0_ff_img .* ~(maskArtery | maskNeighbors);
    saveImage(neighborsMaskSeg, ToolBox, 'NeighborsOnVessels.png')
end

%% Mask Section & Force Barycenter

L = (numY + numX) / 2;
r1 = PW_params.velocitySmallRadiusRatio * L;
r2 = PW_params.velocityBigRadiusRatio * L;

if ~isempty(PW_params.forcebarycenter)
    y_CRA = PW_params.forcebarycenter(1);
    x_CRA = PW_params.forcebarycenter(2);
    maskCircle = diskMask(numX, numY, PW_params.masks_crop_radius, "options", [x_CRA, y_CRA]);
end

xy_barycenter = [x_CRA, y_CRA];
maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMapArtery', maskArtery);
createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMap', maskArtery, maskVein);

%% Create Segmentation Map

if veinsAnalysis
    segmentationMap(:, :, 1) = M0_ff_img .* ~(maskArtery | maskVein) + maskArtery;
    segmentationMap(:, :, 2) = M0_ff_img .* ~(maskArtery | maskVein);
    segmentationMap(:, :, 3) = M0_ff_img .* ~(maskArtery | maskVein) + maskVein;
    saveImage(segmentationMap, ToolBox, 'Segmentation.png')
else
    segmentationMapArtery(:, :, 1) = M0_ff_img .* ~maskArtery + maskArtery;
    segmentationMapArtery(:, :, 2) = M0_ff_img .* ~maskArtery;
    segmentationMapArtery(:, :, 3) = M0_ff_img .* ~maskArtery;
    saveImage(segmentationMapArtery, ToolBox, 'Segmentation.png')
end

%% new masks

labeled = bwlabel(and(maskArtery, not(maskCircle)));
saveImage(rescale(labeled), ToolBox, 'maskLabeled.png')
saveImage(bwskel(maskArtery), ToolBox, 'skeletonArtery.png')
saveImage(bwskel(maskVein), ToolBox, 'skeletonVein.png')

disp('Done!')

end
