function [maskArtery, maskVein, maskSection, maskNeighbors, xy_barycenter] = createMasks(M0_ff_video, f_AVG_mean)

    TB = getGlobalToolBox;
    params = TB.getParams;

    if ~isfolder(fullfile(TB.path_png, 'mask'))
        mkdir(TB.path_png, 'mask')
        mkdir(TB.path_eps, 'mask')
        mkdir(fullfile(TB.path_png, 'mask'), 'steps')
        mkdir(fullfile(TB.path_eps, 'mask'), 'steps')
    end

    folder_steps = fullfile('mask', 'steps');

    %% 0) Initialisation

    % 0) 1) Parameters Initialisation

    [numX, numY, numFrames] = size(M0_ff_video);

    diaphragmRadius = params.json.Mask.DiaphragmRadius;
    forceBarycenter = params.json.Mask.ForceBarycenter;
    blur = params.json.Mask.Blur;
    cropChoroidRadius = params.json.Mask.CropChoroidRadius;

    % Parameters for arteries
    vesselParams.threshold = params.json.Mask.VascularThreshold;
    vesselParams.classes = params.json.Mask.VascularClasses;

    diasysAnalysis = params.json.Mask.DiaSysAnalysis;

    % Parameters for arteries
    arteryParams.threshold = params.json.Mask.ArterialThreshold;
    arteryParams.classes = params.json.Mask.ArterialClasses;

    % Parameters for veins
    veinParams.threshold = params.json.Mask.VenousThreshold;
    veinParams.classes = params.json.Mask.VenousClasses;

    minPixelSize = params.json.Mask.MinPixelSize;
    CRACRV_Threshold = params.json.Mask.CRACRVThreshold;
    forceVesselWidth = params.json.Mask.ForceVesselWidth;

    bgWidth = params.json.Velocity.LocalBackgroundWidth;
    r1 = params.json.SizeOfField.SmallRadiusRatio;
    r2 = params.json.SizeOfField.BigRadiusRatio;

    % 0) 2) Test the input for usual cases of wrong aquisition data

    if max(f_AVG_mean, [], 'all') <= 0
        figure(3); imshow(rescale(f_AVG_mean));
        error("Measurement error from the Moment 1 input. Check the input or re-do the measurement.")
    end

    %% 1) First Masks and Correlation

    maskDiaphragm = diskMask(numX, numY, diaphragmRadius);

    M0_ff_img = squeeze(mean(M0_ff_video, 3));
    M0_ff_video_centered = M0_ff_video - mean(M0_ff_video, [1 2]);
    saveImage(M0_ff_img, TB, 'all_10_M0.png', isStep = true)

    if ~isfile(fullfile(TB.path_gif, sprintf("%s_M0.gif", TB.folder_name)))
        writeGifOnDisc(rescale(M0_ff_video), "M0")
    end

    saveImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, TB, 'all_11_maskDiaphragm.png', isStep = true)

    % 1) 1) Compute vesselness response

    [maskVesselnessFrangi] = frangiVesselness(M0_ff_img, 'all_12', TB);
    [maskVesselnessGabor, M0_Gabor] = gaborVesselness(M0_ff_img, 'all_13', TB);

    maskVesselness = (maskVesselnessFrangi | maskVesselnessGabor) & maskDiaphragm;

    % 1) 2) Compute the barycenters and the circle mask

    if ~isempty(forceBarycenter)
        y_CRA = forceBarycenter(1);
        x_CRA = forceBarycenter(2);
        x_CRV = numX / 2;
        y_CRV = numY / 2;
    else
        saveImage(f_AVG_mean, TB, 'all_13_fAVG.png', isStep = true)

        if blur ~= 0
            averaged_fAVG = imgaussfilt(f_AVG_mean, blur, 'Padding', 0) .* maskDiaphragm;
        else
            averaged_fAVG = f_AVG_mean;
        end

        [y_CRA, x_CRA] = find(averaged_fAVG == max(averaged_fAVG, [], 'all'));
        [y_CRV, x_CRV] = find(averaged_fAVG == min(averaged_fAVG, [], 'all'));
    end

    maskCircle = diskMask(numX, numY, cropChoroidRadius, 'center', [x_CRA / numX, y_CRA / numY]);
    maskCircle = maskCircle | diskMask(numX, numY, cropChoroidRadius, 'center', [x_CRV / numX, y_CRV / numY]);

    maskVesselnessClean = removeDisconnected(maskVesselness, maskVesselness, maskCircle, 'all_14_VesselMask', TB);

    %  1) 3) Compute first correlation

    cVascular = [0 0 0];
    % compute signal in 3 dimentions for correlation in all vessels
    vascularSignal = sum(M0_ff_video .* maskVesselnessClean, [1 2]);
    vascularSignal = vascularSignal ./ nnz(maskVesselnessClean);

    t = linspace(0, numFrames * TB.stride / TB.fs / 1000, numFrames);
    tLabel = 'Time(s)';
    yLabel = 'Power Doppler (a.u.)';

    graphSignal('all_15_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, Title = 'Vascular Signal', xlabel = tLabel, ylabel = yLabel);

    % compute local-to-average signal wave zero-lag correlation
    vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);
    R_VascularSignal = mean(M0_ff_video_centered .* vascularSignal_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularSignal_centered, [], 3));
    saveImage(R_VascularSignal, TB, 'all_15_Correlation.png', isStep = true)

    mR_vascular = sum(R_VascularSignal .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));

    RGBcorr = labDuoImage(M0_ff_img, R_VascularSignal - mR_vascular);
    saveImage(RGBcorr, TB, 'all_15_Correlation_rgb.png', isStep = true)

    % 1) 4) Segment Vessels

    cArtery = [255 22 18] / 255;
    cVein = [18 23 255] / 255;

    cmapArtery = [0 0 0; cArtery];
    cmapVein = [0 0 0; cVein];
    cmapVessels = [cVein; cArtery];

    if vesselParams.threshold >= -1 && vesselParams.threshold <= 1
        % IF Manual Thresholds have been set between -1 and 1 then they are used

        maskArtery = (R_VascularSignal > vesselParams.threshold) .* maskVesselnessClean;
        maskVein = (R_VascularSignal < vesselParams.threshold) .* maskVesselnessClean;
        graphThreshHistogram(R_VascularSignal, vesselParams.threshold, maskVesselnessClean, cmapVessels, 'all_16')

    else
        % ELSE automatic Otsu segmentation is performed
        % Number of classes for Vessels: 4
        % 1 & 2 = Veins & CoroidalVessels, 3 = CoroidalVessel, 4 = Arteries
        [maskArtery, maskVein] = autoOtsuThresholding(R_VascularSignal, maskVesselnessClean, vesselParams.classes, 'all_16');
    end

    saveImage(maskArtery, TB, 'artery_17_FirstMask.png', isStep = true, cmap = cmapArtery)
    saveImage(maskVein, TB, 'vein_17_FirstMask.png', isStep = true, cmap = cmapVein)

    % Remove small blobs
    maskArtery = bwareaopen(maskArtery, minPixelSize);
    maskVein = bwareaopen(maskVein, minPixelSize);

    saveImage(maskArtery, TB, 'artery_18_FirstMaskClean.png', isStep = true, cmap = cmapArtery)
    saveImage(maskVein, TB, 'vein_18_FirstMaskClean.png', isStep = true, cmap = cmapVein)

    RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArtery;
    RGBM0(:, :, 2) = rescale(M0_ff_img);
    RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVein;
    saveImage(RGBM0, TB, 'all_19_RGB.png', isStep = true)

    %% 2)  Improvements of the first mask

    if params.json.Mask.ImproveMask

        % 2) 0) Computation of the M0 in Diastole and in Systole

        [M0_Systole_img, M0_Diastole_img, M0_Systole_video, ~] = compute_diasys(M0_ff_video, maskArtery);
        saveImage(rescale(M0_Systole_img), TB, 'artery_20_systole_img.png', isStep = true)
        saveImage(rescale(M0_Diastole_img), TB, 'vein_20_diastole_img.png', isStep = true)

        % 2) 1) New Vesselness Mask

        Systole_Frangi = frangiVesselness(M0_Systole_img, 'artery_20', TB);
        Diastole_Frangi = frangiVesselness(M0_Diastole_img, 'vein_20', TB);
        Systole_Gabor = gaborVesselness(M0_Systole_img, 'artery_20', TB);
        Diastole_Gabor = gaborVesselness(M0_Diastole_img, 'vein_20', TB);
        maskVesselness = (Systole_Frangi | Diastole_Frangi | Systole_Gabor | Diastole_Gabor) & maskDiaphragm;
        maskVesselnessClean = removeDisconnected(maskVesselness, maskVesselness, maskCircle, 'all_20_VesselMask', TB);

        % 2) 2) Diastole-Systole Image

        diasysArtery = M0_Systole_img - M0_Diastole_img;
        mDiasys = sum(diasysArtery .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));
        diasysVein = mDiasys - diasysArtery;
        saveImage(diasysArtery, TB, 'artery_21_diasys_img.png', isStep = true)
        saveImage(diasysVein, TB, 'vein_21_diasys_img.png', isStep = true)

        RGBdiasys = labDuoImage(rescale(M0_Gabor), (diasysArtery - mDiasys));
        saveImage(RGBdiasys, TB, 'vessel_40_diasys_rgb.png', isStep = true)
        saveImage(RGBdiasys, TB, 'DiaSysRGB.png')

        if diasysAnalysis % Systole/Diastole Analysis

            % 2) 3) Diastole-Systole based Segmentation
            maskArtery = processDiaSysSignal(diasysArtery, maskVesselnessClean, arteryParams, cmapArtery, 'artery_23');
            [~, maskVein] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23');

        else % Second Correlation Analysis

            % 2) 3) Artery-Vein correlation based Segmentation
            maskArtery = processVascularSignal(M0_Systole_video, maskArtery, maskVesselnessClean, arteryParams, cmapArtery, 'artery_23', TB);
            %         maskVein = processVascularSignal(M0_Diastole_video, maskVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23', TB);
            [~, maskVein] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23');

        end

        %% 3) Mask Clearing

        % 3) 0) Morphological Operations
        results = cell(2, 1);

        parfor i = 1:2

            if i == 1
                % Process artery mask
                results{i} = clearMasks(maskArtery, 'artery_30', cmapArtery, TB);
            else
                % Process vein mask
                results{i} = clearMasks(maskVein, 'vein_30', cmapVein, TB);
            end

        end

        maskArtery = results{1};
        maskVein = results{2};

        % 3) 1) Final Blob removal
        maskVessel = maskArtery | maskVein;
        maskArtery = removeDisconnected(maskArtery, maskVessel, maskCircle, 'artery_31_VesselMask', TB);
        maskVein = removeDisconnected(maskVein, maskVessel, maskCircle, 'vein_31_VesselMask', TB);

        % 3) 1 prime) HoloNet intervention
        % holonet_vessels = getHolonetprediction(M0_ff_img);
        % maskVessel = maskVessel & holonet_vessels;
        % maskArtery = maskArtery & holonet_vessels;
        % maskVein = maskVein & holonet_vessels;

        % 3) 2) Force Create Masks in case they exist

        if isfile(fullfile(TB.path_main, 'mask', 'forceMaskArtery.png'))
            maskArtery = mat2gray(mean(imread(fullfile(TB.path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;

            if size(maskArtery, 1) ~= maskCircle
                maskArtery = imresize(maskArtery, [numX, numY], "nearest");
            end

        end

        if isfile(fullfile(TB.path_main, 'mask', 'forceMaskVein.png'))
            maskVein = mat2gray(mean(imread(fullfile(TB.path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;

            if size(maskVein, 1) ~= maskCircle
                maskVein = imresize(maskVein, [numX, numY], "nearest");
            end

        end

        % 3) 3) Segmentation Scores Calculation

        segmentationScores(maskArtery, maskVein);

        % 3) 4) Segmention force width

        if forceVesselWidth > 0
            dilationSE = strel('disk', forceVesselWidth);
            maskArtery = imdilate(bwskel(maskArtery), dilationSE);
            maskVein = imdilate(bwskel(maskVein), dilationSE);

            maskArtery = imdilate(bwskel(maskArtery), dilationSE);
            maskVein = imdilate(bwskel(maskVein), dilationSE);
        end

    else

        % 3) 0) Morphological Operations
        results = cell(2, 1);

        parfor i = 1:2

            if i == 1
                % Process artery mask
                results{i} = clearMasks(maskArtery, 'artery_30', cmapArtery, TB);
            else
                % Process vein mask
                results{i} = clearMasks(maskVein, 'vein_30', cmapVein, TB);
            end

        end

        maskArtery = results{1};
        maskVein = results{2};

        maskVessel = maskArtery | maskVein;
        maskArtery = removeDisconnected(maskArtery, maskVessel, maskCircle, 'artery_31_VesselMask', TB);
        maskVein = removeDisconnected(maskVein, maskVessel, maskCircle, 'vein_31_VesselMask', TB);

    end

    % 3) 5) Create Vessel and Background Mask

    maskVessel = maskArtery | maskVein;
    maskBackground = not(maskVessel);

    %% 4) FINAL FIGURES

    % 4) 1) RGB Figures
    cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
    cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
    cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

    M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
    M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
    M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

    M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
    saveImage(M0_RGB, TB, 'vessel_40_RGB.png', isStep = true)
    saveImage(M0_RGB, TB, 'RGB_img.png')

    % 4) 2) Neighbours Mask

    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', bgWidth)) - (maskArtery | maskVein);

    cmapNeighbors = cmapLAB(256, [0 1 0], 0, [1 1 1], 1);

    M0_Neighbors = setcmap(M0_ff_img, maskNeighbors, cmapNeighbors);

    neighborsMaskSeg = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + ...
        M0_AV + M0_Neighbors + ...
        rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors);
    saveImage(neighborsMaskSeg, TB, 'neighbors_img.png')

    % 4) 3) CRA and CRV Masks

    f_AVG_std = std2(f_AVG_mean);
    maskCRA = f_AVG_mean > (CRACRV_Threshold * f_AVG_std);
    maskCRV = f_AVG_mean < (-CRACRV_Threshold * f_AVG_std);

    saveImage(maskCRA, TB, 'maskCRA.png')
    saveImage(maskCRV, TB, 'maskCRV.png')

    % 4) 4) Save all images

    saveImage(maskArtery, TB, 'maskArtery.png')
    saveImage(maskVein, TB, 'maskVein.png')
    saveImage(maskVessel, TB, 'maskVessel.png')
    saveImage(maskNeighbors, TB, 'maskNeighbors.png')
    saveImage(maskBackground, TB, 'maskBackground.png')
    saveImage(bwskel(maskArtery), TB, 'skeletonArtery.png')
    saveImage(bwskel(maskVein), TB, 'skeletonVein.png')

    % 4) 5) Mask Section & Force Barycenter

    xy_barycenter = [x_CRA, y_CRA];
    maskSection = createMaskSection(TB, M0_ff_img, r1, r2, xy_barycenter, 'vesselMapArtery', maskArtery, thin = 0.01);
    createMaskSection(TB, M0_ff_img, r1, r2, xy_barycenter, 'vesselMap', maskArtery, maskVein, thin = 0.01);

    close all
end
