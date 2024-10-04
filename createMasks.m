function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(flatfieldM0, videoM0, videoM1M0, path, ToolBox)

    %% loading parameters and compute useful variables
    PW_params = Parameters_json(path);
    mkdir(ToolBox.PW_path_png, 'mask')

    [numX, numY, numFrames] = size(flatfieldM0);

    meanIm = squeeze(mean(flatfieldM0, 3));
    meanM0 = squeeze(mean(videoM0, 3));
    meanM1M0 = squeeze(mean(videoM1M0, 3));
    figure(666), imagesc(meanIm), axis image, colormap gray;
    
    videoM0Centered = videoM0 - meanM0;

    % compute vesselness response
    localContrastedIm = adapthisteq(rescale(meanIm), 'NumTiles', PW_params.arteryMask_vesselnessContrastNumSlides');
    vesselnessIm = vesselness_filter(localContrastedIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
    figure(667), imagesc(vesselnessIm), axis image, colormap gray;

    %% Create CRA and CRV Mask

    stdM1M0 = std2(meanM1M0);
    maskCRA = meanM1M0 > (PW_params.CRACRV_Threshold * stdM1M0);
    maskCRV = meanM1M0 < (-PW_params.CRACRV_Threshold * stdM1M0);

    %%  Compute first correlation to find arteries

    % compute pulse in 3 dimentions for correlation in all vessels
    vesselPulse = squeeze(mean(videoM0 .* (vesselnessIm > 0), [1 2]));
    vesselPulseInit = vesselPulse - mean(vesselPulse, "all");
    vesselPulseInit3d = zeros(numX, numY, numFrames);

    parfor xx = 1:numX

        for yy = 1:numY
            vesselPulseInit3d(xx, yy, :) = vesselPulseInit;
        end

    end

    % compute local-to-average pulse wave zero-lag correlation
    R_VesselPulse = squeeze(mean((videoM0Centered .* vesselPulseInit3d), 3)) .* (vesselnessIm > 0);

    % Create first artery mask
    firstMaskArtery = (R_VesselPulse > PW_params.arteryMask_CorrelationMatrixThreshold * mean2(R_VesselPulse(R_VesselPulse > 0)));

    if PW_params.masks_cleaningCoroid
        firstMaskArtery = firstMaskArtery & bwareafilt(firstMaskArtery , 1, 4);
    end

    firstMaskArtery = bwareaopen(firstMaskArtery, PW_params.masks_minSize);

    figure(10), imagesc(vesselnessIm);
    title('Vesselness map');
    colorbar('southoutside');
    figure(11), imshow(firstMaskArtery)

    clear vesselPulse vesselPulseInit vesselPulseInit3d R_VesselPulse;
    %% Compute correlation to segment veins and arteries

    % compute pulse in 3 dimentions for correlation in main arteries
    pulse = squeeze(mean(videoM0 .* firstMaskArtery, [1 2]));
    pulseInit = pulse - mean(pulse, "all");
    pulseInit3d = zeros(numX, numY, numFrames);

    parfor xx = 1:numX

        for yy = 1:numY
            pulseInit3d(xx, yy, :) = pulseInit;
        end

    end

    % compute local-to-average pulse wave zero-lag correlation
    correlationMatrix = squeeze(mean((videoM0Centered .* pulseInit3d), 3));
    
    % Create correlation matrix to segment vein and arteries
    correlationMatrixArtery = correlationMatrix ./ max(correlationMatrix, [], 'all');
    correlationMatrixVein = correlationMatrix ./ min(correlationMatrix, [], 'all');
    correlationMatrixVein(correlationMatrixVein <- 1) = -1;

    figure(20), imagesc(correlationMatrixArtery);
    title('Artery correlation map');
    colorbar('southoutside');
    imwrite(rescale(correlationMatrixArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrix_artery.png')))

    figure(21), imagesc(correlationMatrixVein);
    title('Vein correlation map');
    colorbar('southoutside');
    imwrite(rescale(correlationMatrixVein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrix_vein.png')))

    clear pulse pulseInit pulseInit3d;

    %% Circles Sectioning

    blurred_mask = imgaussfilt(double(meanM0 .* meanM1M0).*correlationMatrixArtery, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0);
    [ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));

    [x, y] = meshgrid(1:numY, 1:numX);
    cercleMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

    radiusTreshold = PW_params.masks_radius_treshold;
    cercleTresholdMask = [];
    minDistToEdge = min([ToolBox.x_barycentre, ToolBox.y_barycentre, numX - ToolBox.x_barycentre, numY - ToolBox.y_barycentre]);

    if radiusTreshold == 0

        cercleTresholdMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= minDistToEdge;
    elseif radiusTreshold > 0
        cercleTresholdMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= (radiusTreshold * minDistToEdge);
    else % if radius_treshold==-1
    end

    %% Region growing to get a clear segmentation

    % cleaning correlation matrix
    correlationMatrixArtery(vesselnessIm == 0) = min(correlationMatrixArtery, [], 'all');
    correlationMatrixVein(vesselnessIm == 0) = min(correlationMatrixVein, [], 'all');

    % Create vesselness response weighted by correlation
    vesselnessArtery = vesselnessIm .* rescale(correlationMatrixArtery);
    vesselnessVein = vesselnessIm .* rescale(correlationMatrixVein);

    % Create seed artery mask and condition artery mask for region growing
    correlationMatrixArteryCondition = correlationMatrixArtery .* rescale(vesselnessIm);

    levelArtery = graythresh(correlationMatrixArteryCondition);
    seedsArtery = imbinarize(correlationMatrixArteryCondition, levelArtery ./ PW_params.RG_ArterySeedsThreshold);
    seedsArtery = bwareaopen(seedsArtery, PW_params.masks_minSize);

    conditionArtery = imbinarize(correlationMatrixArteryCondition, levelArtery ./ PW_params.RG_ArteryConditionThreshold);

    % Create seed vein mask and contidion vein mask for region growing
    correlationMatrix_veinCondition = correlationMatrixVein .* rescale(vesselnessIm);

    levelVein = graythresh(correlationMatrix_veinCondition);
    seedsVein = imbinarize(correlationMatrix_veinCondition, levelVein) | (vesselnessVein > 0.5 * max(vesselnessVein, [], 'all'));
    seedsVein = bwareaopen(seedsVein, PW_params.masks_minSize);

    condition_vein = vesselnessVein > PW_params.RG_veinConditionThreshold * mean2(vesselnessVein(seedsVein));

    % Cleaning condition
    conditionArtery = conditionArtery & bwareafilt(conditionArtery | condition_vein | cercleMask, 1, 4);
    conditionArtery = bwareaopen(conditionArtery, PW_params.masks_minSize);
    condition_vein = condition_vein & ~conditionArtery;
    condition_vein = bwareaopen(condition_vein, PW_params.masks_minSize);

    % Cleaning seeds
    seedsArtery = seedsArtery & conditionArtery;
    seedsVein = seedsVein & condition_vein;

    % Compute region growing to segment vein and arteries
    %[maskArtery,RG_video_artery] = region_growing_for_vessel(vesselness_artery, seeds_artery, condition_artery,path);
    %[maskVein,RG_video_vein]  = region_growing_for_vessel(vesselness_vein, seeds_vein, condition_vein, path);
    [maskVessel, rgVideoVessel] = region_growing_for_vessel(vesselnessIm, seedsArtery | seedsVein, condition_vein | conditionArtery, path);

    maskVessel = bwareafilt(maskVessel, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);
    maskArtery = maskVessel .* correlationMatrixArtery;
    maskArtery = maskArtery > PW_params.arteryMask_ArteryCorrThreshold;
    maskVein = maskVessel .* correlationMatrixVein;
    maskVein = maskVein > 0;

    %maskVein = maskVessel & ~maskArtery;

    if PW_params.masks_showIntermediateFigures
        figure(30), imagesc(vesselnessArtery);
        title('Artery Vesselness map');
        colorbar('southoutside');
        figure(31), imagesc(vesselnessVein);
        title('Vein Vesselness map');
        colorbar('southoutside');
        figure(32), imagesc(uint8(cat(3, uint8(meanIm) + uint8(seedsArtery) * 255, uint8(meanIm), uint8(meanIm) + uint8(seedsVein) * 255)));
        title('Artery/Vein initialisation for region growing');
        figure(33), imagesc(uint8(cat(3, uint8(meanIm) + uint8(conditionArtery) * 255, uint8(meanIm), uint8(meanIm) + uint8(condition_vein) * 255)));
        title('Artery/Vein condition for region growing');
        figure(34), imshow(maskVessel);
    end

    clear vesselnessVein vesselnessArtery floor_vein floor_artery levelVein levelArtery seedsVein seedsArtery ;
    %% Cleaning coroid from masks

    if ~isempty(cercleTresholdMask)
        maskArtery = maskArtery & cercleTresholdMask;
    end

    % if PW_params.masks_cleaningCoroid
    %
    %     maskArtery = maskArtery & bwareafilt(maskArtery | maskVein | cercle_mask,1,4);
    %     maskVein = maskVein & bwareafilt(maskArtery | maskVein | cercle_mask,1,4);
    %
    %     if PW_params.masks_showIntermediateFigures
    %     figure(40), imagesc(uint8( cat( 3, uint8(meanIm)+ uint8(maskArtery)*255, uint8(meanIm) + uint8(cercle_mask) , uint8(meanIm) + uint8(maskVein)*255 )));
    %     title('Artery/Vein cleaned segmentation');
    %     end
    %
    %     clear cercle_mask;
    % end

    maskArtery = bwareaopen(maskArtery, PW_params.masks_minSize);
    maskArtery = imdilate(maskArtery, strel('disk', 3));
    maskArtery = imclose(maskArtery, strel('disk', 5));

    maskArtery = bwareafilt(maskArtery, PW_params.arteryMask_magicwand_nb_of_area_artery, 4);

    maskVein = bwareaopen(maskVein, PW_params.masks_minSize);
    maskVein = imdilate(maskVein, strel('disk', 3));
    maskVein = imclose(maskVein, strel('disk', 5)) & ~maskArtery;

    maskVessel = maskArtery | maskVein;

    if PW_params.masks_showIntermediateFigures
        figure(40), imagesc(uint8(cat(3, uint8(meanIm) + uint8(maskArtery) * 255, uint8(meanIm), uint8(meanIm) + uint8(maskVein) * 255)));
        title('Artery/Vein region growing segmentation');
    end

    %import manual mask if provided
    try

        if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'maskVessel_New.png')) % import manually tuned mask if needed
            maskVessel = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'maskVessel_New.png')), 3)) > 0;
        end

        if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'maskArtery_New.png')) % import manually tuned mask if needed
            maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'maskArtery_New.png')), 3)) > 0;
        end

        if isfile(fullfile(ToolBox.PW_path_main, 'mask', 'maskVein_New.png')) % import manually tuned mask if needed
            maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'maskVein_New.png')), 3)) > 0;
        end

    catch
        warning("There was an error with the manual import of the masks, check if the dimensions of the mask match the ones of the video.\n Normal procedure will follow.")        
    end

    %% Create Background Mask

    maskBackground = not(maskVessel);

    %% Create Mask Section

    radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (numY + numX) / 2;
    radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (numY + numX) / 2;

    circleMask1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius1;
    circleMask2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius2;

    maskSection = xor(circleMask1, circleMask2);

    %% Create Colormap Artery/Vein

    meanIm = mat2gray(meanIm);
    [hueArtery, satArtery, val] = createHSVmap(meanIm, maskArtery - maskArtery .* maskSection, 0, 0);
    [hueVein, satVein, ~] = createHSVmap(meanIm, maskVein - maskVein .* maskSection - maskVein .* maskArtery, 0.7, 0.7);
    [hueSectionA, satSectionA, ~] = createHSVmap(meanIm, maskSection .* maskArtery, 0.15, 0.15);
    [hueSectionV, satSectionV, ~] = createHSVmap(meanIm, maskSection .* maskVein, 0.5, 0.5);
    satSectionArtery = satSectionA;
    satSectionVein = satSectionV; %+~maskSectionArtery.*(~maskArtery);
    val = val .* (~maskSection) + val .* maskSection + maskSection .* (~(maskArtery + maskVein));
    vesselImageRGB = hsv2rgb(hueArtery + hueVein + hueSectionA + hueSectionV, satArtery + satVein + satSectionArtery + satSectionVein, val);

    figure(101)
    imshow(vesselImageRGB)

    %% Create Colormap Artery Only
    [hueArtery, satArtery, val] = createHSVmap(meanIm, maskArtery - maskArtery .* maskSection, 0, 0);
    [hueSectionA, satSectionA, ~] = createHSVmap(meanIm, maskSection .* maskArtery, 0.15, 0.15);
    satSectionArtery = satSectionA;
    val = val .* (~maskSection) + val .* maskSection + maskSection .* (~maskArtery);
    VesselImageRGB_Artery = hsv2rgb(hueArtery + hueSectionA, satArtery + satSectionArtery, val);

    figure(102)
    imshow(VesselImageRGB_Artery)

    %% Create Segmentation Map

    segmentationMap = zeros(numX, numY, 3);
    segmentationMap(:, :, 1) = meanIm - (maskArtery + maskVein) .* meanIm + maskArtery;
    segmentationMap(:, :, 2) = meanIm - (maskArtery + maskVein) .* meanIm;
    segmentationMap(:, :, 3) = meanIm - (maskArtery + maskVein) .* meanIm + maskVein;
    imwrite(segmentationMap, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arteryVeinSegmentation.png')), 'png');

    %% Saving masks as PNG
    foldername = ToolBox.main_foldername;

    imwrite(mat2gray(double(vesselnessIm)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselness.png')), 'png');
    imwrite(mat2gray(single(maskArtery)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskArtery_New.png')), 'png');
    imwrite(mat2gray(single(maskVein)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVein_New.png')), 'png');
    imwrite(mat2gray(single(maskVessel)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVessel_New.png')), 'png');
    imwrite(mat2gray(single(maskBackground)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskBackground_New.png')), 'png');
    %vesselMap = uint8( cat( 3, uint8(meanIm)+ uint8(maskArtery)*255, uint8(meanIm) , uint8(meanIm) + uint8(maskVein)*255 ));
    imwrite(vesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMap.png')), 'png');
    imwrite(VesselImageRGB_Artery, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselMapArtery.png')), 'png');
    imwrite(mat2gray(maskCRA), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRA.png')), 'png');
    imwrite(mat2gray(maskCRV), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRV.png')), 'png');

    %% new masks

    labeled = bwlabel(and(maskArtery, not(circleMask1)));
    imwrite(rescale(labeled), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskLabeled.png')), 'png');
    imwrite(bwskel(maskArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSkeleton.png')), 'png');

    %% Saving AVI

    % w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_RG_video_artery.avi')));
    % tmp = mat2gray(RG_video_artery);
    % open(w)
    % for j = 1:size(RG_video_artery,3)
    %     writeVideo(w,tmp(:,:,j)) ;
    % end
    % close(w);
    %
    % w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_RG_video_vein.avi')));
    % tmp = mat2gray(RG_video_vein);
    % open(w)
    % for j = 1:size(RG_video_vein,3)
    %     writeVideo(w,tmp(:,:,j)) ;
    % end
    % close(w);

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_RG_video_vessel.avi')));
    tmp = mat2gray(rgVideoVessel);
    open(w)

    for j = 1:size(rgVideoVessel, 3)
        writeVideo(w, tmp(:, :, j));
    end

    close(w);

end
