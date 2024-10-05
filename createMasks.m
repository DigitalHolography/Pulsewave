function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(M0_disp_video, f_RMS_video, f_AVG_mean, path, ToolBox)

    %% loading parameters and compute useful variables
    PW_params = Parameters_json(path);
    mkdir(ToolBox.PW_path_png, 'mask')

    [numX, numY, ~] = size(M0_disp_video);

    M0_disp_img = squeeze(mean(M0_disp_video, 3));
    figure(666), imagesc(M0_disp_img), axis image, colormap gray;

    f_RMS_mean = squeeze(mean(f_RMS_video, 3));
    f_RMS_centered_video = f_RMS_video - f_RMS_mean;

    % compute vesselness response
    M0_contrasted_img = adapthisteq(rescale(M0_disp_img), 'NumTiles', PW_params.arteryMask_vesselnessContrastNumSlides');
    vesselness_img = vesselness_filter(M0_contrasted_img, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
    figure(667), imagesc(vesselness_img), axis image, colormap gray;

    %% Create CRA and CRV Mask

    f_AVG_std = std2(f_AVG_mean);
    maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
    maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);

    %%  Compute first correlation to find arteries

    % compute pulse in 3 dimentions for correlation in all vessels
    vascularPulse = squeeze(mean(f_RMS_video .* (vesselness_img > 0), [1 2]));
    [~, ~, vascularPulse_video] = meshgrid(ones(1, numY), ones(1, numX), vascularPulse - mean(vascularPulse, "all"));

    % compute local-to-average pulse wave zero-lag correlation
    R_VascularPulse = squeeze(mean((f_RMS_centered_video .* vascularPulse_video), 3)) .* (vesselness_img > 0);

    % Create first artery mask
    firstMaskArtery = (R_VascularPulse > PW_params.arteryMask_CorrelationMatrixThreshold * mean2(R_VascularPulse(R_VascularPulse > 0)));

    if PW_params.masks_cleaningCoroid
        firstMaskArtery = firstMaskArtery & bwareafilt(firstMaskArtery , 1, 4);
    end

    firstMaskArtery = bwareaopen(firstMaskArtery, PW_params.masks_minSize);

    figure(10), imagesc(vesselness_img);
    title('Vesselness map');
    colorbar('southoutside');
    figure(11), imshow(firstMaskArtery)

    clear vascularPulse vascularPulseInit vascularPulseInit3d R_VascularPulse;
    %% Compute correlation to segment veins and arteries

    % compute pulse in 3 dimentions for correlation in main arteries
    arterialPulse = squeeze(mean(f_RMS_video .* firstMaskArtery, [1 2]));
    [~, ~, arterialPulse_video] = meshgrid(ones(1, numY), ones(1, numX), arterialPulse - mean(arterialPulse, "all"));


    % compute local-to-average pulse wave zero-lag correlation
    R_ArterialPulse = squeeze(mean((f_RMS_video .* arterialPulse_video), 3));
    
    % Create correlation matrix to segment vein and arteries
    R_Artery = R_ArterialPulse ./ max(R_ArterialPulse, [], 'all');
    R_Vein = R_ArterialPulse ./ min(R_ArterialPulse, [], 'all');
    R_Vein(R_Vein <- 1) = -1;

    figure(20), imagesc(R_Artery);
    title('Artery correlation map');
    colorbar('southoutside');
    imwrite(rescale(R_Artery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrix_artery.png')))

    figure(21), imagesc(R_Vein);
    title('Vein correlation map');
    colorbar('southoutside');
    imwrite(rescale(R_Vein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrix_vein.png')))

    clear arterialPulse arterialPulseInit arterialPulse_video;

    %% Circles Sectioning

    blurred_mask = imgaussfilt(double(f_RMS_mean .* f_AVG_mean).*R_Artery, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0);
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
    R_Artery(vesselness_img == 0) = min(R_Artery, [], 'all');
    R_Vein(vesselness_img == 0) = min(R_Vein, [], 'all');

    % Create vesselness response weighted by correlation
    vesselnessArtery = vesselness_img .* rescale(R_Artery);
    vesselnessVein = vesselness_img .* rescale(R_Vein);

    % Create seed artery mask and condition artery mask for region growing
    correlationMatrixArteryCondition = R_Artery .* rescale(vesselness_img);

    levelArtery = graythresh(correlationMatrixArteryCondition);
    seedsArtery = imbinarize(correlationMatrixArteryCondition, levelArtery ./ PW_params.RG_ArterySeedsThreshold);
    seedsArtery = bwareaopen(seedsArtery, PW_params.masks_minSize);

    conditionArtery = imbinarize(correlationMatrixArteryCondition, levelArtery ./ PW_params.RG_ArteryConditionThreshold);

    % Create seed vein mask and contidion vein mask for region growing
    R_VeinCondition = R_Vein .* rescale(vesselness_img);

    levelVein = graythresh(R_VeinCondition);
    seedsVein = imbinarize(R_VeinCondition, levelVein) | (vesselnessVein > 0.5 * max(vesselnessVein, [], 'all'));
    seedsVein = bwareaopen(seedsVein, PW_params.masks_minSize);

    conditionVein = vesselnessVein > PW_params.RG_veinConditionThreshold * mean2(vesselnessVein(seedsVein));

    % Cleaning condition
    conditionArtery = conditionArtery & bwareafilt(conditionArtery | conditionVein | cercleMask, 1, 4);
    conditionArtery = bwareaopen(conditionArtery, PW_params.masks_minSize);
    conditionVein = conditionVein & ~conditionArtery;
    conditionVein = bwareaopen(conditionVein, PW_params.masks_minSize);

    % Cleaning seeds
    seedsArtery = seedsArtery & conditionArtery;
    seedsVein = seedsVein & conditionVein;

    % Compute region growing to segment vein and arteries
    %[maskArtery,RGVideoArtery] = region_growing_for_vessel(vesselness_artery, seeds_artery, condition_artery,path);
    %[maskVein,RGVideoVein]  = region_growing_for_vessel(vesselness_vein, seeds_vein, condition_vein, path);
    [maskVessel, rgVideoVessel] = region_growing_for_vessel(vesselness_img, seedsArtery | seedsVein, conditionVein | conditionArtery, path);

    maskVessel = bwareafilt(maskVessel, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);
    maskArtery = maskVessel .* R_Artery;
    maskArtery = maskArtery > PW_params.arteryMask_ArteryCorrThreshold;
    maskVein = maskVessel .* R_Vein;
    maskVein = maskVein > 0;

    %maskVein = maskVessel & ~maskArtery;

    if PW_params.masks_showIntermediateFigures
        figure(30), imagesc(vesselnessArtery);
        title('Artery Vesselness map');
        colorbar('southoutside');
        figure(31), imagesc(vesselnessVein);
        title('Vein Vesselness map');
        colorbar('southoutside');
        figure(32), imagesc(uint8(cat(3, uint8(M0_disp_img) + uint8(seedsArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(seedsVein) * 255)));
        title('Artery/Vein initialisation for region growing');
        figure(33), imagesc(uint8(cat(3, uint8(M0_disp_img) + uint8(conditionArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(conditionVein) * 255)));
        title('Artery/Vein condition for region growing');
        figure(34), imshow(maskVessel);
    end

    clear vesselnessVein vesselnessArtery levelVein levelArtery seedsVein seedsArtery ;
    %% Cleaning coroid from masks

    if ~isempty(cercleTresholdMask)
        maskArtery = maskArtery & cercleTresholdMask;
    end

    if PW_params.masks_cleaningCoroid

        maskArtery = maskArtery & bwareafilt(maskArtery | maskVein | cercle_mask,1,4);
        maskVein = maskVein & bwareafilt(maskArtery | maskVein | cercle_mask,1,4);

        if PW_params.masks_showIntermediateFigures
        figure(40), imagesc(uint8( cat( 3, uint8(M0_disp_img)+ uint8(maskArtery)*255, uint8(M0_disp_img) + uint8(cercle_mask) , uint8(M0_disp_img) + uint8(maskVein)*255 )));
        title('Artery/Vein cleaned segmentation');
        end

        clear cercle_mask;
    end

    maskArtery = bwareaopen(maskArtery, PW_params.masks_minSize);
    maskArtery = imdilate(maskArtery, strel('disk', 3));
    maskArtery = imclose(maskArtery, strel('disk', 5));

    maskArtery = bwareafilt(maskArtery, PW_params.arteryMask_magicwand_nb_of_area_artery, 4);

    maskVein = bwareaopen(maskVein, PW_params.masks_minSize);
    maskVein = imdilate(maskVein, strel('disk', 3));
    maskVein = imclose(maskVein, strel('disk', 5)) & ~maskArtery;

    maskVessel = maskArtery | maskVein;

    if PW_params.masks_showIntermediateFigures
        figure(40), imagesc(uint8(cat(3, uint8(M0_disp_img) + uint8(maskArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(maskVein) * 255)));
        title('Artery/Vein region growing segmentation');
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

    M0_disp_img = mat2gray(M0_disp_img);
    [hueArtery, satArtery, val] = createHSVmap(M0_disp_img, maskArtery - maskArtery .* maskSection, 0, 0);
    [hueVein, satVein, ~] = createHSVmap(M0_disp_img, maskVein - maskVein .* maskSection - maskVein .* maskArtery, 0.7, 0.7);
    [hueSectionA, satSectionA, ~] = createHSVmap(M0_disp_img, maskSection .* maskArtery, 0.15, 0.15);
    [hueSectionV, satSectionV, ~] = createHSVmap(M0_disp_img, maskSection .* maskVein, 0.5, 0.5);
    satSectionArtery = satSectionA;
    satSectionVein = satSectionV; %+~maskSectionArtery.*(~maskArtery);
    val = val .* (~maskSection) + val .* maskSection + maskSection .* (~(maskArtery + maskVein));
    vesselImageRGB = hsv2rgb(hueArtery + hueVein + hueSectionA + hueSectionV, satArtery + satVein + satSectionArtery + satSectionVein, val);

    figure(101)
    imshow(vesselImageRGB)

    %% Create Colormap Artery Only
    [hueArtery, satArtery, val] = createHSVmap(M0_disp_img, maskArtery - maskArtery .* maskSection, 0, 0);
    [hueSectionA, satSectionA, ~] = createHSVmap(M0_disp_img, maskSection .* maskArtery, 0.15, 0.15);
    satSectionArtery = satSectionA;
    val = val .* (~maskSection) + val .* maskSection + maskSection .* (~maskArtery);
    VesselImageRGB_Artery = hsv2rgb(hueArtery + hueSectionA, satArtery + satSectionArtery, val);

    figure(102)
    imshow(VesselImageRGB_Artery)

    %% Create Segmentation Map

    segmentationMap = zeros(numX, numY, 3);
    segmentationMap(:, :, 1) = M0_disp_img - (maskArtery + maskVein) .* M0_disp_img + maskArtery;
    segmentationMap(:, :, 2) = M0_disp_img - (maskArtery + maskVein) .* M0_disp_img;
    segmentationMap(:, :, 3) = M0_disp_img - (maskArtery + maskVein) .* M0_disp_img + maskVein;
    imwrite(segmentationMap, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'arteryVeinSegmentation.png')), 'png');

    %% Saving masks as PNG
    foldername = ToolBox.main_foldername;

    imwrite(mat2gray(double(vesselness_img)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselness.png')), 'png');
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

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_RGVideoVessel.avi')));
    tmp = mat2gray(rgVideoVessel);
    open(w)

    for j = 1:size(rgVideoVessel, 3)
        writeVideo(w, tmp(:, :, j));
    end

    close(w);

end
