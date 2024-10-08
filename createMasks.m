function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(M0_disp_video, ~, f_AVG_mean, path, ToolBox)

%% 0) loading parameters and compute useful variables
PW_params = Parameters_json(path);
mkdir(ToolBox.PW_path_png, 'mask')
close all

[numX, numY, ~] = size(M0_disp_video);
[X, Y] = meshgrid(1:numX, 1:numY);
maskDiaphragm = sqrt((X - numX/2) .^ 2 + (Y - numY/2) .^ 2) <= 0.4 * (numY + numX) / 2;

M0_disp_img = squeeze(mean(M0_disp_video, 3));
M0_disp_centered_video = M0_disp_video - M0_disp_img;
figure(1), imagesc(M0_disp_img), axis image, colormap gray;

%% 1) First Masks and Correlation

%% 1) 1) Compute vesselness response
M0_contrasted_img = adapthisteq(rescale(M0_disp_img), 'NumTiles', PW_params.arteryMask_vesselnessContrastNumSlides');
vesselness_img = vesselness_filter(M0_contrasted_img, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);

maskVessel = imbinarize(vesselness_img);
figure(10), imagesc(vesselness_img), axis image, colormap gray;
imwrite(maskVessel, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'vesselnessBinaryMask.png')))

clear M0_contrasted_img

%%  1) 2) Compute first correlation to find arteries

% compute pulse in 3 dimentions for correlation in all vessels
vascularPulse = mean(M0_disp_video .* maskVessel .* maskDiaphragm, [1 2]);
vascularPulse = vascularPulse ./ nnz(maskVessel .* maskDiaphragm);
vascularPulse_centered = vascularPulse - mean(vascularPulse, 3);

% compute local-to-average pulse wave zero-lag correlation
R_VascularPulse = mean(M0_disp_centered_video .* vascularPulse_centered .* maskDiaphragm , 3) ./ (std((M0_disp_centered_video), [], 3) * std(vascularPulse_centered, [], 3));
vesselCorrelation = R_VascularPulse .* maskVessel;

% Create first artery mask
figure, histogram(vesselCorrelation(vesselCorrelation ~= 0))
title("Histogram of the Vascular Pulse Correlation")

firstThresholds = multithresh(vesselCorrelation, 3);
quantizedVesselCorrelation = imquantize(vesselCorrelation, firstThresholds); % 1 = Veins, 2 = Background, 3 = CoroidalVessels, 4 = Arteries
firstMaskArtery = quantizedVesselCorrelation == 4;
if min(vesselCorrelation, [], 'all') == 0
    firstMaskVein = quantizedVesselCorrelation == 2;
else
    firstMaskVein = quantizedVesselCorrelation == 1;
end

firstMaskArtery = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVein = bwareaopen(firstMaskVein, PW_params.masks_minSize);

figure(11), imagesc(vesselness_img), axis image, title('Vesselness map'), colorbar('southoutside');
figure(12), imagesc(firstMaskArtery), axis image, title('First Artery Mask'), colorbar('southoutside');
figure(13), imagesc(firstMaskVein), axis image, title('First Vein Mask'), colorbar('southoutside');
figure(14), imagesc(vesselCorrelation), axis image, title('Vessel Correlation Mask'), colorbar('southoutside');

clear vascularPulse R_VascularPulse quantizedVesselCorrelation firstThresholds;
%% Compute correlation to segment veins and arteries

% compute pulse in 3 dimentions for correlation in main arteries
arterialPulse = mean(M0_disp_video .* firstMaskArtery .* maskVessel .* maskDiaphragm, [1 2]);
arterialPulse = arterialPulse ./ nnz(firstMaskArtery .* maskVessel .* maskDiaphragm);
arterialPulse_centered = arterialPulse - mean(arterialPulse, 3);

% compute pulse in 3 dimentions for correlation in main arteries
venousPulse = mean(M0_disp_video .* firstMaskVein, [1 2]);
venousPulse = venousPulse ./ nnz(firstMaskVein);
venousPulse_centered = venousPulse - mean(venousPulse, 3);

% compute local-to-average pulse wave zero-lag correlation
R_ArterialPulse = mean((M0_disp_centered_video .* arterialPulse_centered .* maskDiaphragm .* maskVessel), 3) ./ (std((M0_disp_centered_video), [], 3) * std(arterialPulse_centered, [], 3));
R_VenousPulse = mean((M0_disp_centered_video .* venousPulse_centered .* maskDiaphragm .* maskVessel), 3) ./ (std((M0_disp_centered_video), [], 3) * std(arterialPulse_centered, [], 3));

figure(20), imagesc(R_ArterialPulse);
axis image
title('Arterial Pulse correlation map');
colorbar('southoutside');
imwrite(rescale(R_ArterialPulse), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrixArtery.png')))

figure(21), imagesc(R_VenousPulse);
axis image
title('Venous Pulse correlation map');
colorbar('southoutside');
imwrite(rescale(R_VenousPulse), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'correlationMatrixVein.png')))

clear arterialPulse arterialPulse_centered venousPulse venousPulse_centered;

%% 3) Circles Sectionning

%% 3) 1) Circles Sectioning

blurred_mask = imgaussfilt(double(M0_disp_img .* f_AVG_mean) .* R_ArterialPulse, PW_params.gauss_filt_size_for_barycentre * numX, 'Padding', 0);
[ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));

cercleMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (numY + numX) / 2;

radiusTreshold = PW_params.masks_radius_treshold;
cercleTresholdMask = [];
minDistToEdge = min([ToolBox.x_barycentre, ToolBox.y_barycentre, numX - ToolBox.x_barycentre, numY - ToolBox.y_barycentre]);

if radiusTreshold == 0
    cercleTresholdMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= minDistToEdge;
elseif radiusTreshold > 0
    cercleTresholdMask = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= (radiusTreshold * minDistToEdge);
else % if radius_treshold==-1
end

%% 3) 2) Create CRA and CRV Mask

f_AVG_std = std2(f_AVG_mean);
maskCRA = f_AVG_mean > (PW_params.CRACRV_Threshold * f_AVG_std);
maskCRV = f_AVG_mean < (-PW_params.CRACRV_Threshold * f_AVG_std);
clear f_AVG_std f_AVG_mean

%% 4) Region growing to get a clear segmentation

% Create vesselness response weighted by correlation
vesselnessArtery = vesselness_img .* rescale(R_ArterialPulse);
vesselnessVein = vesselness_img .* rescale(R_VenousPulse);

% Create seed artery mask and condition artery mask for region growing
R_ArteryCondition = R_ArterialPulse .* rescale(vesselness_img);
levelArtery = multithresh(R_ArteryCondition, 3);
maskArtery = imquantize(R_ArteryCondition, levelArtery);
seedsArtery = maskArtery == 4;
seedsArtery = bwareaopen(seedsArtery, PW_params.masks_minSize);

conditionArtery = imbinarize(R_ArterialPulse, levelArtery(end) ./ PW_params.RG_ArteryConditionThreshold);

% Create seed vein mask and contidion vein mask for region growing
R_VeinCondition = R_VenousPulse .* rescale(vesselness_img);

levelVein = multithresh(R_VeinCondition, 2);
maskVein = imquantize(R_VeinCondition, levelVein);
seedsVein = maskVein == 3 | (vesselnessVein > 0.5 * max(vesselnessVein, [], 'all'));
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
[maskArtery, ~] = region_growing_for_vessel(vesselnessArtery, seedsArtery, conditionArtery,path);
[maskVein, ~]  = region_growing_for_vessel(vesselnessVein, seedsVein, conditionVein, path);
% [maskVessel, ~] = region_growing_for_vessel(vesselness_img, seedsArtery, conditionArtery, path);

maskArtery = bwareafilt(maskArtery, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);
maskVein = bwareafilt(maskVein, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);


%maskVein = maskVessel & ~maskArtery;

if PW_params.masks_showIntermediateFigures
    figure(30), imagesc(vesselnessArtery);
    axis image
    title('Artery Vesselness map');
    colorbar('southoutside');
    figure(31), imagesc(vesselnessVein);
    axis image
    title('Vein Vesselness map');
    colorbar('southoutside');
    figure(32), imagesc(uint8(cat(3, uint8(M0_disp_img) + uint8(seedsArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(seedsVein) * 255)));
    axis image
    title('Artery/Vein initialisation for region growing');
    figure(33), imagesc(uint8(cat(3, uint8(M0_disp_img) + uint8(conditionArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(conditionVein) * 255)));
    axis image
    title('Artery/Vein condition for region growing');
    figure(34), imagesc(maskVessel);
    axis image
end

clear vesselnessVein vesselnessArtery levelVein levelArtery seedsVein seedsArtery ;
%% Cleaning coroid from masks

if ~isempty(cercleTresholdMask)
    maskArtery = maskArtery & cercleTresholdMask;
end

if PW_params.masks_cleaningCoroid

    maskArtery = maskArtery & bwareafilt(maskArtery | maskVein | cercleMask, 1, 4);
    maskVein = maskVein & bwareafilt(maskArtery | maskVein | cercleMask, 1, 4);

    if PW_params.masks_showIntermediateFigures
        figure(40), imagesc(uint8( cat( 3, uint8(M0_disp_img)+ uint8(maskArtery)*255, uint8(M0_disp_img) + uint8(cercleMask) , uint8(M0_disp_img) + uint8(maskVein)*255 )));
        title('Artery/Vein cleaned segmentation');
    end

    clear cercleMask;
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

circleMask1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius1;
circleMask2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= radius2;

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
