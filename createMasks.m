function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(M0_disp_video, ~, f_AVG_mean, path, ToolBox)

%% 0) loading parameters and compute useful variables
PW_params = Parameters_json(path);
mkdir(ToolBox.PW_path_png, 'mask')
mkdir(fullfile(ToolBox.PW_path_png, 'mask'), 'step')
close all

[numX, numY, ~] = size(M0_disp_video);
[X, Y] = meshgrid(1:numX, 1:numY);
maskDiaphragm = sqrt((X - numX / 2) .^ 2 + (Y - numY / 2) .^ 2) <= 0.45 * (numY + numX) / 2;

M0_disp_img = squeeze(mean(M0_disp_video, 3));
M0_disp_centered_video = M0_disp_video - M0_disp_img;
imwrite(rescale(M0_disp_img), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'M0_disp_img.png')))

%% 1) First Masks and Correlation

%% 1) 1) Compute vesselness response
vesselnessM0 = vesselness_filter(M0_disp_img, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
maskVesselness = logical(imbinarize(vesselnessM0) .* maskDiaphragm);

imwrite(rescale(vesselnessM0), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_1_1_Vesselness.png')))
imwrite(maskVesselness, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_1_2_vesselMask.png')))

%%  1) 2) Compute first correlation to find arteries

% compute pulse in 3 dimentions for correlation in all vessels
vascularPulse = mean(M0_disp_video .* maskVesselness, [1 2]);
vascularPulse = vascularPulse ./ nnz(maskVesselness);
vascularPulse_centered = vascularPulse - mean(vascularPulse, 3);

% compute local-to-average pulse wave zero-lag correlation
R_VascularPulse = mean(M0_disp_centered_video .* vascularPulse_centered, 3) ./ (std((M0_disp_centered_video), [], 3) * std(vascularPulse_centered, [], 3));

% Create first artery mask
figure, histogram(R_VascularPulse(R_VascularPulse ~= 0))
title("Histogram of the Vascular Pulse Correlation")

firstThresholds = multithresh(R_VascularPulse(maskVesselness), 3);
firstThresholds = [-1 firstThresholds];
quantizedVesselCorrelation = imquantize(R_VascularPulse - ~maskVesselness * 2, firstThresholds);
% 1 = Background, 2,3 = Veins & CoroidalVessels, 4 = CoroidalVessel, 5 = Arteries

firstMaskArtery = quantizedVesselCorrelation == 5;
firstMaskVein = quantizedVesselCorrelation == 3 | quantizedVesselCorrelation == 2;

firstMaskArtery = bwareaopen(firstMaskArtery, PW_params.masks_minSize);
firstMaskVein = bwareaopen(firstMaskVein, PW_params.masks_minSize);

imwrite(firstMaskArtery, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_1_First.png')))
imwrite(firstMaskVein, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_1_First.png')))
imwrite(rescale(R_VascularPulse), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_1_vesselCorrelation.png')))
imwrite(rescale(quantizedVesselCorrelation), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_1_quantizedVesselCorrelation.png')))

clear vascularPulse quantizedVesselCorrelation firstThresholds;
%% 2) 1) Compute correlation to segment veins and arteries

% compute pulse in 3 dimentions for correlation in main arteries
arterialPulse = mean(M0_disp_video .* firstMaskArtery .* maskDiaphragm, [1 2]);
arterialPulse = arterialPulse ./ nnz(firstMaskArtery .* maskDiaphragm);
arterialPulse_centered = arterialPulse - mean(arterialPulse, 3);

% compute pulse in 3 dimentions for correlation in main arteries
venousPulse = mean(M0_disp_video .* firstMaskVein, [1 2]);
venousPulse = venousPulse ./ nnz(firstMaskVein);
venousPulse_centered = venousPulse - mean(venousPulse, 3);

% compute local-to-average pulse wave zero-lag correlation
R_ArterialPulse = mean((M0_disp_centered_video .* arterialPulse_centered), 3) ./ (std((M0_disp_centered_video), [], 3) * std(arterialPulse_centered, [], 3));
R_VenousPulse = mean((M0_disp_centered_video .* venousPulse_centered), 3) ./ (std((M0_disp_centered_video), [], 3) * std(arterialPulse_centered, [], 3));

imwrite(rescale(R_ArterialPulse), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_2_CorrelMatrix.png')))
imwrite(rescale(R_VenousPulse), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_2_CorrelMatrix.png')))

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

%% 3) 2) Segmentation

% New Vesselness for arteries as they thend to have a larger caliber
vesselness4Artery = vesselness_filter(R_ArterialPulse, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta) .* maskDiaphragm;
maskVessel4Artery = logical(imbinarize(vesselness4Artery) .* maskDiaphragm);

% Choroidal Vessels Cleaning
se = strel('disk', 2);
maskVessel4Artery = maskVessel4Artery & bwareafilt(imopen(maskVessel4Artery, se) | cercleMask, 1, 4);
maskVessel4Vein = maskVesselness & bwareafilt(imopen(maskVesselness, se) | cercleMask, 1, 4);

imwrite(rescale(maskVessel4Vein), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_3_1_Vesselness4Vein.png')))
imwrite(rescale(vesselness4Artery), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_1_Vesselness4Artery.png')))

% Correlation Matrix times vesselness to increase contrast in vessels
R_ArteryVessel = R_ArterialPulse .* vesselness4Artery .* maskVessel4Artery;
R_VeinVessel = R_VenousPulse .* vesselnessM0 .* maskVesselness;

imwrite(rescale(R_ArteryVessel), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_1_R_ArteryVessel.png')))
imwrite(rescale(R_VeinVessel), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_1_R_VeinVessel.png')))

% Segmentation based on Otsu's criteria for arteries and veins
% Number of classes for Arteries: 3
% 1: Background
% 2: Other Vessels
% 3: Arteries
numClassesArtery = 3;

levelArtery = multithresh(R_ArteryVessel(maskVessel4Artery), numClassesArtery - 2);
levelArtery = [-1 levelArtery];
maskArteryQ = imquantize(R_ArteryVessel - 2 * ~maskVessel4Artery, levelArtery);
maskArtery = maskArteryQ == max(maskArteryQ, [], 'all');

imwrite(rescale(maskArteryQ), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_2_Quantize.png')))
imwrite(rescale(maskArtery), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_2_Thresh.png')))

% Number of classes for Veins: 4
% 1: Background
% 2: Arteries
% 3: Choroidal Vessels
% 3: Veins
numClassesVein = 4;

levelVein = multithresh(R_VeinVessel(maskVesselness), numClassesVein - 2);
levelVein = [-1 levelVein];
maskVeinQ = imquantize(R_VeinVessel - 2 * ~maskVesselness, levelVein);
maskVein = maskVeinQ == max(maskVeinQ, [], 'all');

imwrite(rescale(maskVeinQ), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_2_Quantize.png')))
imwrite(rescale(maskVein), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_2_Thresh.png')))

%% Clearing 3) 3)

% Remove the small bits
maskArteryClean = bwareaopen(maskArtery, PW_params.masks_minSize, 4);
maskVeinClean = bwareaopen(maskVein, PW_params.masks_minSize, 4);

imwrite(rescale(maskArteryClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_3_AreaOpen.png')))
imwrite(rescale(maskVeinClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_3_AreaOpen.png')))

% Remove the small squares from the correlation
maskArteryClean = imclose(maskArteryClean, se);
maskVeinClean = imclose(maskVeinClean, se);

imwrite(rescale(maskArteryClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_3_ImClose.png')))
imwrite(rescale(maskVeinClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_3_ImClose.png')))

% maskArtery = bwareafilt(maskArtery, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);
% maskArteryClean = bwareafilt(maskArteryClean, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);

RGBM0(:, :, 1) = rescale(M0_disp_img) + maskArteryClean;
RGBM0(:, :, 2) = rescale(M0_disp_img);
RGBM0(:, :, 3) = rescale(M0_disp_img) + maskVeinClean;

imwrite(RGBM0, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_3_3_RGB.png')))

clear vesselnessVein vesselnessArtery ;
%% 3) 4) Cleaning Choroid

% if ~isempty(cercleTresholdMask)
%     maskArteryClean = maskArteryClean & cercleTresholdMask;
%     maskVeinClean = maskVeinClean & cercleTresholdMask;
%     imwrite(rescale(maskArteryClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, '3_maskArtery1.png')))
%     imwrite(rescale(maskVeinClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, '3_maskArtery1.png')))
% end

% Remove the small squares from the correlation
se = strel('disk', 3);
maskArteryClean = imclose(maskArteryClean, se);
maskVeinClean = imclose(maskVeinClean, se);

imwrite(rescale(maskArteryClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_3_ImDilate.png')))
imwrite(rescale(maskVeinClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_3_ImDilate.png')))

% maskArteryClean = bwareafilt(maskArteryClean, PW_params.arteryMask_magicwand_nb_of_area_artery, 4);

% if PW_params.masks_cleaningCoroid
%
%     maskArteryClean = maskArteryClean & bwareafilt(maskArteryClean | maskVeinClean | cercleMask, 1, 4);
%     maskVeinClean = maskVeinClean & bwareafilt(maskArteryClean | maskVeinClean | cercleMask, 1, 4);
%
%     imwrite(rescale(maskArteryClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_3_4_CleaningChoroid.png')))
%     imwrite(rescale(maskVeinClean), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_3_4_CleaningChoroid.png')))
%     imwrite(rescale(uint8( cat( 3, uint8(M0_disp_img)+ uint8(maskArtery)*255, uint8(M0_disp_img) + uint8(cercleMask) , uint8(M0_disp_img) + uint8(maskVein)*255 ))), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_3_4_CleaningChoroid.png')))
%
%     clear cercleMask;
% end

maskArtery = maskArteryClean;
maskVein = maskVeinClean;
maskVessel = maskArtery | maskVein;

%% SECRET CODE BIT: SENSITIVITY & SPECIFICITY
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
catch ME
    disp(ME)
end

%% Export Masks as the final step

imwrite(maskArtery, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskArtery_4_Final.png')))
imwrite(maskVein, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVein_4_Final.png')))
imwrite(maskVessel, fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_4_Final.png')))
imwrite(rescale(uint8(cat(3, uint8(M0_disp_img) + uint8(maskArtery) * 255, uint8(M0_disp_img), uint8(M0_disp_img) + uint8(maskVein) * 255))), fullfile(ToolBox.PW_path_png, 'mask', 'step', sprintf("%s_%s", ToolBox.main_foldername, 'maskVessel_4_ArteryVeinRGSegmentation.png')))

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

imwrite(mat2gray(double(vesselnessM0)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'vesselness.png')), 'png');
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
