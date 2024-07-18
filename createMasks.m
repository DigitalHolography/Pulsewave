function [mask_artery, mask_vein, mask_vessel, maskBackground, maskCRA, maskCRV, maskSection] = createMasks(videoM0, videoM1M0, path, ToolBox)

%% loading parameters and compute useful variables
PW_params = Parameters_json(path);

[Nx, Ny, N_frame] = size(videoM0);

meanIm = squeeze(mean(videoM0, 3));
meanM1M0 = squeeze(mean(videoM1M0, 3));

figure(666), imagesc(meanIm);

for frame_idx = 1:N_frame
    videoM0(:, :, frame_idx) = flat_field_correction(squeeze(videoM0(:, :, frame_idx)), PW_params.flatField_gwRatio * (Ny + Nx) / 2, 0);
end

meanIm = squeeze(mean(videoM0, 3));

figure(667), imagesc(meanIm);

[x, y] = meshgrid(1:Ny, 1:Nx);
cercle_mask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius * (Ny + Nx) / 2;

radius_treshold = PW_params.masks_radius_treshold;
cercle_treshold_mask = [];
min_dist_to_edge=min([ToolBox.x_barycentre,ToolBox.y_barycentre,Nx-ToolBox.x_barycentre,Ny-ToolBox.y_barycentre]);
if radius_treshold == 0
    
    cercle_treshold_mask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= min_dist_to_edge;
elseif  radius_treshold > 0
    cercle_treshold_mask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= (radius_treshold * min_dist_to_edge);
else % if radius_treshold==-1
end

videoM0_zero = videoM0 - meanIm;

% compute vesselness response
vesselnessIm = vesselness_filter(meanIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);

%%  Compute first correlation to find arteries

% compute pulse in 3 dimentions for correlation in all vessels
pulse = squeeze(mean(videoM0 .* (vesselnessIm > 0), [1 2]));
pulse_init = pulse - mean(pulse, "all");
pulse_init_3d = zeros(Nx, Ny, N_frame);

for xx = 1:Nx
    
    for yy = 1:Ny
        pulse_init_3d(xx, yy, :) = pulse_init;
    end
    
end

% compute local-to-average pulse wave zero-lag correlation
correlationMatrix_artery = squeeze(mean((videoM0_zero .* pulse_init_3d), 3)) .* (vesselnessIm > 0);

% Create first artery mask
firstMaskArtery = (correlationMatrix_artery > PW_params.arteryMask_CorrelationMatrixThreshold * mean2(correlationMatrix_artery(correlationMatrix_artery > 0)));

if PW_params.masks_cleaningCoroid
    firstMaskArtery = firstMaskArtery & bwareafilt(firstMaskArtery | cercle_mask, 1, 4);
end


firstMaskArtery = bwareaopen(firstMaskArtery, PW_params.masks_minSize);

if PW_params.masks_showIntermediateFigures
    figure(10), imagesc(vesselnessIm);
    title('Vesselness map');
    colorbar('southoutside');
    figure(11), imshow(firstMaskArtery)
end

clear pulse pulse_init pulse_init_3d correlationMatrix_artery;
%% Compute correlation to segment veins and arteries

% compute pulse in 3 dimentions for correlation in main arteries
pulse = squeeze(mean(videoM0 .* firstMaskArtery, [1 2]));
pulse_init = pulse - mean(pulse, "all");
pulse_init_3d = zeros(Nx, Ny, N_frame);

for xx = 1:Nx
    
    for yy = 1:Ny
        pulse_init_3d(xx, yy, :) = pulse_init;
    end
    
end

% compute local-to-average pulse wave zero-lag correlation
CorrelationMatrix = squeeze(mean((videoM0_zero .* pulse_init_3d), 3));

% compute mean correction correlation to find coroid maximums
meanVideoM0 = zeros(Nx, Ny, N_frame);

for frame_idx = 1:N_frame
    meanVideoM0(:, :, frame_idx) = meanIm;
end

% Create correlation matrix to segment vein and arteries
correlationMatrix_artery = CorrelationMatrix ./ max(CorrelationMatrix, [], 'all');
correlationMatrix_vein = CorrelationMatrix ./ min(CorrelationMatrix, [], 'all');
correlationMatrix_vein(correlationMatrix_vein <- 1) = -1;

if PW_params.masks_showIntermediateFigures
    figure(20), imagesc(correlationMatrix_artery);
    title('Artery correlation map');
    colorbar('southoutside');
    figure(21), imagesc(correlationMatrix_vein);
    title('Vein correlation map');
    colorbar('southoutside');
end

clear pulse pulse_init pulse_init_3d;
%% Region growing to get a clear segmentation

% cleaning correlation matrix
correlationMatrix_artery(vesselnessIm == 0) = min(correlationMatrix_artery, [], 'all');
correlationMatrix_vein(vesselnessIm == 0) = min(correlationMatrix_vein, [], 'all');

% Create vesselness response weighted by correlation
vesselness_artery = vesselnessIm .* rescale(correlationMatrix_artery);
vesselness_vein = vesselnessIm .* rescale(correlationMatrix_vein);

% Create seed artery mask and condition artery mask for region growing
correlationMatrix_arteryCondition = correlationMatrix_artery .* rescale(vesselnessIm);

level_artery = graythresh(correlationMatrix_arteryCondition);
seeds_artery = imbinarize(correlationMatrix_arteryCondition, level_artery ./ PW_params.RG_ArterySeedsThreshold);
seeds_artery = bwareaopen(seeds_artery, PW_params.masks_minSize);

condition_artery = imbinarize(correlationMatrix_arteryCondition, level_artery ./ PW_params.RG_ArteryConditionThreshold);

% Create seed vein mask and contidion vein mask for region growing
correlationMatrix_veinCondition = correlationMatrix_vein .* rescale(vesselnessIm);

level_vein = graythresh(correlationMatrix_veinCondition);
seeds_vein = imbinarize(correlationMatrix_veinCondition, level_vein) | (vesselness_vein > 0.5 * max(vesselness_vein, [], 'all'));
seeds_vein = bwareaopen(seeds_vein, PW_params.masks_minSize);

condition_vein = vesselness_vein > PW_params.RG_veinConditionThreshold * mean2(vesselness_vein(seeds_vein));

% Cleaning condition
condition_artery = condition_artery & bwareafilt(condition_artery | condition_vein | cercle_mask, 1, 4);
condition_artery = bwareaopen(condition_artery, PW_params.masks_minSize);
condition_vein = condition_vein & ~condition_artery;
condition_vein = bwareaopen(condition_vein, PW_params.masks_minSize);

% Cleaning seeds
seeds_artery = seeds_artery & condition_artery;
seeds_vein = seeds_vein & condition_vein;

% Compute region growing to segment vein and arteries
%[mask_artery,RG_video_artery] = region_growing_for_vessel(vesselness_artery, seeds_artery, condition_artery,path);
%[mask_vein,RG_video_vein]  = region_growing_for_vessel(vesselness_vein, seeds_vein, condition_vein, path);
[mask_vessel, RG_video_vessel] = region_growing_for_vessel(vesselnessIm, seeds_artery | seeds_vein, condition_vein | condition_artery, path);


if isfile(fullfile(ToolBox.PW_path_main,'mask', 'maskVessel_New.png')) % import manually tuned mask if needed
    mask_vessel = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main,'mask', 'maskVessel_New.png')),3))>0;
else
    mask_vessel = bwareafilt(mask_vessel, PW_params.arteryMask_magicwand_nb_of_area_vessels, 4);
end

if isfile(fullfile(ToolBox.PW_path_main,'mask', '_maskArtery_New.png')) % import manually tuned mask if needed
    mask_artery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main,'mask', '_maskArtery_New.png')),3))>0;
else
    mask_artery = mask_vessel .* correlationMatrix_artery;
    mask_artery = mask_artery > PW_params.arteryMask_ArteryCorrThreshold;
end

if isfile(fullfile(ToolBox.PW_path_main,'mask', '_maskVein_New.png')) % import manually tuned mask if needed
    mask_vein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main,'mask', '_maskVein_New.png')),3))>0;
else
    mask_vein = mask_vessel .* correlationMatrix_vein;
    mask_vein = mask_vein > 0;
end
%mask_vein = mask_vessel & ~mask_artery;


if PW_params.masks_showIntermediateFigures
    figure(30), imagesc(vesselness_artery);
    title('Artery Vesselness map');
    colorbar('southoutside');
    figure(31), imagesc(vesselness_vein);
    title('Vein Vesselness map');
    colorbar('southoutside');
    figure(32), imagesc(uint8(cat(3, uint8(meanIm) + uint8(seeds_artery) * 255, uint8(meanIm), uint8(meanIm) + uint8(seeds_vein) * 255)));
    title('Artery/Vein initialisation for region growing');
    figure(33), imagesc(uint8(cat(3, uint8(meanIm) + uint8(condition_artery) * 255, uint8(meanIm), uint8(meanIm) + uint8(condition_vein) * 255)));
    title('Artery/Vein condition for region growing');
    figure(34), imshow(mask_vessel);
end

clear vesselness_vein vesselness_artery floor_vein floor_artery level_vein level_artery seeds_vein seeds_artery ;
%% Cleaning coroid from masks

if ~isempty(cercle_treshold_mask)
    mask_artery = mask_artery & cercle_treshold_mask;
end

% if PW_params.masks_cleaningCoroid
%
%     mask_artery = mask_artery & bwareafilt(mask_artery | mask_vein | cercle_mask,1,4);
%     mask_vein = mask_vein & bwareafilt(mask_artery | mask_vein | cercle_mask,1,4);
%
%     if PW_params.masks_showIntermediateFigures
%     figure(40), imagesc(uint8( cat( 3, uint8(meanIm)+ uint8(mask_artery)*255, uint8(meanIm) + uint8(cercle_mask) , uint8(meanIm) + uint8(mask_vein)*255 )));
%     title('Artery/Vein cleaned segmentation');
%     end
%
%     clear cercle_mask;
% end



mask_artery = bwareaopen(mask_artery, PW_params.masks_minSize);
mask_artery = imdilate(mask_artery, strel('disk', 3));
mask_artery = imclose(mask_artery, strel('disk', 5));

mask_artery = bwareafilt(mask_artery, PW_params.arteryMask_magicwand_nb_of_area_artery, 4);

mask_vein = bwareaopen(mask_vein, PW_params.masks_minSize);
mask_vein = imdilate(mask_vein, strel('disk', 3));
mask_vein = imclose(mask_vein, strel('disk', 5)) & ~mask_artery;

mask_vessel = mask_artery | mask_vein;




if PW_params.masks_showIntermediateFigures
    figure(40), imagesc(uint8(cat(3, uint8(meanIm) + uint8(mask_artery) * 255, uint8(meanIm), uint8(meanIm) + uint8(mask_vein) * 255)));
    title('Artery/Vein region growing segmentation');
end

%% Create CRA and CRV Mask

stdM1M0 = std2(meanM1M0);
maskCRA = meanM1M0 > (PW_params.CRACRV_Threshold * stdM1M0);
maskCRV = meanM1M0 < (-PW_params.CRACRV_Threshold * stdM1M0);

%% Creat Background Mask

maskBackground = not(mask_vessel);

%% Create Mask Section

radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (Ny + Nx) / 2;
radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (Ny + Nx) / 2;

circle_mask1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius1;
circle_mask2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius2;

maskSection = xor(circle_mask1, circle_mask2);

%% Create Colormap Artery/Vein

meanIm = mat2gray(meanIm);
[hue_artery, sat_artery, val] = createHSVmap(meanIm, mask_artery - mask_artery .* maskSection, 0, 0);
[hue_vein, sat_vein, ~] = createHSVmap(meanIm, mask_vein - mask_vein .* maskSection - mask_vein .* mask_artery, 0.7, 0.7);
[hue_sectionA, sat_sectionA, ~] = createHSVmap(meanIm, maskSection .* mask_artery, 0.15, 0.15);
[hue_sectionV, sat_sectionV, ~] = createHSVmap(meanIm, maskSection .* mask_vein, 0.5, 0.5);
sat_section_artery = sat_sectionA;
sat_section_vein = sat_sectionV; %+~maskSectionArtery.*(~mask_artery);
val = val .* (~maskSection) + val .* maskSection + maskSection .* (~(mask_artery + mask_vein));
VesselImageRGB = hsv2rgb(hue_artery + hue_vein + hue_sectionA + hue_sectionV, sat_artery + sat_vein + sat_section_artery + sat_section_vein, val);

figure(101)
imshow(VesselImageRGB)

%% Create Colormap Artery Only
[hue_artery, sat_artery, val] = createHSVmap(meanIm, mask_artery - mask_artery .* maskSection, 0, 0);
[hue_sectionA, sat_sectionA, ~] = createHSVmap(meanIm, maskSection .* mask_artery, 0.15, 0.15);
sat_section_artery = sat_sectionA;
val = val .* (~maskSection) + val .* maskSection + maskSection .* (~mask_artery);
VesselImageRGB_Artery = hsv2rgb(hue_artery + hue_sectionA, sat_artery + sat_section_artery, val);

figure(102)
imshow(VesselImageRGB_Artery)

%% Saving masks as PNG
foldername = ToolBox.main_foldername;

imwrite(mat2gray(double(vesselnessIm)), fullfile(ToolBox.PW_path_png, [foldername, '_vesselness.png']), 'png');
imwrite(mat2gray(single(mask_artery)), fullfile(ToolBox.PW_path_png, [foldername, '_maskArtery_New.png']), 'png');
imwrite(mat2gray(single(mask_vein)), fullfile(ToolBox.PW_path_png, [foldername, '_maskVein_New.png']), 'png');
imwrite(mat2gray(single(mask_vessel)), fullfile(ToolBox.PW_path_png, [foldername, '_maskVessel_New.png']), 'png');
imwrite(mat2gray(single(maskBackground)), fullfile(ToolBox.PW_path_png, [foldername, '_maskBackground_New.png']), 'png');
%vesselMap = uint8( cat( 3, uint8(meanIm)+ uint8(mask_artery)*255, uint8(meanIm) , uint8(meanIm) + uint8(mask_vein)*255 ));
imwrite(VesselImageRGB, fullfile(ToolBox.PW_path_png, [foldername, '_vesselMap.png']), 'png');
imwrite(VesselImageRGB_Artery, fullfile(ToolBox.PW_path_png, [foldername, '_vesselMapArtery.png']), 'png');
imwrite(mat2gray(maskCRA), fullfile(ToolBox.PW_path_png, [foldername, '_maskCRA.png']), 'png');
imwrite(mat2gray(maskCRV), fullfile(ToolBox.PW_path_png, [foldername, '_maskCRV.png']), 'png');

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
tmp = mat2gray(RG_video_vessel);
open(w)

for j = 1:size(RG_video_vessel, 3)
    writeVideo(w, tmp(:, :, j));
end

close(w);

end
