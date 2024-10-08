function [mask_artery, mask_vein, mask_vessel,maskBackground,maskCRA,maskCRV,maskSectionArtery] = createMasksNewOld(videoM0,videoM1M0, path, ToolBox)

%% loading parameters and compute useful variables
PW_params = Parameters_json(path);

[N,M,L] = size(videoM0);

for pp = 1:L
    videoM0(:,:,pp) = flat_field_correction(squeeze(videoM0(:,:,pp)), PW_params.flatField_gwRatio*N, PW_params.flatField_border);
end

meanIm = squeeze(mean(videoM0, 3));
meanM1M0 = squeeze(mean(videoM1M0, 3));

[x, y] = meshgrid(1:M,1:N);
cercle_mask = sqrt((x - ToolBox.x_barycentre).^2 + (y - ToolBox.y_barycentre).^2) <= PW_params.masks_radius*(M+N)/2;

videoM0_zero = videoM0 - meanIm;

% compute vesselness response
vesselnessIm = vesselness_filter(meanIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);

%% Create CRA and CRV Mask

stdM1M0 = std2(meanM1M0);
maskCRA = meanM1M0>(PW_params.CRACRV_Threshold*stdM1M0);
maskCRV = meanM1M0<(-PW_params.CRACRV_Threshold*stdM1M0);


%%  Compute first correlation to find arteries

% compute pulse in 3 dimentions for correlation in all vessels
pulse = squeeze(mean(videoM0 .* (vesselnessIm>0), [1 2]));
pulse_init = pulse - mean(pulse, "all");
pulse_init_3d = zeros(N,M,L);
for nn = 1:N
    for mm = 1:M
        pulse_init_3d(nn,mm,:) = pulse_init;
    end
end

% compute local-to-average pulse wave zero-lag correlation
correlationMatrix_artery = squeeze(mean((videoM0_zero .* pulse_init_3d), 3)).*(vesselnessIm>0);

% Create first artery mask
firstMaskArtery = (correlationMatrix_artery >    1.5*mean2(correlationMatrix_artery(correlationMatrix_artery>0)));
if PW_params.masks_cleaningCoroid 
    firstMaskArtery = firstMaskArtery & bwareafilt(firstMaskArtery | cercle_mask,1,4);
end
firstMaskArtery = bwareaopen(firstMaskArtery,PW_params.masks_minSize);

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
pulse_init_3d = zeros(N,M,L);
for nn = 1:N
    for mm = 1:M
        pulse_init_3d(nn,mm,:) = pulse_init;
    end
end

% compute local-to-average pulse wave zero-lag correlation
CorrelationMatrix = squeeze(mean((videoM0_zero.* pulse_init_3d), 3));

% compute mean correction correlation to find coroid maximums
meanVideoM0 = zeros(N,M,L);
for ll = 1:L
    meanVideoM0(:,:,L) = meanIm;
end

% Create correlation matrix to segment vein and arteries
correlationMatrix_artery = CorrelationMatrix./max(CorrelationMatrix,[],'all');
correlationMatrix_vein = CorrelationMatrix./min(CorrelationMatrix,[],'all');
correlationMatrix_vein(correlationMatrix_vein<-1) = -1;

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
correlationMatrix_artery(vesselnessIm == 0) = min(correlationMatrix_artery,[],'all');
correlationMatrix_vein(vesselnessIm == 0) = min(correlationMatrix_vein,[],'all');

% Create vesselness response weighted by correlation
vesselness_artery = vesselnessIm.*rescale(correlationMatrix_artery);
vesselness_vein = vesselnessIm.*rescale(correlationMatrix_vein);

% Create seed artery mask and contidion artery mask for region growing
correlationMatrix_arteryCondition = correlationMatrix_artery.*rescale(vesselnessIm);

level_artery = graythresh(correlationMatrix_arteryCondition);
seeds_artery = imbinarize(correlationMatrix_arteryCondition, level_artery);
seeds_artery = bwareaopen(seeds_artery,PW_params.masks_minSize);

condition_artery = imbinarize(correlationMatrix_arteryCondition,level_artery./PW_params.RG_ArteryConditionThreshold);

% Create seed vein mask and contidion vein mask for region growing
correlationMatrix_veinCondition = correlationMatrix_vein.*rescale(vesselnessIm);

level_vein = graythresh(correlationMatrix_veinCondition);
seeds_vein = imbinarize(correlationMatrix_veinCondition, level_vein)  | (vesselness_vein > 0.5*max(vesselness_vein,[],'all'));
seeds_vein = bwareaopen(seeds_vein,PW_params.masks_minSize);

condition_vein = vesselness_vein > PW_params.RG_veinConditionThreshold*mean2(vesselness_vein(seeds_vein)); 

% Cleaning condition
condition_artery = condition_artery & bwareafilt(condition_artery | condition_vein |cercle_mask,1,4);
condition_artery = bwareaopen(condition_artery,PW_params.masks_minSize);
condition_vein = condition_vein & ~condition_artery;
condition_vein = bwareaopen(condition_vein,PW_params.masks_minSize);

% Cleaning seeds
seeds_artery = seeds_artery & condition_artery;
seeds_vein = seeds_vein & condition_vein;

% Compute region growing to segment vein and arteries
%[mask_artery,RG_video_artery] = region_growing_for_vessel(vesselness_artery, seeds_artery, condition_artery,path);
%[mask_vein,RG_video_vein]  = region_growing_for_vessel(vesselness_vein, seeds_vein, condition_vein, path);
[mask_vessel,RG_video_vessel] = region_growing_for_vessel(vesselnessIm, seeds_artery | seeds_vein, condition_vein | condition_artery,path);

% mask_artery = mask_vessel.*correlationMatrix_artery;
% mask_artery = mask_artery>PW_params.arteryMask_ArteryCorrThreshold;
mask_artery = condition_artery;
mask_artery = imdilate(mask_artery,strel('disk',1));
mask_artery = bwareaopen(mask_artery,PW_params.masks_minSize);
mask_artery = imdilate(mask_artery,strel('disk',2));

%mask_vein = imdilate(mask_vessel,strel('disk',2)) & ~mask_artery& ~maskCRA;
mask_vein = imdilate(mask_vessel,strel('disk',2)) & condition_vein;
% mask_vein = mask_vessel.*correlationMatrix_vein;
% mask_vein = mask_vein > 0;
mask_vein = bwareaopen(mask_vein,PW_params.masks_minSize);
%mask_vein = imdilate(mask_vein,strel('disk',2));

if PW_params.masks_showIntermediateFigures
    figure(30), imagesc(vesselness_artery);
    title('Artery Vesselness map');
    colorbar('southoutside');
    figure(31), imagesc(vesselness_vein);
    title('Vein Vesselness map');
    colorbar('southoutside');
    figure(32), imagesc( uint8( cat( 3,uint8(meanIm) + uint8(seeds_artery)*255, uint8(meanIm) , uint8(meanIm) + uint8(seeds_vein)*255 )));
    title('Artery/Vein initialisation for region growing');
    figure(33), imagesc( uint8( cat( 3, uint8(meanIm) + uint8(condition_artery)*255, uint8(meanIm)  , uint8(meanIm) + uint8(condition_vein)*255 )));
    title('Artery/Vein condition for region growing')
    figure(34), imagesc(uint8( cat( 3, uint8(meanIm)+ uint8(mask_artery)*255, uint8(meanIm) , uint8(meanIm) + uint8(mask_vein)*255 )));
    title('Artery/Vein region growing segmentation');
    figure(35), imshow(mask_vessel);
end

clear vesselness_vein vesselness_artery floor_vein floor_artery level_vein level_artery seeds_vein seeds_artery ;
%% Cleaning coroid from masks

% if PW_params.masks_cleaningCoroid 
% 
%     mask_artery = mask_artery & bwareafilt(mask_artery | mask_vein | cercle_mask,1,4);
%     mask_artery = imclose(mask_artery,strel('disk',5));
% 
%     mask_vein = mask_vein & bwareafilt(mask_artery | mask_vein | cercle_mask,1,4);
%     mask_vein = imclose(mask_vein,strel('disk',5));
% 
%     mask_vessel = mask_artery | mask_vein;
% 
%     if PW_params.masks_showIntermediateFigures
%     figure(40), imagesc(uint8( cat( 3, uint8(meanIm)+ uint8(mask_artery)*255, uint8(meanIm) + uint8(cercle_mask) , uint8(meanIm) + uint8(mask_vein)*255 )));
%     title('Artery/Vein cleaned segmentation');
%     end
% 
%     clear cercle_mask;
% end
%% Creat Background Mask

maskBackground = not(mask_vessel);
%% Create Mask Section 

ecart = 0.01;

radius1 = (PW_params.radius_ratio-ecart)* (M+N)/2;
radius2 = (PW_params.radius_ratio+ecart)* (M+N)/2;

cercle_mask1 = sqrt((x - ToolBox.x_barycentre).^2 + (y - ToolBox.y_barycentre).^2) <= radius1;
cercle_mask2 = sqrt((x - ToolBox.x_barycentre).^2 + (y - ToolBox.y_barycentre).^2) <= radius2;

maskSectionArtery = xor(cercle_mask1,cercle_mask2);


%% Create Colormap ARtery/Vein 
meanIm = mat2gray(meanIm);
[hue_artery,sat_artery,val] = createHSVmap(meanIm,mask_artery-mask_artery.*maskSectionArtery,0,0);
[hue_vein,sat_vein,~] = createHSVmap(meanIm,mask_vein-mask_vein.*maskSectionArtery-mask_vein.*mask_artery,0.7,0.7);
[hue_section,sat_section,~] = createHSVmap(meanIm,maskSectionArtery,0.35,0.35);
sat_section = sat_section.*mask_artery ; %+~maskSectionArtery.*(~mask_artery);
val = val.*(~maskSectionArtery)+val.*maskSectionArtery+ maskSectionArtery.*(~mask_artery);
VesselImageRGB =  hsv2rgb(hue_artery+hue_vein+hue_section, sat_artery+sat_vein+sat_section, val);

figure(15)
imshow(VesselImageRGB)


%% Saving masks as PNG 
foldername = ToolBox.main_foldername;

imwrite(mat2gray(double(vesselnessIm)),fullfile(ToolBox.PW_path_png,[foldername,'_vesselness.png']),'png') ;
imwrite(mat2gray(single(mask_artery)),fullfile(ToolBox.PW_path_png,[foldername,'_maskArtery_New.png']),'png') ;
imwrite(mat2gray(single(mask_vein)),fullfile(ToolBox.PW_path_png,[foldername,'_maskVein_New.png']),'png') ;
imwrite(mat2gray(single(mask_vessel)),fullfile(ToolBox.PW_path_png,[foldername,'_maskVessel_New.png']),'png') ;
imwrite(mat2gray(single(maskBackground)),fullfile(ToolBox.PW_path_png,[foldername,'_maskBackground_New.png']),'png') ;
%vesselMap = uint8( cat( 3, uint8(meanIm)+ uint8(mask_artery)*255, uint8(meanIm) , uint8(meanIm) + uint8(mask_vein)*255 ));
imwrite(VesselImageRGB,fullfile(ToolBox.PW_path_png,[foldername,'_vesselMap.png']),'png') ;
imwrite(mat2gray(maskCRA),fullfile(ToolBox.PW_path_png,[foldername,'_maskCRA.png']),'png') ;
imwrite(mat2gray(maskCRA),fullfile(ToolBox.PW_path_png,[foldername,'_maskCRV.png']),'png') ;

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

w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_RG_video_vessel.avi')));
tmp = mat2gray(RG_video_vessel);
open(w)
for j = 1:size(RG_video_vessel,3)
    writeVideo(w,tmp(:,:,j)) ;  
end
close(w);


end
