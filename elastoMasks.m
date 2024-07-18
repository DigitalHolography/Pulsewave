%% charge video
V = VideoReader('X:\230705_BRZ0546\230705_BRZ0546_OD1_1_0\avi\230705_BRZ0546_OD1_1_0_M0.avi');
videoM0 = zeros(V.Height, V.Width, V.NumFrames);

for n = 1:V.NumFrames
    videoM0(:, :, n) = rgb2gray(read(V, n));
end

%% loading parameters and compute useful variables
arteryMask_vesselness_sigma = 2;
arteryMask_vesselness_beta = 0.9;
arteryMask_ArteryCorrThreshold = 0.05;
masks_radius = 0.15;
masks_minSize = 150;
% masks_cleaningCoroid = false;
masks_showIntermediateFigures = true;
RG_veinConditionThreshold = 0.45;
RG_ArteryConditionThreshold = 3;

[N, M, L] = size(videoM0);
% tmpvideo = zeros(N*2-1,M*2-1,L);
% for i = 1:L
% tmpvideo(:,:,i) = interp2(videoM0(:,:,i),1);
% end
% videoM0 = tmpvideo;
% [N,M,L] = size(videoM0);

for pp = 1:L
    videoM0(:, :, pp) = flat_field_correction(squeeze(videoM0(:, :, pp)), 0.07 * (M + N) / 2, 0);

end

meanIm = squeeze(mean(videoM0, 3));
blurred_mask = imgaussfilt(double(meanIm), 0.05 * size(meanIm, 1), 'Padding', 0);
[y_barycentre, x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));
[x, y] = meshgrid(1:M, 1:N);
cercle_mask = sqrt((x - x_barycentre) .^ 2 + (y - y_barycentre) .^ 2) <= masks_radius * (M + N) / 2;

videoM0_zero = videoM0 - meanIm;

% compute vesselness response
vesselnessIm = vesselness_filter(meanIm, arteryMask_vesselness_sigma, arteryMask_vesselness_beta);

%%  Compute first correlation to find arteries

% compute pulse in 3 dimentions for correlation in all vessels
pulse = squeeze(mean(videoM0 .* (vesselnessIm > 0), [1 2]));
pulse_init = pulse - mean(pulse, "all");
pulse_init_3d = zeros(N, M, L);

for nn = 1:N

    for mm = 1:M
        pulse_init_3d(nn, mm, :) = pulse_init;
    end

end

% compute local-to-average pulse wave zero-lag correlation
correlationMatrix_artery = squeeze(mean((videoM0_zero .* pulse_init_3d), 3)) .* (vesselnessIm > 0);

% Create first artery mask
firstMaskArtery = (correlationMatrix_artery > 1.5 * mean2(correlationMatrix_artery(correlationMatrix_artery > 0)));

% if masks_cleaningCoroid
%     firstMaskArtery = firstMaskArtery & bwareafilt(firstMaskArtery | cercle_mask, 1, 4);
% end

firstMaskArtery = bwareaopen(firstMaskArtery, masks_minSize);

if masks_showIntermediateFigures
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
pulse_init_3d = zeros(N, M, L);

for nn = 1:N

    for mm = 1:M
        pulse_init_3d(nn, mm, :) = pulse_init;
    end

end

% compute local-to-average pulse wave zero-lag correlation
CorrelationMatrix = squeeze(mean((videoM0_zero .* pulse_init_3d), 3));

% compute mean correction correlation to find coroid maximums
meanVideoM0 = zeros(N, M, L);

for ll = 1:L
    meanVideoM0(:, :, L) = meanIm;
end

% Create correlation matrix to segment vein and arteries
correlationMatrix_artery = CorrelationMatrix ./ max(CorrelationMatrix, [], 'all');
correlationMatrix_vein = CorrelationMatrix ./ min(CorrelationMatrix, [], 'all');
correlationMatrix_vein(correlationMatrix_vein <- 1) = -1;

if masks_showIntermediateFigures
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

% Create seed artery mask and contidion artery mask for region growing
correlationMatrix_arteryCondition = correlationMatrix_artery .* rescale(vesselnessIm);

level_artery = graythresh(correlationMatrix_arteryCondition);
seeds_artery = imbinarize(correlationMatrix_arteryCondition, level_artery ./ 2) | (vesselness_artery > 0.5 * max(vesselness_artery, [], 'all'));
seeds_artery = bwareaopen(seeds_artery, masks_minSize);

condition_artery = imbinarize(correlationMatrix_arteryCondition, level_artery ./ RG_ArteryConditionThreshold);

% Create seed vein mask and contidion vein mask for region growing
correlationMatrix_veinCondition = correlationMatrix_vein .* rescale(vesselnessIm);

level_vein = graythresh(correlationMatrix_veinCondition);
seeds_vein = imbinarize(correlationMatrix_veinCondition, level_vein) | (vesselness_vein > 0.5 * max(vesselness_vein, [], 'all'));
seeds_vein = bwareaopen(seeds_vein, masks_minSize);

condition_vein = vesselness_vein > RG_veinConditionThreshold * mean2(vesselness_vein(seeds_vein));

% Cleaning condition
condition_artery = bwareaopen(condition_artery, masks_minSize);
condition_vein = condition_vein & ~condition_artery;
condition_vein = bwareaopen(condition_vein, masks_minSize);

% Cleaning seeds
seeds_artery = seeds_artery & condition_artery;
seeds_vein = seeds_vein & condition_vein;

% Compute region growing to segment vein and arteries
%[mask_artery,RG_video_artery] = region_growing_for_vessel(vesselness_artery, seeds_artery, condition_artery,path);
%[mask_vein,RG_video_vein]  = region_growing_for_vessel(vesselness_vein, seeds_vein, condition_vein, path);
[mask_vessel, RG_video_vessel] = region_growing_for_vessel(vesselnessIm, seeds_artery | seeds_vein, condition_vein | condition_artery, path);

mask_artery = mask_vessel .* correlationMatrix_artery;
mask_artery = mask_artery > arteryMask_ArteryCorrThreshold;

% mask_vein = mask_vessel & ~mask_artery;
mask_vein = mask_vessel .* correlationMatrix_vein;
mask_vein = mask_vein > 0;

if masks_showIntermediateFigures
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

% clear vesselness_vein vesselness_artery floor_vein floor_artery level_vein level_artery seeds_vein seeds_artery ;
%% Cleaning coroid from masks

mask_artery = bwareaopen(mask_artery, masks_minSize);
mask_artery = imdilate(mask_artery, strel('disk', 2));
mask_artery = imclose(mask_artery, strel('disk', 5));

mask_vein = bwareaopen(mask_vein, masks_minSize);
mask_vein = imdilate(mask_vein, strel('disk', 2));
mask_vein = imclose(mask_vein, strel('disk', 5));

mask_vessel = mask_artery | mask_vein;

if masks_showIntermediateFigures
    figure(40), imagesc(uint8(cat(3, uint8(meanIm) + uint8(mask_artery) * 255, uint8(meanIm), uint8(meanIm) + uint8(mask_vein) * 255)));
    title('Artery/Vein region growing segmentation');
end
