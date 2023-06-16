
function mask = createVesselMask(video,path)
meanIm = squeeze(mean(video, 3));
PW_params = Parameters(path);



mask= vesselness_filter(meanIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
% figure(306)
% imagesc(mask)
% colormap gray
% 
%mask= vesselness_filter(meanIm ,3, 0.8);
%mask= vesselness_filter(mask, 2, 0.8);
figure(307)
imagesc(mask)
colormap gray

%mask_skel= bwskel(logical(mask),"MinBranchLength",500);
%mask= vesselness_filter(mask, 2, 0.8);
% figure(3077)
% imagesc(mask_skel)
% colormap gray
% 
% se = strel('disk',100);
% mask = imbothat(mask,se);
% 
% figure(308)
% imagesc(mask)
% colormap gray



%mask = mask > (max(mask,[],'all') * 0.2);
mask = mask > (max(mask,[],'all')* PW_params.arteryMask_bin_threshold);

% figure(3088)
% imagesc(mask)
% colormap gray


% se = strel('disk',2);
% mask = imerode(mask,se);
%mask = magicwand(mask,meanIm, 0.2, 8, 3);
mask = magicwand(mask,meanIm, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);
se = strel('disk',2);
mask = imdilate(mask,se);

figure(310)
imagesc(mask)
colormap gray

end