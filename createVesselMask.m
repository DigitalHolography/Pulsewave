
function mask = createVesselMask(video,path)
mask = squeeze(mean(video, 3));
PW_params = Parameters(path);



 mask= vesselness_filter(mask, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
% figure(306)
% imagesc(mask)
% colormap gray
% 
% mask= vesselness_filter(mask, 3, 0.8);
% figure(307)
% imagesc(mask)
% colormap gray


%se = strel('disk',1);
% mask = imtophat(mask,se);
%mask = imopen(mask,se);
% mask = imclose(mask,se);
%mask = mask > (max(mask,[],'all') * 0.2);
mask = mask > (max(mask,[],'all')* PW_params.arteryMask_bin_threshold);
% 
% figure(308)
% imagesc(mask)
% colormap gray

se = strel('disk',1);
mask = imerode(mask,se);

%mask = magicwand(mask, 0.2, 8, 6);
mask = magicwand(mask, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);
se = strel('disk',2);
mask = imdilate(mask,se);

% figure(310)
% imagesc(mask)
% colormap gray

end