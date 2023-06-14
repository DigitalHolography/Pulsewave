
function mask = createVesselMask(video,path)
mask = squeeze(mean(video, 3));
PW_params = Parameters(path);


mask= vesselness_filter(mask, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
mask = mask > (max(mask,[],'all') * PW_params.arteryMask_bin_threshold);
mask = magicwand(mask, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);


% figure(307)
% imagesc(mask)
% colormap gray
% 


end