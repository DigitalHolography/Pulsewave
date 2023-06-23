function mask = createVesselMask(video,path)

meanIm = squeeze(mean(video, 3));
PW_params = Parameters(path);



mask= vesselness_filter(meanIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);



mask = mask > (max(mask,[],'all')* PW_params.arteryMask_vesselness_bin_threshold);


 se = strel('disk',2);
% mask = imerode(mask,se);
%mask = magicwand(mask,meanIm, 0.2, 8, 3);
mask = imerode(mask,se);
mask = magicwand(mask,meanIm, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);
se = strel('disk',2);
mask = imdilate(mask,se);

figure(310)
imagesc(mask)
colormap gray

%%
% Imref = squeeze(mean(video, 3));
% 
% mask = flat_field_correction(Imref, ceil(size(Imref,1)*0.07), 0.25);
% 
% % mask = log10(abs(mask));
% %FIXME : add threshold parameter
% mask = imbinarize(mat2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.35);
% mask= magicwand(mask,Imref,0.2, 8, 2);
% figure(3030)
% imagesc(mask)
% colormap gray


end