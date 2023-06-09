
function mask = createVesselMask(video)
mask = squeeze(mean(video, 3));


sigma = 3;
beta = 0.8;
mask= vesselness_filter(mask, sigma, beta);
mask = mask > (max(mask,[],'all') * 0.2);
mask = magicwand(mask, 0.2, 8, 2);


% figure(307)
% imagesc(mask)
% colormap gray
% 


end