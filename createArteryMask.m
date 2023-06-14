function [mask_artery_retina_choroid,mask_artery] = createArteryMask(video,path)
%FIXME : check the link with VesselMask (order between flat_field_correction and imbinarize)
%FIXME : add threshold parameter
PW_params = Parameters(path);

mask_artery = squeeze(mean(video, 3));


mask_artery= vesselness_filter(mask_artery, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);
mask_artery = mask_artery > (max(mask_artery,[],'all') * PW_params.arteryMask_bin_threshold);
mask_artery = magicwand(mask_artery, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);


%mask_artery = imbinarize(mat2gray(mask_artery), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);

figure(307)
imagesc(mask_artery)
colormap gray

% for ii = 1:size(video,3)
%     video(:,:,ii) = flat_field_correction(squeeze(video(:,:,ii)), ceil(size(video,1)*0.07), 0.25);
% end
% I  = std(video .* mask_artery, 0,3);
% 
% figure(3079)
% imagesc(I)



pulse = squeeze(mean(video .* mask_artery, [1 2]));
% zero-mean initial average pulse estimate : pulse_init
pulse_init = pulse - mean(pulse, "all");
C = video;
% image flattening
for ii = 1:size(video,3)
    video(:,:,ii) = flat_field_correction(squeeze(video(:,:,ii)), ceil(size(video,1)*0.07), 0.25);
end
% zero-mean local pulse : C
C(:,:,:) = video(:,:,:) - squeeze(mean(video, 3));
C = C.* mask_artery;
pulse_init_3d = zeros(size(video));
for mm = 1:size(video, 1)
    for pp = 1:size(video, 2)
        pulse_init_3d(mm,pp,:) = pulse_init;
    end
end
% compute local-to-average pulse wave zero-lag correlation 
C = C .* pulse_init_3d;
C = squeeze(mean(C, 3));
figure(204)
imagesc(C);

% find max. values of the correlation
mask_artery_retina_choroid = C > (max(C(:)) * PW_params.arteryMask_choroidCorrThreshold);
mask_artery = mask_artery_retina_choroid;

figure(204)
imagesc(mask_artery_retina_choroid);

mask_artery = magicwand(mask_artery, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_artery);


figure(204)
imagesc(mask_artery);

list_fig_close = [307,204];
for ii=1:length(list_fig_close)
    close(list_fig_close(ii));
end

end