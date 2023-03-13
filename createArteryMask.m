function [mask_artery_retina_choroid,mask_artery] = createArteryMask(video)
%FIXME : check the link with VesselMask (order between flat_field_correction and imbinarize)
%FIXME : add threshold parameter

mask_artery = squeeze(std(video, 0, 3));

figure(307)
imagesc(mask_artery)
colormap gray

mask_artery = imbinarize(mat2gray(mask_artery), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);

% C = {};
% C(1:2,1) = {'artery_mask_1stPass'; mask};

% figure(202)
% imagesc(mask);
% title('segmented arteries');
% axis off
% axis equal
% colormap gray

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
pulse_init_3d = zeros(size(video));
for mm = 1:size(video, 1)
    for pp = 1:size(video, 2)
        pulse_init_3d(mm,pp,:) = pulse_init;
    end
end
% compute local-to-average pulse wave zero-lag correlation 
C = C .* pulse_init_3d;
C = squeeze(mean(C, 3));
% find max. values of the correlation
mask_artery_retina_choroid = C > (max(C(:)) * 0.2);
mask_artery = mask_artery_retina_choroid;

figure(204)
imagesc(mask_artery_retina_choroid);

mask_artery = magicwand(mask_artery, 0.2, 8, 2);