function artery_mask = createArteryMask(video)
%FIXME : add threshold parameter
mask = squeeze(std(video, 0, 3));
mask = imbinarize(mat2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
pulse = squeeze(mean(video .* mask, [1 2]));
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
artery_mask = C > (max(C(:)) * 0.2);
end