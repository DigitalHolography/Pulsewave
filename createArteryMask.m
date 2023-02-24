function artery_mask = createArteryMask(video)
%FIXME : check the link with VesselMask (order between flat_field_correction and imbinarize)
%FIXME : add threshold parameter
mask = squeeze(std(video, 0, 3));
mask = imbinarize(mat2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);

figure(203)
imagesc(mask);
title('segmented arteries');
axis off
axis equal
colormap gray

imwrite(mat2gray(single(mask)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArtery_0.png')),'png') ;

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

figure(204)
imagesc(artery_mask);
title('segmented arteries');
axis off
axis equal
colormap gray

imwrite(mat2gray(single(artery_mask)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArtery_1.png')),'png') ;

%Magicwand to delete choroid signal
artery_mask = magicwand(artery_mask, 0.4, 8);

figure(205)
imagesc(artery_mask);
title('segmented arteries');
axis off
axis equal
colormap gray

imwrite(mat2gray(single(artery_mask)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArtery_2.png')),'png') ;

end