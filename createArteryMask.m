function mask_artery = createArteryMask(video)
%FIXME : check the link with VesselMask (order between flat_field_correction and imbinarize)
%FIXME : add threshold parameter
mask = squeeze(std(video, 0, 3));
mask = imbinarize(mat2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);

% C = {};
% C(1:2,1) = {'artery_mask_1stPass'; mask};

figure(202)
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
mask_artery = C > (max(C(:)) * 0.2);


%% FIXME interpolation
tmp = interp2(squeeze(video(:,:,1)),1);
video_interp = zeros(size(tmp,1),size(tmp,2),size(video,3)); 
for jj=1:size(video,3)
    video_interp(:,:,jj) = interp2(video(:,:,jj),1);
end
%FIXME C_interp = zeros(size(video_interp));
C_interp = video_interp;
% image flattening
for ii = 1:size(video_interp,3)
    video_interp(:,:,ii) = flat_field_correction(squeeze(video_interp(:,:,ii)), ceil(size(video_interp,1)*0.07), 0.25);
end
% zero-mean local pulse : C
C_interp(:,:,:) = video_interp(:,:,:) - squeeze(mean(video_interp, 3));
pulse_init_3d_interp = zeros(size(video_interp));
for mm = 1:size(video_interp, 1)
    for pp = 1:size(video_interp, 2)
        pulse_init_3d_interp(mm,pp,:) = pulse_init;
    end
end
% compute local-to-average pulse wave zero-lag correlation 
C_interp = C_interp .* pulse_init_3d_interp;
C_interp = squeeze(mean(C_interp, 3));
% find max. values of the correlation
mask_artery_hd = C_interp > (max(C_interp(:)) * 0.2);

figure(203)
imagesc(mask_artery_hd);
title('segmented arteries HD');
axis off
axis equal
colormap gray

figure(204)
imagesc(mask_artery);
title('segmented arteries');
axis off
axis equal
colormap gray

imwrite(mat2gray(single(mask_artery)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArtery_1.png')),'png') ;

%Magicwand to delete choroid signal
mask_artery = magicwand(mask_artery, 0.4, 8);
mask_artery_hd = magicwand(mask_artery_hd, 0.4, 8);

% C(1:2,3) = {'artery_mask_3rdPass'; artery_mask};

figure(205)
imagesc(mask_artery_hd);
title('segmented arteries HD (after wand)');
axis off
axis equal
colormap gray

figure(206)
imagesc(mask_artery);
title('segmented arteries (after wand)');
axis off
axis equal
colormap gray


imwrite(mat2gray(single(mask_artery)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArtery_2.png')),'png') ;

imwrite(mat2gray(single(mask_artery_hd)),fullfile("C:\Users\Rakushka\Downloads",strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArteryHD_2.png')),'png') ;


end