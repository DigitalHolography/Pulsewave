function [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV] = createMasks(videoM0,videoM1M0, path, ToolBox)

PW_params = Parameters(path);
%video(:,:,:) = video(:,:,:) ./ mean(video(:,:,:), [1 2]);
stdIm = std(videoM0,0,3);
meanIm = squeeze(mean(videoM0, 3));
meanM1M0 = squeeze(mean(videoM1M0, 3));


%% Creat Vessel Mask

maskVessel = vesselness_filter(meanIm, PW_params.arteryMask_vesselness_sigma, PW_params.arteryMask_vesselness_beta);


maskVessel = maskVessel > (max(maskVessel,[],'all')* PW_params.arteryMask_vesselness_bin_threshold);

maskCorr = maskVessel;

se = strel('disk',1);
maskVessel = imerode(maskVessel,se);
maskVessel = magicwand(maskVessel,meanIm, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_vessels);
se = strel('disk',2);
maskVessel = imdilate(maskVessel,se);



%% Creat Artery Mask

maskArtery = imbinarize(mat2gray(stdIm), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);

pulse = squeeze(mean(videoM0 .* maskArtery, [1 2]));
% zero-mean initial average pulse estimate : pulse_init
pulse_init = pulse - mean(pulse, "all");
C = videoM0;

% image flattening
for ii = 1:size(videoM0,3)
    videoM0(:,:,ii) = flat_field_correction(squeeze(videoM0(:,:,ii)), ceil(size(videoM0,1)*0.07), 0.25);
end

% zero-mean local pulse : C
C(:,:,:) = videoM0(:,:,:) - meanIm;
C = C.* maskCorr  ;
pulse_init_3d = zeros(size(videoM0));
for mm = 1:size(videoM0, 1)
    for pp = 1:size(videoM0, 2)
        pulse_init_3d(mm,pp,:) = pulse_init;
    end
end
% compute local-to-average pulse wave zero-lag correlation 
C = C .* pulse_init_3d;
C = squeeze(mean(C, 3));

maskArtery = C > (max(C(:)) * PW_params.arteryMask_ArteryCorrThreshold);
maskArtery = magicwand(maskArtery,meanIm, 0.2, 8, PW_params.arteryMask_magicwand_nb_of_area_artery);


%% Create Vein Mask 

maskVein = double(maskVessel) - double(maskArtery); 
maskVein = maskVein > 0; 
se = strel('disk',3);
maskVein = imerode(maskVein,se);
maskVein = imdilate(maskVein,se);
%% Create CRA and CRV Mask

stdM1M0 = std2(meanM1M0);
maskCRA = meanM1M0>(3.0*stdM1M0);
maskCRV = meanM1M0<(-3.0*stdM1M0);

%% Creat Background Mask

maskBackground = not(maskVessel);

%% Saving masks as PNG 
foldername = ToolBox.main_foldername;
imwrite(mat2gray(single(maskArtery)),fullfile(ToolBox.PW_path_png,[foldername,'_maskArtery.png']),'png') ;
imwrite(mat2gray(single(maskVein)),fullfile(ToolBox.PW_path_png,[foldername,'_maskVein.png']),'png') ;
imwrite(mat2gray(single(maskVessel)),fullfile(ToolBox.PW_path_png,[foldername,'_maskVessel.png']),'png') ;
imwrite(mat2gray(single(maskBackground)),fullfile(ToolBox.PW_path_png,[foldername,'_maskBackground.png']),'png') ;



end
