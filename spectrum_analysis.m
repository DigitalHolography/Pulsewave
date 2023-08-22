function [] = spectrum_analysis(maskArtery,maskBackground ,SH_cube, ToolBox)




f1 = 1;
f2 = 15;
f3 = 30;
fs = ToolBox.fs/2 ;
cubeFrameLength = size(SH_cube,4);
batch_size = size(SH_cube, 3);
gw = 3 ;
SH_ColorVideoRGB = zeros(size(SH_cube,1),size(SH_cube,2),3,size(SH_cube,4));

%% integration intervals
low_n1 = round(f1 * batch_size / fs) + 1;
low_n2 = round(f2 * batch_size / fs);
high_n1 = low_n2 + 1;
high_n2 = round(f3 * batch_size / fs);

MeanFreqLow = mean(SH_cube(:, :, low_n1:low_n2,:),4);
MeanFreqHigh= mean(SH_cube(:, :, high_n1:high_n2,:),4);
MeanImLow = mat2gray(squeeze(sum(abs(MeanFreqLow), 3)));
MeanImHigh = mat2gray(squeeze(sum(abs(MeanFreqHigh), 3)));

multiband_img = cat(3, MeanImLow, MeanImHigh);
DCR_imgs = decorrstretch(multiband_img, 'tol', [0.002 0.999]);
%imgAVG = imfuse(multiband_img(:,:,2), multiband_img(:,:,1), 'ColorChannels', 'red-cyan');
imgAVG = imfuse(DCR_imgs(:,:,2), DCR_imgs(:,:,1), 'ColorChannels', 'red-cyan');
low_high = stretchlim(imgAVG, [0, 1]);
gamma_composite = 0.8;
imgAVG = imadjust(imgAVG, low_high, low_high, gamma_composite);
imgAVG = imsharpen(imgAVG, 'Radius', 10, 'Amount', 0.6);


figure(15)
imshow(imgAVG)





for ii = 1:cubeFrameLength

    %% integration
    freq_low = squeeze(sum(abs(SH_cube(:, :, low_n1:low_n2,ii)), 3));
    freq_high = squeeze(sum(abs(SH_cube(:, :, high_n1:high_n2,ii)), 3));


%     %% normalization
%     freq_low = freq_low ./ imgaussfilt(freq_low, gw);
%     freq_high = freq_high./ imgaussfilt(freq_high, gw);


    M_freq_low = squeeze(freq_low);
    M_freq_high = squeeze(freq_high);

    avg_M0_low = mean(M_freq_low, 3);
    avg_M0_high = mean(M_freq_high, 3);

    avg_M0_low = mat2gray(avg_M0_low);
    avg_M0_high = mat2gray(avg_M0_high);

    %% composite generation
    multiband_img = cat(3, avg_M0_low, avg_M0_high);
    DCR_imgs = decorrstretch(multiband_img, 'tol', [0.02 0.998]);
    img = imfuse(DCR_imgs(:,:,2), DCR_imgs(:,:,1), 'ColorChannels', 'red-cyan');
    low_high = stretchlim(img, [0, 1]);
    gamma_composite = 0.8;
    img = imadjust(img, low_high, low_high, gamma_composite);
    img = imsharpen(img, 'Radius', 10, 'Amount', 0.6);
%     tmp = zeros(2 * size(img, 1) -1, 2 * size(img, 2) - 1, size(img, 3));
%     for mm = 1:size(img, 3)
%         tmp(:,:,mm) = interp2(single(img(:,:,mm)), 1);
%     end



   figure(15)
   imshow(img)

   if mod(ii,20)==0
        disp(['Frame : ',num2str(ii) , '/',cubeFrameLength])
   end

   SH_ColorVideoRGB(:,:,:,ii) = mat2gray(img);
   


end
%% save video

% avi
w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_SH_ColorVideo'))) ;
open(w)
for jj = 1:size(SH_ColorVideoRGB,4)
    writeVideo(w,squeeze(SH_ColorVideoRGB(:,:,:,jj))) ;
end
close(w);

imwrite(imgAVG,fullfile(ToolBox.PW_path_png,[foldername,'_ColorDoppler.png']),'png') ;



end