function pulseAnalysis(Ninterp, fullVideo, fullVideoM1M0, one_pulse_video, one_cycle_dir, filename, sys_index_list, maskArtery)
%
% one_pulse_video :
% one_cycle_dir :
% filename :
% sys_index_list :

%% FIXME : const;
borderAmount = 0; % for flatFieldCorrection() and createBorderMask
dataCube = fullVideo;
dataCubeM1M0 = fullVideoM1M0;

%RMS freq to velocity
theta = 2.5/20;
% theta = 38 * pi / 180;
% theta = 0.03;

opticalIndex = 1.35;
lambda = 852e-9;
scalingFactorVelocity = 1000 * 1000 * lambda / (3 *opticalIndex * theta); % 1000 for kHz -> Hz and 1000 for m -> mm
scalingFactorVelocity2 = 1000 * 1000 * lambda / opticalIndex * (3/theta)^(1/2); % 1000 for kHz -> Hz and 1000 for m -> mm
%scalingFactorVelocityCRA2  = 1000 * 1000 * lambda / ((pi*(theta+sin(2*theta)/2))^(1/2)); % 1000 for kHz -> Hz and 1000 for m -> mm

% for robust rendering : 
% 1-flat-field correction, 2-background substraction
for pp = 1:size(dataCube,3)
    dataCube(:,:,pp) = flat_field_correction(squeeze(dataCube(:,:,pp)), 0.07*size(dataCube,1), 0.25);
end
%

%%

% FIXME maskArtery : best computed from reference video
% maskArtery = createArteryMask(fullVideo);
figure(5)
imagesc(maskArtery);
title('segmented arteries');
axis off
axis equal
colormap gray

%
maskVessel = createVesselMask(fullVideo);
figure(4)
imagesc(maskVessel);
title('segmented vessels');
axis off
axis equal
colormap gray

%
fullArterialPulse = fullVideo .* maskArtery;
fullArterialPulse = squeeze(sum(fullArterialPulse, [1 2]))/nnz(maskArtery);

% maskBorder = createBorderMask(maskVessel,borderAmount);
% maskBackground = not(maskVessel) .* maskBorder;
maskBackground = not(maskVessel);

%
fullBackgroundSignal = fullVideo .* maskBackground;
fullBackgroundSignal = squeeze(sum(fullBackgroundSignal, [1 2]))/nnz(maskBackground);

[cache_exists, batch_stride, Fs] = getTimelineParamsFromCache(one_cycle_dir);
strXlabel = createXlabelTime(cache_exists);
if cache_exists
    fullTime = linspace(0,size(fullVideo,3)*batch_stride/Fs,size(fullVideo,3));
else
    fullTime = 1:size(fullVideo,3);
end

figure(2)
plot(fullTime,fullArterialPulse,'-k', fullTime,fullBackgroundSignal,':k','LineWidth',2) ;
title('arterial pulse waveform and background signal'); % averaged outside of segmented vessels
legend('arterial pulse','background') ;
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

% Renormalization to avoid bias
% renormalize avg. background to avg. as local pulse
avgFullBackgroundSignal = mean(fullBackgroundSignal(:));
avgFullArterialPulse = mean(fullArterialPulse(:));
% fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal * avgFullArterialPulse / avgFullBackgroundSignal;
fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal;
fullArterialPulse = fullArterialPulseMinusBackground;

% same with dataCube
% dataCube = dataCube - fullBackgroundSignal * avgFullArterialPulse / avgFullBackgroundSignal;
% alternatively : 
% dataCube = dataCube - fullBackgroundSignal

% le plus simple fonctionne bien : soustraire le bkg.
A = ones(size(dataCube));
% B = ones(size(dataCube));
for pp = 1:size(dataCube,3)
%     A(:,:,pp) = A(:,:,pp) * fullBackgroundSignal(pp) * avgFullArterialPulse / avgFullBackgroundSignal;
      A(:,:,pp) = A(:,:,pp) * fullBackgroundSignal(pp);
%     B(:,:,pp) = B(:,:,pp) * fullBackgroundSignal(pp);    
end
dataCube = dataCube - A ;
% dataCube = dataCube - A .* maskVessel;
% dataCube = dataCube - B .* maskBackground;


% remove outliers
% 1st pass
disp('remove outliers... 1st pass.');
[idxOutPw,fullArterialPulseRmOut] = discardPulseWaveOutliers(fullArterialPulse,3);
[idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundSignal,3);
dataReliabilityIndex1 = ceil(100*(1-0.5*(length(idxOutPw)/length(fullArterialPulse) + length(idxOutBkg)/length(fullBackgroundSignal))));
disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);

% smooth trendline data by iterative local linear regression.
% fullArterialPulseClean = smoothdata(fullArterialPulse,'rloess');
fullArterialPulseClean = smoothdata(fullArterialPulseRmOut,'lowess');
fullBackgroundSignalClean = smoothdata(fullBackgroundSignalRmOut,'lowess');

figure(8)
plot(fullTime,fullArterialPulse,':k', ...
    fullTime,fullArterialPulseClean,'-k', ...
    'LineWidth',2) ;
title('arterial pulse minus background vs. filtered pulse');
legend('<p(t)> - <b(t)>','local linear regression');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

fullArterialPulseDerivative = diff(fullArterialPulse);%
fullArterialPulseCleanDerivative = diff(fullArterialPulseClean);%
% dataCubeDiff = diff(dataCube,1,3);
% fullArterialPulseDerivative = dataCubeDiff .* maskArtery;
% fullArterialPulseDerivative = squeeze(sum(fullArterialPulseDerivative, [1 2]))/nnz(maskArtery);

figure(6)
plot( ...
    fullTime(1:length(fullArterialPulseDerivative)),fullArterialPulseDerivative,':k', ...
    fullTime(1:length(fullArterialPulseCleanDerivative)),fullArterialPulseCleanDerivative,'-k', ...
    'LineWidth',2) ;
text(fullTime(sys_index_list)-0.3,fullArterialPulseCleanDerivative(sys_index_list)+0.03,num2str((1:numel(sys_index_list))'))
title('derivative of the arterial pulse waveform');
legend('\delta <p(t)> - <b(t)>','from smoothed data');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

%create a function findSystoleIndexFromTrendline()
[~, sys_index_list] = findpeaks(fullArterialPulseCleanDerivative, ...
    1:length(fullArterialPulseCleanDerivative), ...
    'MinPeakHeight', max(fullArterialPulseCleanDerivative) * 0.7);

%%
% now cleanup dataCube to create_one_cycle()
% strategy : use average pulse profiles to detect and
% find  noisy frames @ >3 std from zero-mean
% replace them with cleaned data
noise = sqrt(abs(abs(fullArterialPulseMinusBackground).^2 - abs(fullArterialPulseClean).^2));
idxOutNoise = find(noise>4*std(noise));

%
dataReliabilityIndex2 = ceil(100*(1-(length(idxOutNoise)/length(fullArterialPulseClean) )));
disp(['data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);

% remove bad frames and tackle renormalization issue later
% if ~isempty(idxOutNoise)    
%     for pp = 1:length(idxOutNoise)
% % regularized frame avg
%         avgFrameBackgroundSignalClean = squeeze(sum(fullBackgroundSignalClean(idxOutNoise(pp)) * maskBackground,[1 2]) / nnz(maskVessel));
%         avgFrameArterialPulseClean = squeeze(sum(fullArterialPulseClean(idxOutNoise(pp)) * maskVessel,[1 2]) / nnz(maskVessel));
% %         
%         dataCube(:,:,idxOutNoise(pp)) = ...
%         fullArterialPulseClean(idxOutNoise(pp)) * maskVessel + ...
%         fullBackgroundSignalClean(idxOutNoise(pp)) * maskBackground;
%     end
% end

% figure(23)
% imagesc(squeeze(mean(fullVideo,3))) ;
% colormap gray
% title('raw RMS frequency map');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image
% c = colorbar('southoutside');
% c.Label.String = 'RMS Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% 
% dMap = squeeze(mean(fullVideo,3));
% dMap = flat_field_correction(dMap, ceil(0.07*size(dMap,1)), .33);
% 
% figure(24)
% imagesc(dMap) ;
% colormap gray
% title('flattened RMS frequency map');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image
% c = colorbar('southoutside');
% c.Label.String = 'RMS Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% range0(1:2) = clim;
% % range0 = [1.2*range0(1),0.8*range0(2)];
% figure(23)
% clim(range0);

% FIXME : compute true regularied cube by replacing bad frames
fullArterialPulseRegularized = squeeze(sum(dataCube .* maskArtery, [1 2])) / nnz(maskArtery);

% NanInFullArterialPulseRegularized = isnan(fullArterialPulseRegularized)

figure(22)
% %     
plot( ...
    fullTime(1:length(fullArterialPulseRegularized)),fullArterialPulseRegularized,'-k', ...
    fullTime(1:length(fullArterialPulse)),fullArterialPulse,':k', ...
    'LineWidth',2) ;
yline(0,':',{''},LineWidth=2) ;
title('Regularized pulse wave signal averaged in retinal arteries');
legend('Regularized arterial pulse','<p(t)> - <b(t)>');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;
%%


figure(9)
plot( ...
    fullTime(1:length(fullArterialPulseClean)),fullArterialPulseClean,':k', ...
    fullTime(1:length(noise)),noise,'-k', ...
    'LineWidth',2) ;
yline(0,':',{''},LineWidth=2) ;
yline(std(noise),':',{'1 std'},LineWidth=2) ;
yline(2*std(noise),':',{'2 std'},LineWidth=2) ;
yline(3*std(noise),':',{'3 std'},LineWidth=2) ;
title('signal vs. noise');
legend('filtered arterial pulse','residual');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

% FIXME FIXME
% 1- select reliable pulses
% 2- get reliable bias measure

[maskCRA, ~, ~] = createCentralRetinalArteryVeinMask(dataCubeM1M0);

figure(502)
imagesc(maskArtery-maskCRA.*maskArtery)

[onePulseVideo2, selectedPulseIdx] = create_one_cycle(dataCube, maskArtery-maskCRA.*maskArtery, sys_index_list, Ninterp);
avgArterialPulse =  onePulseVideo2 .* (maskArtery-maskCRA.*maskArtery);
avgArterialPulse = squeeze(sum(avgArterialPulse, [1 2]))/nnz(maskArtery-maskCRA.*maskArtery);



% ATTN : Substract baseline signal from avgArterialPulse to onePulseVideo2
% minAvgArterialPulse = min(avgArterialPulse(:));
% onePulseVideo2 = onePulseVideo2 - minAvgArterialPulse;

% FIXME : use selectedPulseIdx in fig. 8 to highlight selected pulses

nb_frames = size(onePulseVideo2,3) ;
blur_time_sys = ceil(nb_frames/100);
blur_time_dia = ceil(nb_frames/100);

%
if cache_exists % .mat with cache from holowaves is present, timeline can be computed
    average_cycle_length = 0;
    nb_of_averaged_cycles = 0;
    if size(sys_index_list,2) == 1
        average_cycle_length = nb_frames ;
        nb_of_averaged_cycles = 1;
    else
        for ii = 2:size(sys_index_list,2)
            average_cycle_length = average_cycle_length + (sys_index_list(ii)-sys_index_list(ii-1)) ;
            nb_of_averaged_cycles = nb_of_averaged_cycles + 1;
        end
        average_cycle_length = average_cycle_length / (length(sys_index_list)-1);
    end
    T = linspace(0,batch_stride/Fs*average_cycle_length,nb_frames);
else % mat with cache from holowaves is not present, timeline cannot be computed
    T = 1:nb_frames;
end

avgArterialPulseVelocity = avgArterialPulse * scalingFactorVelocity;
avgArterialPulseVelocity2 = avgArterialPulse * scalingFactorVelocity2;

% avgArterialPulseVelocityCRA = avgArterialPulseCRA * scalingFactorVelocityCRA2;
% avgArterialPulseVelocityCRA2 = avgArterialPulseCRA * scalingFactorVelocityCRA2;

% figure(1)
% plot( ...
%     T(1:length(avgArterialPulse)),avgArterialPulseVelocity,'-k', ...
%     'LineWidth',2) ;
% title('average blood flow velocity estimate in retinal arteries');
% legend(' arterial pulse');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% ylabel('blood flow velocity (mm/s)');
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;

figure(111)
plot( ...
    T(1:length(avgArterialPulse)),avgArterialPulseVelocity2,'-k', ...
    'LineWidth',2) ;
title('average blood flow velocity 2 estimate in retinal arteries');
legend(' arterial pulse');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel('blood flow velocity (mm/s)');
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

% figure(301)
% plot( ...
%     T(1:length(avgArterialPulseCRA)),avgArterialPulseVelocityCRA,'-k', ...
%     'LineWidth',2) ;
% title('average blood flow velocity estimate in CRA');
% legend(' arterial pulse');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% ylabel('blood flow velocity (mm/s)');
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;

% figure(311)
% plot( ...
%     T(1:length(avgArterialPulseCRA)),avgArterialPulseVelocityCRA2,'-k', ...
%     'LineWidth',2) ;
% title('average blood flow velocity 2 estimate in CRA');
% legend(' arterial pulse');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% ylabel('blood flow velocity (mm/s)');
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;

% save average arterial pulse wave velocity to txt file
tmp = [T(1:length(avgArterialPulse))',avgArterialPulseVelocity];
%size(tmp)
fileID = fopen(fullfile(one_cycle_dir, strcat(filename,'_avgPulse.txt')),'w') ;
fprintf(fileID,'%f %f \r\n',tmp');
fclose(fileID);


disp('arterial resistivity...');
[ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma] = construct_resistivity_index(onePulseVideo2, maskArtery);
ARImap = ARImap.*maskArtery;

% export fig
figure(15)
imshow(ARImapRGB) ;
title(strcat('Arterial resistivity. avg. index value : ', num2str(ARI)));
axis image
axis off
set(gca,'LineWidth', 2);
fontsize(gca,12,"points") ;
c = colorbar('southoutside');
c.Label.String = 'arterial resistivity index';
c.Label.FontSize = 12;
% cmap = colormap(gray);
% cmap(:,2) = 0;
% cmap(:,3) = 0;
cmap = double(ones(size(256,3)));

% x = 1:256; 
% y = sigmoid(x,128,0.2);

for ii = 0 : 255
    cmap(ii+1, 1) = 1;
    cmap(ii+1, 2) = (double(1 - ii/255))^gamma;
    cmap(ii+1, 3) = (double(1 - ii/255))^gamma;
end

colormap(cmap);


% export RImap
imwrite(ARImapRGB,fullfile(one_cycle_dir,strcat(filename,'_ARI_map_img.png')),'png');
print('-f15','-dpng',fullfile(one_cycle_dir,strcat(filename,'_ARI_map.png')));

w = VideoWriter(fullfile(one_cycle_dir,strcat(filename,'_ARIvideoRGB.avi')));
open(w)
ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
end
close(w);

% figure(100)
% imagesc(ARImapRGB);
% figure(101)
% imshow(squeeze(ARIvideoRGB(:,:,:,1)));
% figure(102)
% imagesc(ARImap);

disp('done.');

%%
disp('arterial pulse wave analysis...');

%%
%zero-mean local-to-average arterial pulse cross-correlation
% reference video : average arterial pulse everywhere
% avgArterialPulse_3d = zeros(size(one_pulse_video));
% for mm = 1:size(one_pulse_video, 1)
%     for pp = 1:size(one_pulse_video, 2)
%         avgArterialPulse_3d(mm,pp,:) = avgArterialPulse;
%     end
% end
% [max_C_3, id_max] = zmXCorrVideos(one_pulse_video,avgArterialPulse_3d);
%%

% %     Time lags from cross-correlation between the average arterial pulse and local pulse waveforms
% figure(99);
% imagesc(id_max);
% title('local-to-average arterial pulse lag');
% colormap default
% colorbar;
% axis image
% axis off

% %     cross-correlation value between the average arterial pulse and local pulse waveforms
% figure(77)
% imagesc(log(abs(max_C_3)))
% title('local-to-average arterial pulse maximum cross-correlation value');
% colormap gray;
% colorbar ;
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image


% %% artery mask by lag-correlation with average pulse
% arteries = zeros(size(id_max));
% % arteries (or(id_max < 0.1*nb_frames, id_max > 0.9 * nb_frames)) = 1 ;
% arteries (or(id_max <= ceil(0.01*nb_frames), id_max >= floor(0.95 * nb_frames))) = 1 ;
% %     arteries = createArteryMask(one_pulse_video) ;
% figure(10)
% %denoising
% %arteries = medfilt2(arteries,[3 3]);
% imagesc(arteries) ;
% colormap gray
% title('segmented arteries by lag-correlation with average pulse');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image

% %% vein mask by lag-correlation with average pulse
% veins = zeros(size(id_max)) ;
% veins(and(id_max>=ceil(0.1*nb_frames),id_max<=floor(0.3*nb_frames))) = 1 ;
% %denoising
% %veins = medfilt2(veins,[3 3]);
% figure(11)
% imagesc(veins) ;
% colormap gray
% title('segmented veins by cross-correlation with average arterial pulse');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image
















%% calculation of average pulse in arteries
I_arteries = one_pulse_video .* maskArtery;
% avgArterialPulse = squeeze(sum(I_arteries, [1 2]))/nnz(maskArtery);
% avgArterialPulse = avgArterialPulse - min(avgArterialPulse) ;
% [~,idx] = min(avgArterialPulse,[],1) ;
% avgArterialPulse = circshift(avgArterialPulse,-idx);

% %% calculation of pulse wave in veins
% I_veins = one_pulse_video .* veins ;
% pulse_veins = squeeze(sum(I_veins, [1 2]))/nnz(veins);
% pulse_veins = pulse_veins - min(pulse_veins) ;
% pulse_veins = circshift(pulse_veins,-idx) ;
% 
% max_plot = max(max(avgArterialPulse(:)),max(pulse_veins(:))) ;
% avgArterialPulse = avgArterialPulse ./ max_plot ;
% pulse_veins = pulse_veins ./ max_plot ;

%% plot pulses in veins and arteries

% find peak systole index
sys_index_list_one_cycle = find_systole_index(onePulseVideo2);
if cache_exists % .mat with cache from holowaves is present, timeline can be computed
    [~,idx_sys] = max(avgArterialPulse) ;
    T_syst = batch_stride/Fs*average_cycle_length * (idx_sys-1) / nb_frames ;
else % no .mat present, hence no timeline, T == frames instead of time in s.
    [~,idx_sys] = max(avgArterialPulse) ;
    T_syst = idx_sys ;
    T = 1:nb_frames;
end

%% diastolic Doppler frequency heatmap : 10% of frames before minimum of diastole
heatmap_dia = squeeze(mean(one_pulse_video(:,:,floor(0.9*nb_frames):nb_frames),3));
% onePulseVideo2 : no background correction 
% heatmap_dia = squeeze(mean(onePulseVideo2(:,:,floor(0.9*nb_frames):nb_frames),3));
heatmap_dia = flat_field_correction(heatmap_dia, ceil(.07*size(heatmap_dia,1)), .33);
% figure(45)
% imagesc(heatmap_dia) ;
% colormap gray
% title('bottom diastole RMS frequency map');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% c = colorbar('southoutside');
% c.Label.String = 'RMS Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% axis off
% axis image
% range(1:2) = clim;

%% systolic Doppler frequency heatmap : 10% of frames around peak systole
a = max(ceil(idx_sys-0.05*nb_frames),1);
b = min(ceil(idx_sys+0.05*nb_frames),nb_frames);
heatmap_sys = squeeze(mean(one_pulse_video(:,:,a:b),3));
% onePulseVideo2 : no background correction 
% heatmap_sys = squeeze(mean(onePulseVideo2(:,:,a:b),3));
heatmap_sys = flat_field_correction(heatmap_sys, ceil(.07*size(heatmap_sys,1)), .33);
% figure(46)
% imagesc(heatmap_sys) ;
% colormap gray
% title('peak systole RMS frequency map');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% c = colorbar('southoutside');
% c.Label.String = 'RMS Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% axis off
% axis image
% range(3:4) = clim;
% % same color axis for systolic and diastolic Doppler heatmaps
% clim([min(range),max(range)]);
% figure(45)
% clim([min(range),max(range)]);


% %vein mask
% B_dia = squeeze(mean(heatmap_dia .* not(maskVessel), [1 2]));
% heatmap_dia = heatmap_dia - B_dia;
% A_dia = squeeze(mean(heatmap_dia .* maskArtery, [1 2]));
% heatmap_dia = heatmap_dia / A_dia;
% B_sys = squeeze(mean(heatmap_sys .* not(maskVessel), [1 2]));
% heatmap_sys = heatmap_sys - B_sys;
% A_sys = squeeze(mean(heatmap_sys .* maskArtery, [1 2]));
% heatmap_sys = heatmap_sys / A_sys;
% maskVein = heatmap_dia - heatmap_sys;
% maskVein(maskVein<0)=0;
% maskVein = imbinarize(mat2gray(maskVein), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);
% figure(47)
% imagesc(maskVein) ;
% colormap gray
% title('veins Doppler heatmap');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% c = colorbar('southoutside');
% c.Label.String = 'RMS Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% axis off
% axis image
% 
% figure(70)
% plot(T,avgArterialPulse,'k-',T,pulse_veins,'k--',LineWidth=2) ;
% xline(T_syst,':',{''},LineWidth=2) ;
% legend('Arteries','Veins') ;
% fontsize(gca,12,"points") ;
% xlabel('Time (s)','FontSize',14) ;
% ylabel('Normalized intensity','FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% title('average pulse waves in arteries and veins');
% axis tight

% FIXME : replace sys + dia blur by homogenous blur ? 
pulse_arteries_blurred_sys = movavgvar(avgArterialPulse(1:idx_sys), blur_time_sys);
diff_pulse_sys = diff(pulse_arteries_blurred_sys) ;
pulse_arteries_blurred_dia = movavgvar(avgArterialPulse(idx_sys:end), blur_time_dia);
diff_pulse_dia = diff(pulse_arteries_blurred_dia) ;
diff_avgPulse = diff(movavgvar(avgArterialPulse, blur_time_sys));
delta = max(pulse_arteries_blurred_dia(:))-min(pulse_arteries_blurred_dia(:));
thr = max(pulse_arteries_blurred_dia(:))-delta./exp(1);
idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < thr);
[max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
acc = max_diff_pulse/(T(idx_T_diff_max)-T(idx_T_diff_max-1));

figure(80)
plot(T,avgArterialPulse,'k.', ...
    T(1:idx_sys),pulse_arteries_blurred_sys(1:idx_sys),'k-', ...
    T(idx_sys:nb_frames),pulse_arteries_blurred_dia(1:(nb_frames-idx_sys+1)),'k-', ...
    LineWidth=2) ;
xline(T(idx_T_diff_max + 1),':',{},LineWidth=2)
text(T(idx_T_diff_max + 1),min(pulse_arteries_blurred_dia(:))+0.1,' (1)','FontSize',14);
xline(T_syst,':',{},LineWidth=2) ;
text(T_syst,min(pulse_arteries_blurred_dia(:))+0.1,'  (2)','FontSize',14);
xline(T(idx_sys + idx_list_threshold_dia(1)),':',{},LineWidth=2);
text(T(idx_sys + idx_list_threshold_dia(1)),min(pulse_arteries_blurred_dia(:))+0.1,'  (3)','FontSize',14);
%                 yline(1/exp(1),':',LineWidth=2);
legend('arterial pulse','smoothed line') ;
fontsize(gca,12,"points") ;
xlabel('Time (s)','FontSize',14) ;
ylabel('background-corrected RMS frequency (kHz)','FontSize',14) ;
axis tight;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
title('average background-corrected RMS frequency in retinal arteries');

% figure(90)
% plot(T(1:end-1), diff_avgPulse,'k-', LineWidth=2);
% x = 0;
% yline(x,':',LineWidth=2) ;
% fontsize(gca,12,"points") ;
% xlabel('Time (s)','FontSize',14) ;
% ylabel('time derivative (a.u.)','FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% title('Derivative of average arterial pulse wave');
% axis tight

% computation of average arterial pulse wave parameters
T_syst = T(idx_sys);
systole_area = sum(avgArterialPulse(1:idx_sys)) ;
diastole_area = sum(avgArterialPulse(idx_sys:end));
tmp = systole_area;
systole_area = systole_area / (diastole_area + systole_area);
diastole_area = diastole_area / (diastole_area + tmp);
nb_of_detected_systoles = size(sys_index_list,2) ;
[max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
[min_diff_pulse, idx_T_diff_min] = min(diff_pulse_dia);%faire ca mieux
if (~exist("nb_of_averaged_cycles"))
    nb_of_averaged_cycles = 0;
end

% txt file output with measured pulse wave parameters
fileID = fopen(fullfile(one_cycle_dir,strcat(filename,'_pulseWaveParameters.txt')),'w') ;
fprintf(fileID,[...
    'Value of pulse derivative at the maximum systolic increase :\n%d\n' ...
    'Maximal acceleration (m/s^(-2)) :\n%d\n' ...
    'Time of maximum systolic increase (s) :\n%d\n' ...
    'Time of systolic peak (s) :\n%d\n' ...
    'Time of the intersection between the diastolic descent and the threshold 1/e (s) :\n%d\n' ...
    'Number of detected systoles :\n%d\nNumber of averaged cycles :\n%d\n' ...
    'Area under the systolic rise curve :\n%d\n' ...
    'Area under the diastolic descent curve  :\n%d\n' ...
    'Average arterial resistivity index :\n%d\n'], ...
    max_diff_pulse, ...
    acc, ...
    T(idx_T_diff_max+1), ...
    T_syst, ...
    T(idx_sys+idx_list_threshold_dia(1)), ...
    nb_of_detected_systoles, ...
    nb_of_averaged_cycles, ...
    systole_area, ...
    diastole_area, ...
    ARI);
fclose(fileID) ;



% png
% print('-f2','-dpng',fullfile(one_cycle_dir,strcat(filename,'_pulseVsBackground.png'))) ;
% print('-f8','-dpng',fullfile(one_cycle_dir,strcat(filename,'_filteredPulse.png'))) ;
print('-f9','-dpng',fullfile(one_cycle_dir,strcat(filename,'_filteredPulseVsResidual.png'))) ;
print('-f22','-dpng',fullfile(one_cycle_dir,strcat(filename,'_regularizedPulse.png'))) ;
print('-f6','-dpng',fullfile(one_cycle_dir,strcat(filename,'_derivative.png'))) ;
print('-f111','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWave.png'))) ;
print('-f80','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveLabeled.png'))) ;
% print('-f90','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveDerivative.png'))) ;
% print('-f15','-dpng',fullfile(one_cycle_dir,strcat(filename,'_resistivityMap.png'))) ;
% print('-f45','-dpng',fullfile(one_cycle_dir,strcat(filename,'_diastoleHeatMap.png'))) ;
% print('-f46','-dpng',fullfile(one_cycle_dir,strcat(filename,'_systoleHeatMap.png'))) ;
% print('-f23','-dpng',fullfile(one_cycle_dir,strcat(filename,'_rawDopplerHeatMap.png'))) ;
% print('-f24','-dpng',fullfile(one_cycle_dir,strcat(filename,'_flattenedDopplerHeatMap.png'))) ;
% % print('-f77','-dpng',fullfile(one_cycle_dir,strcat(filename,'_zeroLagXcorr.png'))) ;
% % print('-f99','-dpng',fullfile(one_cycle_dir,strcat(filename,'_timeLags.png'))) ;

% eps
% print('-f2','-depsc',fullfile(one_cycle_dir,strcat(filename,'_pulseVsBackground.eps'))) ;
% print('-f8','-depsc',fullfile(one_cycle_dir,strcat(filename,'_filteredPulse.eps'))) ;
print('-f9','-depsc',fullfile(one_cycle_dir,strcat(filename,'_filteredPulseVsResidual.eps'))) ;
print('-f22','-depsc',fullfile(one_cycle_dir,strcat(filename,'_regularizedPulse.eps'))) ;
print('-f6','-depsc',fullfile(one_cycle_dir,strcat(filename,'_derivative.eps'))) ;
print('-f111','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWave.eps'))) ;
print('-f80','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveLabeled.eps'))) ;
% print('-f90','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveDerivative.eps'))) ;
% print('-f15','-depsc',fullfile(one_cycle_dir,strcat(filename,'_resistivityMap.eps'))) ;
% print('-f45','-depsc',fullfile(one_cycle_dir,strcat(filename,'_diastoleHeatMap.eps'))) ;
% print('-f46','-depsc',fullfile(one_cycle_dir,strcat(filename,'_systoleHeatMap.eps'))) ;
% print('-f23','-depsc',fullfile(one_cycle_dir,strcat(filename,'_rawDopplerHeatMap.eps'))) ;
% print('-f24','-depsc',fullfile(one_cycle_dir,strcat(filename,'_flattenedDopplerHeatMap.eps'))) ;
% % print('-f77','-depsc',fullfile(one_cycle_dir,strcat(filename,'_zeroLagXcorr.eps'))) ;
% % print('-f99','-depsc',fullfile(one_cycle_dir,strcat(filename,'_timeLags.eps'))) ;

% masks
imwrite(mat2gray(single(maskArtery)),fullfile(one_cycle_dir,strcat(filename,'_maskArtery.png')),'png') ;
imwrite(mat2gray(single(maskVessel)),fullfile(one_cycle_dir,strcat(filename,'_maskVessel.png')),'png') ;
imwrite(mat2gray(single(maskBackground)),fullfile(one_cycle_dir,strcat(filename,'_maskBackground.png')),'png') ;


displaySuccessMsg();

%% AVG analysis (CRA & CRV pulse)
% figure(123)
% imagesc(squeeze(mean(dataCubeM1M0,3))) ;
% colormap gray
% title('raw AVG frequency map');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image
% c = colorbar('southoutside');
% c.Label.String = 'AVG Doppler frequency (kHz)';
% c.Label.FontSize = 12;
% 
% % dMap = squeeze(mean(fullVideo,3));
% % dMap = flat_field_correction(dMap, ceil(0.07*size(dMap,1)), .33)
% % 
% % figure(124)
% % imagesc(dMap) ;
% % colormap gray
% % title('flattened RMS frequency map');
% % fontsize(gca,12,"points") ;
% % set(gca, 'LineWidth', 2);
% % axis off
% % axis image
% % c = colorbar('southoutside');
% % c.Label.String = 'RMS Doppler frequency (kHz)';
% % c.Label.FontSize = 12;
% % range0(1:2) = clim;
% % % range0 = [1.2*range0(1),0.8*range0(2)];
% % figure(23)
% % clim(range0);
% 
% [maskCRA,maskCRV,maskBackgroundM1M0] = createCentralRetinalArteryVeinMask(dataCubeM1M0);
% 
% % FIXME maskArtery : best computed from reference video
% % maskArtery = createArteryMask(fullVideo);
% figure(105)
% imagesc(maskCRA);
% title('segmented CRA');
% axis off
% axis equal
% colormap gray
% 
% %
% figure(104)
% imagesc(maskCRV);
% title('segmented CRV');
% axis off
% axis equal
% colormap gray
% 
% figure(107)
% imagesc(maskBackgroundM1M0);
% title('segmented BackgroundAVG');
% axis off
% axis equal
% colormap gray

% %fig 102 108
% 
% fullCRAPulse = dataCubeM1M0 .* maskCRA;
% fullCRAPulse = squeeze(sum(fullCRAPulse, [1 2]))/nnz(maskCRA);
% 
% fullCRVPulse = dataCubeM1M0 .* maskCRV;
% fullCRVPulse = squeeze(sum(fullCRVPulse, [1 2]))/nnz(maskCRV);
% 
% %
% fullBackgroundM1M0Signal = dataCubeM1M0 .* maskBackgroundM1M0;
% fullBackgroundM1M0Signal = squeeze(sum(fullBackgroundM1M0Signal, [1 2]))/nnz(maskBackgroundM1M0);
% 
% 
% figure(102)
% plot(fullTime,fullCRAPulse,'-k', fullTime,fullBackgroundM1M0Signal,':k','LineWidth',2) ;
% title('CRA pulse waveform and backgroundAVG signal'); % averaged outside of segmented vessels
% legend('CRA pulse','backgroundAVG') ;
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
% 
% figure(103)
% plot(fullTime,fullCRVPulse,'-k', fullTime,fullBackgroundM1M0Signal,':k','LineWidth',2) ;
% title('CRV pulse waveform and backgroundAVG signal'); % averaged outside of segmented vessels
% legend('CRV pulse','backgroundAVG') ;
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
% 
% % Renormalization to avoid bias
% fullCRAPulseMinusBackground = fullCRAPulse - fullBackgroundM1M0Signal;
% fullCRAPulse = fullCRAPulseMinusBackground;
% fullCRVPulseMinusBackground = fullCRVPulse - fullBackgroundM1M0Signal;
% fullCRVPulse = fullCRVPulseMinusBackground;
% 
% % le plus simple fonctionne bien : soustraire le bkg.
% A = ones(size(dataCubeM1M0));
% % B = ones(size(dataCube));
% for pp = 1:size(dataCubeM1M0,3)
% %     A(:,:,pp) = A(:,:,pp) * fullBackgroundSignal(pp) * avgFullArterialPulse / avgFullBackgroundSignal;
%       A(:,:,pp) = A(:,:,pp) * fullBackgroundM1M0Signal(pp);
% %     B(:,:,pp) = B(:,:,pp) * fullBackgroundSignal(pp);    
% end
% dataCubeM1M0 = dataCubeM1M0 - A ;
% % dataCube = dataCube - A .* maskVessel;
% % dataCube = dataCube - B .* maskBackground;
% 
% 
% % remove outliers
% % 1st pass
% disp('remove outliers... 1st pass.');
% [idxOutPw,fullCRAPulseRmOut] = discardPulseWaveOutliers(fullCRAPulse,3);
% [idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundM1M0Signal,3);
% dataReliabilityIndex1 = ceil(100*(1-0.5*(length(idxOutPw)/length(fullCRAPulse) + length(idxOutBkg)/length(fullBackgroundM1M0Signal))));
% disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);
% 
% [idxOutPw,fullCRVPulseRmOut] = discardPulseWaveOutliers(fullCRVPulse,3);
% [idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundM1M0Signal,3);
% dataReliabilityIndex1 = ceil(100*(1-0.5*(length(idxOutPw)/length(fullCRVPulse) + length(idxOutBkg)/length(fullBackgroundM1M0Signal))));
% disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);
% 
% % smooth trendline data by iterative local linear regression.
% % fullArterialPulseClean = smoothdata(fullArterialPulse,'rloess');
% fullCRAPulseClean = smoothdata(fullCRAPulseRmOut,'lowess');
% fullCRVPulseClean = smoothdata(fullCRVPulseRmOut,'lowess');
% fullBackgroundSignalClean = smoothdata(fullBackgroundSignalRmOut,'lowess');
% 
% 
% 
% 
% 
% figure(108)
% plot(fullTime,fullCRAPulse,':k', ...
%     fullTime,fullCRAPulseClean,'-k', ...
%     'LineWidth',2) ;
% title('CRA pulse minus backgroundAVG vs. filtered pulse');
% legend('<p(t)> - <b(t)>','local linear regression');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
% 
% figure(110)
% plot(fullTime,fullCRVPulse,':k', ...
%     fullTime,fullCRVPulseClean,'-k', ...
%     'LineWidth',2) ;
% title('CRV pulse minus backgroundAVG vs. filtered pulse');
% legend('<p(t)> - <b(t)>','local linear regression');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;


% filename(end-2:end) = 'AVG'; % Save of AVG files here
% 
% % png
% print('-f102','-dpng',fullfile(one_cycle_dir,strcat(filename,'_CRApulseVsBackground.png'))) ;
% print('-f103','-dpng',fullfile(one_cycle_dir,strcat(filename,'_CRVpulseVsBackground.png'))) ;
% print('-f108','-dpng',fullfile(one_cycle_dir,strcat(filename,'_CRAfilteredPulse.png'))) ;
% print('-f110','-dpng',fullfile(one_cycle_dir,strcat(filename,'_CRVfilteredPulse.png'))) ;
% % print('-f9','-dpng',fullfile(one_cycle_dir,strcat(filename,'_filteredPulseVsResidual.png'))) ;
% % print('-f22','-dpng',fullfile(one_cycle_dir,strcat(filename,'_regularizedPulse.png'))) ;
% % print('-f6','-dpng',fullfile(one_cycle_dir,strcat(filename,'_derivative.png'))) ;
% % print('-f1','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWave.png'))) ;
% % print('-f80','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveLabeled.png'))) ;
% % print('-f90','-dpng',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveDerivative.png'))) ;
% % print('-f15','-dpng',fullfile(one_cycle_dir,strcat(filename,'_resistivityMap.png'))) ;
% % print('-f45','-dpng',fullfile(one_cycle_dir,strcat(filename,'_diastoleHeatMap.png'))) ;
% % print('-f46','-dpng',fullfile(one_cycle_dir,strcat(filename,'_systoleHeatMap.png'))) ;
% print('-f123','-dpng',fullfile(one_cycle_dir,strcat(filename,'_rawDopplerHeatMap.png'))) ;
% % print('-f24','-dpng',fullfile(one_cycle_dir,strcat(filename,'_flattenedDopplerHeatMap.png'))) ;
% % % print('-f77','-dpng',fullfile(one_cycle_dir,strcat(filename,'_zeroLagXcorr.png'))) ;
% % % print('-f99','-dpng',fullfile(one_cycle_dir,strcat(filename,'_timeLags.png'))) ;
% 
% % eps
% print('-f102','-depsc',fullfile(one_cycle_dir,strcat(filename,'_CRApulseVsBackground.eps'))) ;
% print('-f103','-depsc',fullfile(one_cycle_dir,strcat(filename,'_CRVpulseVsBackground.eps'))) ;
% print('-f108','-depsc',fullfile(one_cycle_dir,strcat(filename,'_CRAfilteredPulse.eps'))) ;
% print('-f110','-depsc',fullfile(one_cycle_dir,strcat(filename,'_CRVfilteredPulse.eps'))) ;
% % print('-f9','-depsc',fullfile(one_cycle_dir,strcat(filename,'_filteredPulseVsResidual.eps'))) ;
% % print('-f22','-depsc',fullfile(one_cycle_dir,strcat(filename,'_regularizedPulse.eps'))) ;
% % print('-f6','-depsc',fullfile(one_cycle_dir,strcat(filename,'_derivative.eps'))) ;
% % print('-f1','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWave.eps'))) ;
% % print('-f80','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveLabeled.eps'))) ;
% % print('-f90','-depsc',fullfile(one_cycle_dir,strcat(filename,'_avgPulseWaveDerivative.eps'))) ;
% % print('-f15','-depsc',fullfile(one_cycle_dir,strcat(filename,'_resistivityMap.eps'))) ;
% % print('-f45','-depsc',fullfile(one_cycle_dir,strcat(filename,'_diastoleHeatMap.eps'))) ;
% % print('-f46','-depsc',fullfile(one_cycle_dir,strcat(filename,'_systoleHeatMap.eps'))) ;
% print('-f123','-depsc',fullfile(one_cycle_dir,strcat(filename,'_rawDopplerAvgHeatMap.eps'))) ;
% % print('-f24','-depsc',fullfile(one_cycle_dir,strcat(filename,'_flattenedDopplerHeatMap.eps'))) ;
% % % print('-f77','-depsc',fullfile(one_cycle_dir,strcat(filename,'_zeroLagXcorr.eps'))) ;
% % % print('-f99','-depsc',fullfile(one_cycle_dir,strcat(filename,'_timeLags.eps'))) ;
% 
% % masks
% imwrite(mat2gray(single(maskCRA)),fullfile(one_cycle_dir,strcat(filename,'_maskCRA.png')),'png') ;
% imwrite(mat2gray(single(maskCRV)),fullfile(one_cycle_dir,strcat(filename,'_maskCRV.png')),'png') ;
% imwrite(mat2gray(single(maskBackgroundM1M0)),fullfile(one_cycle_dir,strcat(filename,'_maskBackground.png')),'png') ;
% 
% 
% displaySuccessMsg();
%%
return;

end