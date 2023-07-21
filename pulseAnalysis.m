function [v_RMS] = pulseAnalysis(Ninterp, fullVideoM2M0, fullVideoM1M0,sys_index_list, maskArtery,maskVein,maskBackground ,ToolBox,path)
PW_params = Parameters(path);

dataCubeM2M0 = fullVideoM2M0;
dataCubeM1M0 = fullVideoM1M0;


% for robust rendering : 
% 1-flat-field correction, 2-background substraction
for pp = 1:size(dataCubeM2M0,3)
    dataCubeM2M0(:,:,pp) = flat_field_correction(squeeze(dataCubeM2M0(:,:,pp)), PW_params.flatField_gwRatio*size(dataCubeM2M0,1), PW_params.flatField_border);
end


%% calculate raw signals of arteries, background and veins

% fullArterialPulse = fullVideoM2M0 .* maskArtery;
% fullArterialPulse = squeeze(sum(fullArterialPulse, [1 2]))/nnz(maskArtery);
% 
% fullBackgroundSignal = fullVideoM2M0 .* maskBackground;
% fullBackgroundSignal = squeeze(sum(fullBackgroundSignal, [1 2]))/nnz(maskBackground);
% 
% fullVenousSignal = fullVideoM2M0 .* maskVein;
% fullVenousSignal = squeeze(sum(fullVenousSignal, [1 2]))/nnz(maskVein);

fullArterialPulse = dataCubeM2M0 .* maskArtery;
fullArterialPulse = squeeze(sum(fullArterialPulse, [1 2]))/nnz(maskArtery);

fullBackgroundSignal = dataCubeM2M0 .* maskBackground;
fullBackgroundSignal = squeeze(sum(fullBackgroundSignal, [1 2]))/nnz(maskBackground);

fullVenousSignal = dataCubeM2M0 .* maskVein;
fullVenousSignal = squeeze(sum(fullVenousSignal, [1 2]))/nnz(maskVein);



ArterialPulse = dataCubeM2M0 .* maskArtery;
ArterialPulse = squeeze(sum(ArterialPulse, [1 2]))/nnz(maskArtery);

BackgroundSignal = dataCubeM2M0 .* maskBackground;
BackgroundSignal = squeeze(sum(BackgroundSignal, [1 2]))/nnz(maskBackground);

VenousSignal = dataCubeM2M0 .* maskVein;
VenousSignal = squeeze(sum(VenousSignal, [1 2]))/nnz(maskVein);


%% cleaning signals

% Renormalization to avoid bias
% renormalize avg. background to avg. as local pulse
% avgFullBackgroundSignal = mean(fullBackgroundSignal(:));
% avgFullArterialPulse = mean(fullArterialPulse(:));
% fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal * avgFullArterialPulse / avgFullBackgroundSignal;
fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal;
fullVenousSignalMinusBackground = fullVenousSignal - fullBackgroundSignal;
%fullArterialPulse = fullArterialPulseMinusBackground;

% le plus simple fonctionne bien : soustraire le bkg.
A = ones(size(dataCubeM2M0));

for pp = 1:size(dataCubeM2M0,3)
      A(:,:,pp) = A(:,:,pp) * BackgroundSignal(pp);
end
dataCubeM2M0 = dataCubeM2M0 - A ; 

% remove outliers
% 1st pass
disp('remove outliers... 1st pass.');
[idxOutPw,fullArterialPulseRmOut] = discardPulseWaveOutliers(fullArterialPulseMinusBackground,3);
[idxOutVn,fullVenousSignalRmOut] = discardPulseWaveOutliers(fullVenousSignalMinusBackground,3);
[idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundSignal,3);
% FIXME

dataReliabilityIndex1 = ceil(100*(1-PW_params.pulseAnal_dataReliabilityFactor*(length(idxOutPw)/length(fullArterialPulse) + length(idxOutBkg)/length(fullBackgroundSignal))));
disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);

% smooth trendline data by iterative local linear regression.

fullArterialPulseClean = smoothdata(fullArterialPulseRmOut,'lowess');
fullVenousSignalClean = smoothdata(fullVenousSignalRmOut, 'lowess');
% fullBackgroundSignalClean = smoothdata(fullBackgroundSignalRmOut,'lowess');


%% calculate the pulse derivative and finding/cleaning pulses

fullArterialPulseDerivative = diff(fullArterialPulse);
fullArterialPulseCleanDerivative = diff(fullArterialPulseClean);

% now cleanup dataCube to create_one_cycle()
% strategy : use average pulse profiles to detect and
% find  noisy frames @ >3 std from zero-mean
% replace them with cleaned data
noise = sqrt(abs(abs(fullArterialPulseMinusBackground).^2 - abs(fullArterialPulseClean).^2));
idxOutNoise = find(noise>PW_params.pulseAnal_outNoiseThreshold*std(noise));

dataReliabilityIndex2 = ceil(100*(1-(length(idxOutNoise)/length(fullArterialPulseClean) )));
disp(['data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);


% FIXME : compute true regularied cube by replacing bad frames
fullArterialPulseRegularized = squeeze(sum(dataCubeM2M0 .* maskArtery, [1 2])) / nnz(maskArtery);


%% Creation of the avg pluse for In-plane arteries

[onePulseVideo, selectedPulseIdx, cycles_signal] = create_one_cycle(dataCubeM2M0, maskArtery, sys_index_list, Ninterp,path);

avgArterialPulseHz = squeeze(sum(onePulseVideo .* maskArtery, [1 2]))/nnz(maskArtery);
avgArterialPulseVelocityInPlane = avgArterialPulseHz * ToolBox.ScalingFactorVelocityInPlane;


v_RMS = onePulseVideo * ToolBox.ScalingFactorVelocityInPlane;

% avi
w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_one_cycle.avi')));
tmp = mat2gray(onePulseVideo);
open(w)
for j = 1:size(onePulseVideo,3)
    writeVideo(w,tmp(:,:,j)) ;  
end
close(w);
% mp4
w = VideoWriter(fullfile(ToolBox.PW_path_mp4,strcat(ToolBox.main_foldername,'_one_cycle.mp4')),'MPEG-4');
tmp = mat2gray(onePulseVideo);
open(w)
for j = 1:size(onePulseVideo,3)
    writeVideo(w,tmp(:,:,j)) ;  
end
close(w);

%FIXME: M1/M0 and M2/M0 are subject to aliases at 67 kHz

blur_time_sys = ceil(Ninterp/PW_params.pulseAnal_blurScaleFactor);
blur_time_dia = ceil(Ninterp/PW_params.pulseAnal_blurScaleFactor);


average_cycle_length = 0;
nb_of_averaged_cycles = 0;
if size(sys_index_list,2) == 1
    average_cycle_length = Ninterp ;
    nb_of_averaged_cycles = 1;
else
    for ii = 2:size(sys_index_list,2)
        average_cycle_length = average_cycle_length + (sys_index_list(ii)-sys_index_list(ii-1)) ;
        nb_of_averaged_cycles = nb_of_averaged_cycles + 1;
    end
    average_cycle_length = average_cycle_length / (length(sys_index_list)-1);
end
T = linspace(0,ToolBox.stride/ToolBox.fs*average_cycle_length,Ninterp);



% save average arterial pulse wave velocity to txt file
tmp = [T(1:length(avgArterialPulseHz))',avgArterialPulseVelocityInPlane];
%size(tmp)
fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername,'_avgPulse.txt')),'w') ;
fprintf(fileID,'%f %f \r\n',tmp');
fclose(fileID);


<<<<<<< HEAD
meanIm = mat2gray(squeeze(mean(dataCubeM2M0,3)));
tolVal = [0.02, 0.98]; 
meanIm = mat2gray(imadjust(meanIm, stretchlim(meanIm, tolVal)));
=======
%% Arterial resistivity calculation
disp('arterial resistivity...');
[ARImap, ARI, ARImapRGB, ARIvideoRGB, gamma, img_avg] = construct_resistivity_index(onePulseVideo, maskArtery,path);
ARImap = ARImap.*maskArtery;

% avi
w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ARIvideoRGB.avi')));
open(w)
ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
end
close(w);

% mp4
w = VideoWriter(fullfile(ToolBox.PW_path_mp4,strcat(ToolBox.main_foldername,'_ARIvideoRGB.mp4')),'MPEG-4');
open(w)
ARIvideoRGB = im2uint8(mat2gray(ARIvideoRGB));
for jj = 1:size(ARIvideoRGB,4) % ARIvideoRGB is four dimensional: height-by-width-by-3-by-frames
    writeVideo(w,squeeze(ARIvideoRGB(:,:,:,jj))) ;
end
close(w);

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


%% artery mask by lag-correlation with average pulse
% arteries = zeros(size(id_max));
% % arteries (or(id_max < 0.1*Ninterp, id_max > 0.9 * Ninterp)) = 1 ;
% arteries (or(id_max <= ceil(0.01*Ninterp), id_max >= floor(0.95 * Ninterp))) = 1 ;
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


%% vein mask by lag-correlation with average pulse
% veins = zeros(size(id_max)) ;
% veins(and(id_max>=ceil(0.1*Ninterp),id_max<=floor(0.3*Ninterp))) = 1 ;
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
% I_arteries = one_pulse_video .* maskArtery;
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
% sys_index_list_one_cycle = find_systole_index(onePulseVideo2);
>>>>>>> 09984893e24d0e40f6cd744529583ac45f31b4d6

[~,idx_sys] = max(avgArterialPulseHz) ;

%% diastolic Doppler frequency heatmap : 10% of frames before minimum of diastole
heatmap_dia = squeeze(mean(onePulseVideo(:,:,floor(0.9*Ninterp):Ninterp),3));
% onePulseVideo2 : no background correction 
% heatmap_dia = squeeze(mean(onePulseVideo2(:,:,floor(0.9*Ninterp):Ninterp),3));
heatmap_dia = flat_field_correction(heatmap_dia, ceil(PW_params.flatField_gwRatio*size(heatmap_dia,1)), PW_params.flatField_borderDMap);


%% systolic Doppler frequency heatmap : 10% of frames around peak systole
a = max(ceil(idx_sys-0.05*Ninterp),1);
b = min(ceil(idx_sys+0.05*Ninterp),Ninterp);
heatmap_sys = squeeze(mean(onePulseVideo(:,:,a:b),3));
% onePulseVideo2 : no background correction 
% heatmap_sys = squeeze(mean(onePulseVideo2(:,:,a:b),3));
heatmap_sys = flat_field_correction(heatmap_sys, ceil(PW_params.flatField_gwRatio*size(heatmap_sys,1)), PW_params.flatField_borderDMap);

%%
% FIXME : replace sys + dia blur by homogenous blur ? 
pulse_arteries_blurred_sys = movavgvar(avgArterialPulseHz(1:idx_sys), blur_time_sys);
diff_pulse_sys = diff(pulse_arteries_blurred_sys) ;
pulse_arteries_blurred_dia = movavgvar(avgArterialPulseHz(idx_sys:end), blur_time_dia);
diff_pulse_dia = diff(pulse_arteries_blurred_dia) ;
diff_avgPulse = diff(movavgvar(avgArterialPulseHz, blur_time_sys)) * ToolBox.ScalingFactorVelocityInPlane*Ninterp/T(end);
delta = max(pulse_arteries_blurred_dia(:))-min(pulse_arteries_blurred_dia(:));
thr = max(pulse_arteries_blurred_dia(:))-delta./exp(1);
idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < thr);
[max_diff_pulse,idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
if idx_T_diff_max == 1 
    acc = abs(max_diff_pulse/(T(idx_T_diff_max)-T(idx_T_diff_max+1)));
else
    acc = max_diff_pulse/(T(idx_T_diff_max)-T(idx_T_diff_max-1));
end


% computation of average arterial pulse wave parameters
T_syst = T(idx_sys);
systole_area = sum(avgArterialPulseHz(1:idx_sys)) ;
diastole_area = sum(avgArterialPulseHz(idx_sys:end));
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
fileID = fopen(fullfile(ToolBox.PW_path_txt,strcat(ToolBox.main_foldername,'_pulseWaveParameters.txt')),'w') ;
fprintf(fileID,[...
    'Value of pulse derivative at the maximum systolic increase :\n%d\n' ...
    'Maximal acceleration (m/s^(-2)) :\n%d\n' ...
    'Time of maximum systolic increase (s) :\n%d\n' ...
    'Time of systolic peak (s) :\n%d\n' ...
    'Time of the intersection between the diastolic descent and the threshold 1/e (s) :\n%d\n' ...
    'Number of detected systoles :\n%d\nNumber of averaged cycles :\n%d\n' ...
    'Area under the systolic rise curve :\n%d\n' ...
    'Area under the diastolic descent curve  :\n%d\n'], ...
    max_diff_pulse, ...
    acc, ...
    T(idx_T_diff_max+1), ...
    T_syst, ...
    T(idx_sys+idx_list_threshold_dia(1)), ...
    nb_of_detected_systoles, ...
    nb_of_averaged_cycles, ...
    systole_area, ...
    diastole_area);
fclose(fileID) ;

segmentation_map = zeros(size(meanIm,1), size(meanIm,2), 3);
segmentation_map(:,:, 1) = meanIm - (maskArtery+maskVein).*meanIm + maskArtery;
segmentation_map(:,:, 2) = meanIm - (maskArtery+maskVein).*meanIm;
segmentation_map(:,:, 3) = meanIm - (maskArtery+maskVein).*meanIm + maskVein;




%% Display images and figures
strXlabel = 'Time(ms)' ;%createXlabelTime(1);

strYlabel = 'frequency (kHz)';
range0(1:2) = clim;
fullTime = linspace(0,size(fullVideoM2M0,3)*ToolBox.stride/ToolBox.fs,size(fullVideoM2M0,3));

% calculate raw signals of arteries, background and veins

figure(20)
plot(fullTime,fullArterialPulse,'-k', fullTime,fullBackgroundSignal,':k', fullTime, fullVenousSignal, '-.k', 'LineWidth',2) ;
title('arterial pulse waveform and background signal'); % averaged outside of segmented vessels
legend('arterial pulse','background', 'venous signal') ;
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel(strYlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

% cleaning signals
figure(30)
plot(fullTime,fullArterialPulseMinusBackground,':k', ...
    fullTime,fullArterialPulseClean,'-k', ...
    'LineWidth',2) ;
title('arterial pulse minus background vs. filtered pulse');
legend('<p(t)> - <b(t)>','local linear regression');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel(strYlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

figure(31)
plot(fullTime,fullVenousSignalMinusBackground,':k', ...
    fullTime,fullVenousSignalClean,'-k', ...
    'LineWidth',2) ;
title('venous signal minus background vs. filtered pulse');
legend('<p(t)> - <b(t)>','local linear regression');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel(strYlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

% calculate the pulse derivative and finding/cleaning pulses
figure(40)
plot( ...
    fullTime(1:length(fullArterialPulseDerivative)),fullArterialPulseDerivative,':k', ...
    fullTime(1:length(fullArterialPulseCleanDerivative)),fullArterialPulseCleanDerivative,'-k', ...
    'LineWidth',2) ;
text(fullTime(sys_index_list)-0.3,fullArterialPulseCleanDerivative(sys_index_list)+0.03,num2str((1:numel(sys_index_list))'))
title('derivative of the arterial pulse waveform');
legend('\delta <p(t)> - <b(t)>','from smoothed data');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel('A.U.','FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;


c = colorbar('southoutside');

    % Colorbar for raw/flattened image
colorfig = figure(2410);
colorfig.Units = 'normalized';
colormap(c);
colormap gray
clim(range0);
hCB = colorbar('north');
set(gca,'Visible',false)
set(gca,'LineWidth', 3);
hCB.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca,15,"points");
colorTitleHandle = get(hCB,'Title');
titleString = 'RMS Doppler frequency (kHz)';
set(colorTitleHandle ,'String',titleString);

figure(43)   
plot( ...
    fullTime(1:length(fullArterialPulseRegularized)),fullArterialPulseRegularized,'-k', ...
    fullTime(1:length(fullArterialPulse)),fullArterialPulseMinusBackground,':k', ...
    'LineWidth',2) ;
yline(0,':',{''},LineWidth=2) ;
title('Regularized pulse wave signal averaged in retinal arteries');
legend('Regularized arterial pulse','<p(t)> - <b(t)>');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

figure(44)
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


% Colorbar for AVG image
colorfig = figure(1230);
clim(range0);
colorfig.Units = 'normalized';
colormap(c);
colormap gray
hCB = colorbar('north');
set(gca,'Visible',false)
set(gca,'LineWidth', 3);
hCB.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca,15,"points");
colorTitleHandle = get(hCB,'Title');
titleString = 'AVG Doppler frequency (kHz)';
set(colorTitleHandle ,'String',titleString);


% diastolic Doppler frequency heatmap
figure(80)
imagesc(heatmap_dia) ;
colormap gray
title('bottom diastole RMS frequency map');
fontsize(gca,12,"points") ;
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;

% systolic Doppler frequency heatmap
figure(90)
imagesc(heatmap_sys) ;
colormap gray
title('peak systole RMS frequency map');
fontsize(gca,12,"points") ;
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(3:4) = clim;
    % same color axis for systolic and diastolic Doppler heatmaps
clim([min(range),max(range)]);
figure(85)
clim([min(range),max(range)]);

%
figure(100)
plot( ...
    T(1:length(avgArterialPulseHz)),avgArterialPulseVelocityInPlane,'k-', ...
    LineWidth=2) ;
xline(T(idx_T_diff_max + 1),':',{},LineWidth=2)
text(T(idx_T_diff_max + 1),min(avgArterialPulseVelocityInPlane(:))+0.1*(max(avgArterialPulseVelocityInPlane(:))-min(avgArterialPulseVelocityInPlane(:))),' (1)','FontSize',14);%Display at minimum+x%
xline(T_syst,':',{},LineWidth=2) ;
text(T_syst,min(avgArterialPulseVelocityInPlane(:))+0.1*(max(avgArterialPulseVelocityInPlane(:))-min(avgArterialPulseVelocityInPlane(:))),'  (2)','FontSize',14);
xline(T(idx_sys + idx_list_threshold_dia(1)),':',{},LineWidth=2);
text(T(idx_sys + idx_list_threshold_dia(1)),min(avgArterialPulseVelocityInPlane(:))+0.1*(max(avgArterialPulseVelocityInPlane(:))-min(avgArterialPulseVelocityInPlane(:))),'  (3)','FontSize',14);
legend(' arterial pulse');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel('blood flow velocity (mm/s)','FontSize',14);
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;
title('average blood flow velocity estimate in in-plane retinal arteries');

figure(101)
for ii = 1 : size(cycles_signal, 1)
    if ismember(ii, selectedPulseIdx)
        plot( ...
             T,movavgvar(cycles_signal(ii, :),5),'k-', ...
            'LineWidth',1) ;
        hold on
%     else
%         plot( ...
%             T,movavgvar(cycles_signal(ii, :),5),'k--', ...
%             'LineWidth',1) ;
%         hold on
    end

end
title('arterial Doppler signal ');
legend('arterial signal ');
fontsize(gca,12,"points") ;
xlabel(strXlabel,'FontSize',14) ;
ylabel('Doppler signal (kHz)');
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;

figure(102)
plot(T,avgArterialPulseHz,'k.', ...
    T(1:idx_sys),pulse_arteries_blurred_sys(1:idx_sys),'k-', ...
    T(idx_sys:Ninterp),pulse_arteries_blurred_dia(1:(Ninterp-idx_sys+1)),'k-', ...
    LineWidth=2) ;
xline(T(idx_T_diff_max + 1),':',{},LineWidth=2)
text(T(idx_T_diff_max + 1),min(pulse_arteries_blurred_dia(:))+0.1*(max(pulse_arteries_blurred_dia(:))-min(pulse_arteries_blurred_dia(:))),' (1)','FontSize',14);%Display at minimum+x%
xline(T_syst,':',{},LineWidth=2) ;
text(T_syst,min(pulse_arteries_blurred_dia(:))+0.1*(max(pulse_arteries_blurred_dia(:))-min(pulse_arteries_blurred_dia(:))),'  (2)','FontSize',14);
xline(T(idx_sys + idx_list_threshold_dia(1)),':',{},LineWidth=2);
text(T(idx_sys + idx_list_threshold_dia(1)),min(pulse_arteries_blurred_dia(:))+0.1*(max(pulse_arteries_blurred_dia(:))-min(pulse_arteries_blurred_dia(:))),'  (3)','FontSize',14);
%                 yline(1/exp(1),':',LineWidth=2);
legend('arterial pulse','smoothed line') ;
fontsize(gca,12,"points") ;
xlabel('Time (s)','FontSize',14) ;
ylabel('frequency (kHz)','FontSize',14) ;
axis tight;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
title('average background-corrected RMS frequency in retinal arteries');

figure(103)
plot(T(1:end-1), diff_avgPulse,'k-', LineWidth=2);
x = 0;
yline(x,':',LineWidth=2) ;
xline(T(idx_T_diff_max + 1),':',{},LineWidth=2)
text(T(idx_T_diff_max + 1),min(diff_avgPulse(:))+0.1*(max(diff_avgPulse(:))-min(diff_avgPulse(:))),' (1)','FontSize',14);%Display at minimum+x%
xline(T_syst,':',{},LineWidth=2) ;
text(T_syst,min(diff_avgPulse(:))+0.1*(max(diff_avgPulse(:))-min(diff_avgPulse(:))),'  (2)','FontSize',14);
xline(T(idx_sys + idx_list_threshold_dia(1)),':',{},LineWidth=2);
text(T(idx_sys + idx_list_threshold_dia(1)),min(diff_avgPulse(:))+0.1*(max(diff_avgPulse(:))-min(diff_avgPulse(:))),'  (3)','FontSize',14);
legend(' arterial pulse');
fontsize(gca,12,"points") ;
xlabel('Time (s)','FontSize',14) ;
ylabel('time derivative (mm.s^{-2})','FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
title('Derivative of average arterial pulse wave');
axis tight

%% Saving images and figures
% ToolBox.main_foldername(end-2:end) = 'AVG'; % Save of AVG files here
% 
% % png
% print('-f102','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRApulseVsBackground.png'))) ;
% print('-f103','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRVpulseVsBackground.png'))) ;
% print('-f108','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRAfilteredPulse.png'))) ;
% print('-f110','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRVfilteredPulse.png'))) ;
% print('-f24','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_flattenedDopplerHeatMap.png'))) ;
% print('-f77','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_zeroLagXcorr.png'))) ;
% print('-f99','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_timeLags.png'))) ;
% 
% % eps
% print('-f102','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRApulseVsBackground.eps'))) ;
% print('-f103','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRVpulseVsBackground.eps'))) ;
% print('-f108','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRAfilteredPulse.eps'))) ;
% print('-f110','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_CRVfilteredPulse.eps'))) ;
% print('-f77','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_zeroLagXcorr.eps'))) ;
% print('-f99','-depsc',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_timeLags.eps'))) ;

% png
print('-f20','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_pulseVsBackground.png'))) ;
print('-f30','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_filteredPulse.png'))) ;
print('-f44','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_filteredPulseVsResidual.png'))) ;
print('-f43','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_regularizedPulse.png'))) ;
print('-f40','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_derivative.png'))) ;
print('-f100','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_avgPulseWave.png'))) ;
print('-f102','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_avgPulseWaveLabeled.png'))) ;
print('-f103','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_avgPulseWaveDerivative.png'))) ;
print('-f80','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_diastoleHeatMap.png'))) ;
print('-f90','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_systoleHeatMap.png'))) ;
print('-f101','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_all_cycles.png'))) ;
print('-f2410','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_RMS_frequency_colorbar.png')));
print('-f1230','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_frequency_colorbar.png')));



% print('-f77','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_zeroLagXcorr.png'))) ;
% print('-f99','-dpng',fullfile(one_cycle_dir,strcat(ToolBox.main_foldername,'_timeLags.png'))) ;

% eps
print('-f20','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_pulseVsBackground.eps'))) ;
print('-f30','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_filteredPulse.eps'))) ;
print('-f44','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_filteredPulseVsResidual.eps'))) ;
print('-f43','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_regularizedPulse.eps'))) ;
print('-f40','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_derivative.eps'))) ;
print('-f100','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_avgPulseWave.eps'))) ;
print('-f102','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_avgPulseWaveLabeled.eps'))) ;
print('-f103','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_avgPulseWaveDerivative.eps'))) ;
print('-f80','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_diastoleHeatMap.eps'))) ;
print('-f90','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_systoleHeatMap.eps'))) ;
print('-f101','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_all_cycles.eps'))) ;
print('-f2410','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_RMS_frequency_colorbar.eps')));
print('-f1230','-depsc',fullfile(ToolBox.PW_path_eps,strcat(ToolBox.main_foldername,'_AVG_frequency_colorbar.eps')));


imwrite(segmentation_map,fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_artery_vein_segmentation.png')),'png') ;

close all

return;

%% Analysis CRA & CRV Kept in case we want to study CRA/CRV again 
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
% 

end