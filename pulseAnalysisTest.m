function [v_RMS] = pulseAnalysisTest(Ninterp, fullVideoM2M0, fullVideoM1M0,sys_index_list, maskArtery,maskVessel,maskVein,maskBackground ,ToolBox,path)
PW_params = Parameters(path);
PW_params = Parameters(path);

dataCubeM2M0 = fullVideoM2M0;
dataCubeM1M0 = fullVideoM1M0;


% for robust rendering : 
% 1-flat-field correction, 2-background substraction
for pp = 1:size(dataCubeM2M0,3)
    dataCubeM2M0(:,:,pp) = flat_field_correction(squeeze(dataCubeM2M0(:,:,pp)), PW_params.flatField_gwRatio*size(dataCubeM2M0,1), PW_params.flatField_border);
end


%% calculate raw signals of arteries, background and veins

fullArterialPulse = fullVideoM2M0 .* maskArtery;
fullArterialPulse = squeeze(sum(fullArterialPulse, [1 2]))/nnz(maskArtery);

fullBackgroundSignal = fullVideoM2M0 .* maskBackground;
fullBackgroundSignal = squeeze(sum(fullBackgroundSignal, [1 2]))/nnz(maskBackground);

fullVenousSignal = fullVideoM2M0 .* maskVein;
fullVenousSignal = squeeze(sum(fullVenousSignal, [1 2]))/nnz(maskVein);

fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal;

figure(122)
plot(fullArterialPulseMinusBackground)
title('Pulse-background before flatfield');

ArterialPulse = dataCubeM2M0 .* maskArtery;
ArterialPulse = squeeze(sum(ArterialPulse, [1 2]))/nnz(maskArtery);

BackgroundSignal = dataCubeM2M0 .* maskBackground;
BackgroundSignal = squeeze(sum(BackgroundSignal, [1 2]))/nnz(maskBackground);

VenousSignal = dataCubeM2M0 .* maskVein;
VenousSignal = squeeze(sum(VenousSignal, [1 2]))/nnz(maskVein);

ArterialPulseMinusBackground = ArterialPulse - BackgroundSignal;

figure(123)
plot(ArterialPulseMinusBackground)
title('Pulse-background after flatfield');

mask = imdilate(maskArtery,strel('disk',15));
LocaclBGK = zeros(size(dataCubeM2M0));
for n = 1:size(dataCubeM2M0,3)
    LocaclBGK(:,:,n) = regionfill(fullVideoM2M0(:,:,n),mask);
end
LocaclBGK = LocaclBGK.*maskArtery;
LocaclBGK = squeeze(sum(LocaclBGK, [1 2]))/nnz(maskArtery);

figure(124)
plot(fullArterialPulse-LocaclBGK)
title('Pulse-background local before flatfield');


for n = 1:size(dataCubeM2M0,3)
    LocaclBGK2(:,:,n) = regionfill(dataCubeM2M0(:,:,n),mask);
end
LocaclBGK2 = LocaclBGK2.*maskArtery;
LocaclBGK2 = squeeze(sum(LocaclBGK2, [1 2]))/nnz(maskArtery);

figure(125)
plot(ArterialPulse-LocaclBGK2)
title('Pulse-background local after flatfield');

figure(126)
plot(ArterialPulse-fullBackgroundSignal)
title('OLD');
v_RMS = zeros(size(fullVideoM2M0));
print('-f122','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_pulseMinBackground_BeforeFF.png'))) ;
print('-f123','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_pulseMinBackground_AfterFF.png'))) ;
print('-f124','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_LocalpulseMinBackground_BeforeFF.png'))) ;
print('-f125','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_LocalpulseMinBackground_AfterFF.png'))) ;
print('-f126','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_pulseMinBackgroundOLD.png'))) ;

end