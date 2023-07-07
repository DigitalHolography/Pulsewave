function [] = spectrogram(maskArtery,maskBackground ,SH,ToolBox)

SH_arteries = SH.*maskArtery;
SH_background = SH.*maskBackground;

ArterySpectrum = zeros(1,size(SH_arteries,3));
BgSpectrum = zeros(1,size(SH_arteries,3));
frq_shift = linspace(-ToolBox.fs/2,ToolBox.fs/2,size(SH_arteries,3));

for ii = 1:size(SH_arteries,3)
   ArterySpectrum(ii) = sum(SH_arteries(:,:,ii),'all')/nnz(SH_arteries(:,:,ii));
   BgSpectrum(ii) = sum(SH_background(:,:,ii),'all')/nnz(SH_background(:,:,ii));
end

ArterySpectrum = circshift(ArterySpectrum,round(size(ArterySpectrum)/2));
BgSpectrum = circshift(BgSpectrum,round(size(BgSpectrum)/2));
DeltaSpectrum = ArterySpectrum - BgSpectrum;

mamax = 1.1*max (ArterySpectrum);
mimin = 0.9*min(DeltaSpectrum(1:round(0.25*size(DeltaSpectrum,2))));


% spectrumAVG = squeeze(mean(SH,[1 2]));
% spectrumAVG_calcul_omega = mean(SH,[1 2]);
% momentM2 = moment2(spectrumAVG_calcul_omega, f1, f2, fs, size(SH,3), 0);
% momentM0 = moment0(spectrumAVG_calcul_omega, f1, f2, fs, size(SH,3), 0);
% omegaRMS=squeeze(sqrt(momentM2./momentM0));
% omegaRMS_index=omegaRMS*size(SH,3)/fs;
% I_omega=log(spectrumAVG(round(omegaRMS_index)));


figure(19)
p1 = semilogy(frq_shift ,ArterySpectrum,'DisplayName','ArterySpectrum') ;
hold on 
p2 = semilogy(frq_shift ,BgSpectrum,'DisplayName','BgSpectrum') ;
hold on
p3 = semilogy(frq_shift ,DeltaSpectrum,'DisplayName','DeltaSpectrum') ;


ylim([mimin mamax])
xlim([-ToolBox.fs/2 ToolBox.fs/2])
rectangle('Position', [-ToolBox.f1 mimin  2*ToolBox.f1 (mamax-mimin)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [-ToolBox.fs/2 mimin (ToolBox.fs/2-ToolBox.f2) (mamax-mimin)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [ToolBox.f2 mimin (ToolBox.fs/2-ToolBox.f2) (mamax-mimin)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
% om_RMS_line = line([-omegaRMS omegaRMS],[I_omega I_omega]);
% om_RMS_line.Color = 'black';
% om_RMS_line.LineStyle = '--';
% om_RMS_line.LineWidth = 1;
% om_RMS_line.Tag = 'omega RMS';
% text(0, I_omega, '\bf{\omega_{RMS}}', 'HorizontalAlignment','center', 'VerticalAlignment','top')


xline(ToolBox.f1,'--','HandleVisibility','off')
xline(ToolBox.f2,'--','HandleVisibility','off')
xline(-ToolBox.f1,'--','HandleVisibility','off')
xline(-ToolBox.f2,'--','HandleVisibility','off')
xticks([-ToolBox.fs/2 -ToolBox.f2 -ToolBox.f1 ToolBox.f1 ToolBox.f2 ToolBox.fs/2])
xticklabels({'-f_s/2','-f_2','-f_1','f_1','f_2','f_s/2'})
title('Average spectrum')

fontsize(gca,12,"points") ;
xlabel('frequency (kHz)','FontSize',14) ;
ylabel('A.U.','FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 1);
uistack(p1,'top');
uistack(p2,'top');
uistack(p3,'top');

hold off
legend

print('-f19','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Average spectrogramms.png'))) ;
close all


