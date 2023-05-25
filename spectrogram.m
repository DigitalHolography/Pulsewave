function [SpectrogramVideo] = spectrogram(maskArtery, maskVein, maskCRA, SH, one_cycle_dir, filename, k)
% threshhold = round(0.18*size(SH,3));
% SH = SH(:,:,threshhold:(end-threshhold));
one_cycle_dir_png = fullfile(one_cycle_dir, 'png');
SH_arteries = SH.*maskArtery;

maskBG = ((maskVein +maskArtery) == 0);
SH_background = SH.*maskBG;

ArterySpectrum = zeros(1,size(SH_arteries,3));
BgSpectrum = zeros(1,size(SH_arteries,3));
frq_shift = linspace(-33,33,size(SH_arteries,3));

for ii = 1:size(SH_arteries,3)
   ArterySpectrum(ii) = sum(SH_arteries(:,:,ii),'all')/nnz(SH_arteries(:,:,ii));
   BgSpectrum(ii) = sum(SH_background(:,:,ii),'all')/nnz(SH_background(:,:,ii));
end

ArterySpectrum = circshift(ArterySpectrum,round(size(ArterySpectrum)/2));
BgSpectrum = circshift(BgSpectrum,round(size(BgSpectrum)/2));

DeltaSpectrum = ArterySpectrum - BgSpectrum;


mamax = 1.1*max (ArterySpectrum);
mimin = 0.9*min(DeltaSpectrum(1:round(0.25*size(DeltaSpectrum,2))));

figure(19)
semilogy(frq_shift ,ArterySpectrum) ;
hold on 
semilogy(frq_shift ,BgSpectrum) ;
hold on
semilogy(frq_shift ,DeltaSpectrum) ;

grid on ;
title('Average Spectrogram of the doppler shifts '); % averaged outside of segmented vessels
legend('Average Spectrogram in the arteries ','Average Spectrogram in the background ','difference between the artery and the backgroung','Location','southeast' ) ;
%xticklabels({'-1','wall start', '-0.6','-0.4', '-0.2','0','0.2','0.4', '0.6', '0.8','wall end'});
xlabel('frequence in kHz','FontSize',14) ;
axis tight;
%ylim([mimin mamax]);
ylabel('average density of the spectrum(u.a.)','FontSize',14) ;
hold off

print('-f19','-dpng',fullfile(one_cycle_dir_png,strcat(filename,'_Average spectrogramms.png'))) ;


