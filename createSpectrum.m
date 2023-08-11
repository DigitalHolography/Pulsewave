function [ArterySpectrum, BgSpectrum, DeltaSpectrum, sizeSHFreq] = createSpectrum(maskArtery,maskBackground ,SH)

SH_arteries = SH.*maskArtery;
SH_background = SH.*maskBackground;

sizeSHFreq=size(SH_arteries,3);
ArterySpectrum = zeros(1,sizeSHFreq);
BgSpectrum = zeros(1,sizeSHFreq);

% FIX ME artery-background < 0
% Calcul spectrum
for ii = 1:size(SH_arteries,3)
    ArterySpectrum(ii) = sum(SH_arteries(:,:,ii),'all')/nnz(SH_arteries(:,:,ii));
    BgSpectrum(ii) = sum(SH_background(:,:,ii),'all')/nnz(SH_background(:,:,ii));
end

%Shift for display
ArterySpectrum = circshift(ArterySpectrum,round(size(ArterySpectrum)/2));
BgSpectrum = circshift(BgSpectrum,round(size(BgSpectrum)/2));

DeltaSpectrum = ArterySpectrum - BgSpectrum;
DeltaSpectrum(DeltaSpectrum<1)=1;
end