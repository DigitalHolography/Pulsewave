function [ArterySpectrum, BgSpectrum, DeltaSpectrum, sizeSHFreq, Lorenz_Arteries, Lorenz_BKG, M0_artery, M0_background, M1M0_artery, M1M0_background, M2M0_artery, M2M0_background] = createSpectrum(maskArtery, maskBackground, SH)

ToolBox = getGlobalToolBox;
SH_arteries = SH .* maskArtery;
SH_background = SH .* maskBackground;
fs = ToolBox.fs;
f1 = ToolBox.f1;
f2 = ToolBox.f2;
[numY, numX, ~] = size(SH);

sizeSHFreq = size(SH, 3);
ArterySpectrum = zeros(1, sizeSHFreq);
BgSpectrum = zeros(1, sizeSHFreq);

ArterySpectrumM1M0 = zeros(numY, numX, sizeSHFreq);
BgSpectrumM1M0 = zeros(numY, numX, sizeSHFreq);
ArterySpectrumM2M0 = zeros(numY, numX, sizeSHFreq);
BgSpectrumM2M0 = zeros(numY, numX, sizeSHFreq);

Ones = ones(numY, numX);

%% integration interval
% convert frequencies to indices
n1 = ceil(f1 * sizeSHFreq / fs);
n2 = ceil(f2 * sizeSHFreq / fs);

% symetric integration interval
n3 = sizeSHFreq - n2 + 1;
n4 = sizeSHFreq - n1 + 1;

f_range = (n1:n2) .* (fs / sizeSHFreq);
f_range_sym = (-n2:-n1) .* (fs / sizeSHFreq);

% FIX ME artery-background < 0
%% Calcul spectrum M0
for ii = 1:size(SH_arteries, 3)
    ArterySpectrum(ii) = sum(SH_arteries(:, :, ii), 'all') / nnz(SH_arteries(:, :, ii));
    BgSpectrum(ii) = sum(SH_background(:, :, ii), 'all') / nnz(SH_background(:, :, ii));

    % moment = gather(squeeze(sum(abs(SH(:, :, n1:n2)), 3))) + gather(squeeze(sum(abs(SH(:, :, n3:n4)), 3)));
    M0_artery = gather(squeeze(sum(SH_arteries(:, :, n1:n2), 3))) + gather(squeeze(sum(SH_arteries(:, :, n3:n4), 3)));
    M0_background = gather(squeeze(sum(SH_background(:, :, n1:n2), 3))) + gather(squeeze(sum(SH_background(:, :, n3:n4), 3)));

end

%% Calcul spectrum M1/M0

ArterySpectrumM1M0(:, :, n1:n2) = SH_arteries(:, :, n1:n2) .* reshape(f_range, 1, 1, numel(f_range));
ArterySpectrumM1M0(:, :, n3:n4) = SH_arteries(:, :, n3:n4) .* reshape(f_range_sym, 1, 1, numel(f_range_sym));

BgSpectrumM1M0(:, :, n1:n2) = SH_background(:, :, n1:n2) .* reshape(f_range, 1, 1, numel(f_range));
BgSpectrumM1M0(:, :, n3:n4) = SH_background(:, :, n3:n4) .* reshape(f_range_sym, 1, 1, numel(f_range_sym));

% moment = gather(squeeze(sum(abs(SH(:, :, n1:n2)), 3))) + gather(squeeze(sum(abs(SH(:, :, n3:n4)), 3)));
M1_artery = gather(squeeze(sum(ArterySpectrumM1M0(:, :, n1:n2), 3))) + gather(squeeze(sum(ArterySpectrumM1M0(:, :, n3:n4), 3)));
M1_backgroung = gather(squeeze(sum(BgSpectrumM1M0(:, :, n1:n2), 3))) + gather(squeeze(sum(BgSpectrumM1M0(:, :, n3:n4), 3)));

M1M0_artery = (M1_artery ./ (M0_artery + Ones .* not(maskArtery)));
M1M0_background = (M1_backgroung ./ (M0_background + Ones .* not(maskBackground)));

%% Calcul spectrum M2/M0

ArterySpectrumM2M0(:, :, n1:n2) = SH_arteries(:, :, n1:n2) .* (reshape(f_range, 1, 1, numel(f_range)) .^ 2);
ArterySpectrumM2M0(:, :, n3:n4) = SH_arteries(:, :, n3:n4) .* (reshape(f_range_sym, 1, 1, numel(f_range_sym)) .^ 2);

BgSpectrumM2M0(:, :, n1:n2) = SH_background(:, :, n1:n2) .* (reshape(f_range, 1, 1, numel(f_range)) .^ 2);
BgSpectrumM2M0(:, :, n3:n4) = SH_background(:, :, n3:n4) .* (reshape(f_range_sym, 1, 1, numel(f_range_sym)) .^ 2);

% moment = gather(squeeze(sum(abs(SH(:, :, n1:n2)), 3))) + gather(squeeze(sum(abs(SH(:, :, n3:n4)), 3)));
M2_artery = gather(squeeze(sum(ArterySpectrumM2M0(:, :, n1:n2), 3))) + gather(squeeze(sum(ArterySpectrumM2M0(:, :, n3:n4), 3)));
M2_backgroung = gather(squeeze(sum(BgSpectrumM2M0(:, :, n1:n2), 3))) + gather(squeeze(sum(BgSpectrumM2M0(:, :, n3:n4), 3)));

M2M0_artery = sqrt(M2_artery ./ (M0_artery + Ones .* not(maskArtery)));
M2M0_background = sqrt(M2_backgroung ./ (M0_background + Ones .* not(maskBackground)));

q = sizeSHFreq / 2;

Lorenz_BKG = zeros(1, q);
Lorenz_Arteries = zeros(1, q);
% freq = linspace(1, fs / 2 - fs / sizeSHFreq, sizeSHFreq / 2);
Lorenz_Arteries(1, :) = (ArterySpectrum(1:q) + flip(ArterySpectrum((q + 1):sizeSHFreq))) / 2;
Lorenz_BKG(1, :) = (BgSpectrum(1:q) + flip(BgSpectrum((q + 1):sizeSHFreq))) / 2;

%Shift for displayS
ArterySpectrum = circshift(ArterySpectrum, round(size(ArterySpectrum) / 2));
BgSpectrum = circshift(BgSpectrum, round(size(BgSpectrum) / 2));

DeltaSpectrum = ArterySpectrum - BgSpectrum;
DeltaSpectrum(DeltaSpectrum < 1) = 1;
end
