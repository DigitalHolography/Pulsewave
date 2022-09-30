function mask = createVesselMask(video)
mask = squeeze(mean(video, 3));

mask = flat_field_correction(mask, ceil(size(mask,1)*0.07), 0.25);

% mask = log10(abs(mask));
%FIXME : add threshold parameter
mask = imbinarize(mat2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);
end