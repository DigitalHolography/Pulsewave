function [maskCRA,maskCRV] = createCentralRetinalArteryVeinMask(video)

mask = squeeze(mean(video, 3));
% mask = flat_field_correction(mask, ceil(size(mask,1)*0.07), 0.25);

% mask = log10(abs(mask));
%FIXME : add threshold parameter
maskCRA = mask>(5.0*std2(mask));
maskCRV = mask<(-5.0*std2(mask));

% mask thresholds : make a mix between min/max of data and std 

end