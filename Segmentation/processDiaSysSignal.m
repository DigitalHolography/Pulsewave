function [mask] = processDiaSysSignal(diasys, maskClean, formerMask, params, cmap, suffix)

diasys = rescale(diasys);

if params.threshold >= -1 && params.threshold <= 1
    % Manual threshold
    mask = diasys >= params.threshold;
    graphThreshHistogram(diasys, params.threshold, maskClean, cmap, [suffix '_2_3']);
else
    % Automatic Otsu thresholding
    mask = autoOtsuThresholding(diasys, maskClean, params.classes, [suffix '_2_3']);
    mask = mask | formerMask; % Include the initial clean mask
end

mask = logical(mask);

end