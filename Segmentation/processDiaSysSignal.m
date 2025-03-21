function [mask1, mask2] = processDiaSysSignal(diasys, maskClean, params, cmap, suffix)

if params.threshold >= -1 && params.threshold <= 1
    % Manual threshold
    diasys = rescale(diasys);
    mask1 = diasys >= params.threshold;
    mask2 = diasys <= params.threshold;
    graphThreshHistogram(diasys, params.threshold, maskClean, cmap, suffix);
else
    % Automatic Otsu thresholding
    [mask1, mask2] = autoOtsuThresholding(diasys, maskClean, params.classes, suffix);
end

mask1 = logical(mask1);
mask2 = logical(mask2);

end
