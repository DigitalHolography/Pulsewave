function [radialAverage, binCenters] = computeRadialAverage(image)
% Compute the radial average of an image.
%
% Inputs:
%   image: 2D image array.
%   centerX: X-coordinate of the center.
%   centerY: Y-coordinate of the center.
%
% Output:
%   radialAverage: 1D array of radial average intensities.

[numY, numX, ~] = size(image);
[X, Y] = meshgrid(1:numX, 1:numY);
R = sqrt((X - numX / 2) .^ 2 + (Y - numY / 2) .^ 2); % Distance from center

maxRadius = max(R(:));
numBins = ceil(maxRadius);
binEdges = linspace(0, maxRadius, numBins + 1);
binCenters = (binEdges(1:end - 1) + binEdges(2:end)) / 2;

radialAverage = zeros(1, numBins);

for binIdx = 1:numBins
    mask = (R >= binEdges(binIdx)) & (R < binEdges(binIdx + 1));
    radialAverage(binIdx) = mean(image(mask), 'omitnan');
end

end
