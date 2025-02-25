function [rotatedImg, orientation] = rotateSubImage(subImg)
% Rotate the sub-image to align the blood vessel vertically.
% The orientation is taken from the centroid closest to the center of the image.
%
% Input:
%   subImg - 2D array, the sub-image containing the blood vessel.
%
% Output:
%   rotatedImg - 2D array, the rotated image.
%   orientation - Scalar, the orientation angle used for rotation.

%% UNCOMMENT IN CASE
% angles = linspace(0, 180, 181);
% projx = zeros(size(subImg, 1), length(angles));
% projy = zeros(size(subImg, 2), length(angles));
%
% for theta = 1:length(angles)
%     tmpImg = imrotate(subImg, angles(theta), 'bilinear', 'crop');
%     projx(:, theta) = squeeze(sum(tmpImg, 1));
%     projy(:, theta) = squeeze(sum(tmpImg, 2));
% end
%
% projx_bin = (projx == 0);
% list_x = squeeze(sum(projx_bin, 1));
% [~, idc] = max(list_x);
% tilt_angle = idc(1);
% rotatedImg = imrotate(subImg, tilt_angle, 'bilinear', 'crop');
% Convert the image to double for processing
subImg = im2double(subImg);

% Threshold the image to isolate the blood vessel
bw = imbinarize(subImg); % Binarize the image
bw = bwareaopen(bw, 50); % Remove small objects (noise)

% Compute region properties (centroid and orientation)
stats = regionprops(bw, 'Centroid', 'Orientation');

% Check if any regions were detected
if isempty(stats)
    rotatedImg = subImg;
    orientation = 0;
    return
end

% Get the center of the image
[rows, cols] = size(subImg);
imageCenter = [cols / 2, rows / 2];

% Calculate the distance of each centroid to the image center
centroids = cat(1, stats.Centroid); % Concatenate centroids into a matrix
distances = sqrt((centroids(:, 1) - imageCenter(1)).^2 + ...
    (centroids(:, 2) - imageCenter(2)).^2);

% Find the index of the centroid closest to the center
[~, closestIdx] = min(distances);

% Get the orientation of the closest centroid
orientation = stats(closestIdx).Orientation;

% Adjust orientation to make the blood vessel vertical
orientation = orientation + 90; % Add 90 degrees to align vertically

% Rotate the image to make the blood vessel vertical
rotatedImg = imrotate(subImg, -orientation, 'bilinear', 'crop');
end