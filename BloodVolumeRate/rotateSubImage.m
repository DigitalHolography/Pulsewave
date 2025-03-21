function [rotatedImg, orientation] = rotateSubImage(subImg)
% Rotate the sub-image to align the blood vessel vertically.
% The orientation is determined by maximizing the number of zero pixels in the horizontal projection.
%
% Input:
%   subImg - 2D array, the sub-image containing the blood vessel.
%
% Output:
%   rotatedImg - 2D array, the rotated image.
%   orientation - Scalar, the orientation angle used for rotation.

% Define a range of angles to test (0 to 180 degrees in 1-degree increments)
angles = linspace(0, 180, 181);

% Initialize arrays to store horizontal and vertical projections for each angle
projx = zeros(size(subImg, 1), length(angles)); % Horizontal projections

% Loop over each angle and compute projections
for theta = 1:length(angles)
    % Rotate the image by the current angle
    tmpImg = imrotate(subImg, angles(theta), 'bilinear', 'crop');

    % Compute the projection
    projx(:, theta) = squeeze(sum(tmpImg, 1));
end

% Create a binary mask where pixels with zero values in the horizontal projection are set to 1
projx_bin = (projx == 0);

% Sum the binary mask along the rows to count the number of zero pixels for each angle
list_x = squeeze(sum(projx_bin, 1));

% Find the angle with the maximum number of zero pixels in the horizontal projection
[~, idc] = max(list_x);

% Get the corresponding orientation angle
orientation = idc(1);

% Rotate the original image by the determined orientation angle
rotatedImg = imrotate(subImg, orientation, 'bilinear', 'crop');

end
