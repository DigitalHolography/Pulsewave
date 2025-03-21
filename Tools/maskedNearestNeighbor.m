function Aresult = maskedNearestNeighbor(A, d, maskzone, maskapply)
% Function to compute nearest neighbor pixel values of every pixel in
% 'maskapply', considering only the masked region 'maskzone' in the 2D image A.
%
% Inputs:
% A - 2D image (matrix)
% d - Maximum search distance for nearest neighbor
% maskzone - Logical mask defining the region to search neighbors within
% maskapply - Logical mask defining the region where neighbors are to be applied
%
% Output:
% Aresult - 2D image containing nearest neighbor values in the 'maskapply' region

% Initialize the result matrix
Aresult = A;

% Get indices of maskzone and maskapply
[maskzone_y, maskzone_x] = find(maskzone);
[maskapply_y, maskapply_x] = find(maskapply);

% Loop over each pixel in maskapply
for i = 1:length(maskapply_y)
    y_apply = maskapply_y(i);
    x_apply = maskapply_x(i);

    % Initialize the minimum distance and nearest value
    min_dist = inf;
    nearest_val = NaN;

    % Loop over each pixel in maskzone
    for j = 1:length(maskzone_y)
        y_zone = maskzone_y(j);
        x_zone = maskzone_x(j);

        % Calculate the Euclidean distance
        dist = sqrt((y_zone - y_apply) ^ 2 + (x_zone - x_apply) ^ 2);

        % Check if within maximum distance and update the nearest neighbor
        if dist <= d && dist < min_dist
            min_dist = dist;
            nearest_val = A(y_zone, x_zone);
        end

    end

    % Assign the nearest value to the result image
    Aresult(y_apply, x_apply) = nearest_val;
end

% Ensure that pixels outside maskapply retain their original values
Aresult(~maskapply) = A(~maskapply);
end
