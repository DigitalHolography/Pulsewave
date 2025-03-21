function Aresult = maskedNearestNeighbor(A, d, n, maskzone, maskapply)
% Function to compute the average of n nearest neighbor pixel values
% for every pixel in 'maskapply', considering only the masked region
% 'maskzone' in the 2D image A, using MATLAB's knnsearch.
%
% Inputs:
% A - 2D image (matrix)
% d - Maximum search distance for nearest neighbors
% n - Number of nearest neighbors to consider
% maskzone - Logical mask defining the region to search neighbors within
% maskapply - Logical mask defining the region where neighbors are to be applied
%
% Output:
% Aresult - 2D image containing averaged nearest neighbor values in the 'maskapply' region

% Get the coordinates of the maskzone and maskapply
[maskzone_y, maskzone_x] = find(maskzone);
[maskapply_y, maskapply_x] = find(maskapply);

% Combine the coordinates into arrays
maskzone_coords = [maskzone_y, maskzone_x];
maskapply_coords = [maskapply_y, maskapply_x];

% Initialize the result matrix
Aresult = A;

% Use knnsearch to find the n nearest neighbors for each point in maskapply
[indices, distances] = knnsearch(maskzone_coords, maskapply_coords, 'K', n);

% Loop over all pixels in maskapply
for i = 1:size(maskapply_coords, 1)
    % Get the distances and corresponding indices for the current point
    neighbor_indices = indices(i, :);
    neighbor_distances = distances(i, :);

    % Filter neighbors within the maximum distance d
    valid_neighbors = neighbor_distances <= d;

    % If valid neighbors exist, compute the average value
    if any(valid_neighbors)
        valid_indices = neighbor_indices(valid_neighbors);
        neighbor_values = A(sub2ind(size(A), maskzone_coords(valid_indices, 1), maskzone_coords(valid_indices, 2)));
        Aresult(maskapply_coords(i, 1), maskapply_coords(i, 2)) = mean(neighbor_values);
    else
        % If no valid neighbors, retain the original value
        Aresult(maskapply_coords(i, 1), maskapply_coords(i, 2)) = A(maskapply_coords(i, 1), maskapply_coords(i, 2));
    end

end

% Ensure pixels outside maskapply retain their original values
Aresult(~maskapply) = A(~maskapply);
end
