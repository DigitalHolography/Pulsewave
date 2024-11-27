function [outIm] = setcmap(Im, mask, cmap)

%Im grayscale image with value between 0 and 1

% Normalize the grayscale image data to the range [0, 1]
normalizedImage = rescale(Im);

% Map the normalized grayscale data to the colormap
% The grayscale values are mapped to the colormap index
[~, colormap_indices] = min(abs(cmap(:, 1) - normalizedImage(:)'), [], 1); % Finding closest colormap index

% Create the RGB image
mappedIm = reshape(cmap(colormap_indices, :), size(Im, 1), size(Im, 2), 3);
outIm = mappedIm .* mask + Im .* ~mask;

end
