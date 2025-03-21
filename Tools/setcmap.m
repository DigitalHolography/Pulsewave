function [outIm] = setcmap(Im, mask, cmap)
% Im: grayscale image with values between 0 and 1
% mask: binary mask to apply to the output image
% cmap: colormap to apply to the grayscale image

% Number of levels in the colormap
num_levels = size(cmap, 1);

Im = rescale(Im);

% Normalize the grayscale image to the range [1, num_levels]
% This avoids the need for the linspace and min operations
Im_normalized = round(Im * (num_levels - 1)) + 1;

% Ensure that the indices are within the valid range
Im_normalized = max(1, min(Im_normalized, num_levels));

% Map the normalized grayscale data to the colormap
mappedIm = cmap(Im_normalized, :);

% Reshape the mapped image to the original image size
mappedIm = reshape(mappedIm, [size(Im, 1), size(Im, 2), 3]);

% Apply the mask to the output image
outIm = mappedIm .* mask;
end
