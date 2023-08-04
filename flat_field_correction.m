function corrected_image = flat_field_correction(image, gw, borderAmount)
% Apply flat field correction to an image before display.
% Only use this function with [0 1]-valued mat2gray(array) input.
% Avoid using it for zero-mean arrays.

% Check input validity
if nargin < 3
    borderAmount = 0;
    if nargin < 2
        gw = 0.07*size(image,1); % typical blur parameter
    end
end

flag = 0;
% Check if image is normalized between 0 and 1
if min(image,[],'all') < 0 || max(image,[],'all') > 1
    Im_max = max(image,[],'all');
    Im_min = min(image,[],'all');
    image = (image - Im_min)./(Im_max - Im_min);
    flag = 1;
end

% Validate gw parameter
if gw <= 0
    error('Parameter gw must be positive.');
end

% Calculate border indices
if (borderAmount == 0)
    a = 1;
    b = size(image,1);
    c = 1;
    d = size(image,2);
else
    a = ceil(size(image,1)*borderAmount);
    b = floor(size(image,1)*(1-borderAmount));
    c = ceil(size(image,2)*borderAmount);
    d = floor(size(image,2)*(1-borderAmount));
end

% Calculate sums of pixel values in the region of interest
ms = sum(image(a:b,c:d), [1 2]);

% Apply flat field correction
blurred_image = imgaussfilt(image, gw);
blurred_image(blurred_image < 1e-3) = 1e-3;
normalized_image = image ./ blurred_image;
ms2 = sum(normalized_image(a:b,c:d), [1 2]);
correction_factor = ms ./ ms2;
corrected_image = correction_factor .* image;

if flag == 1
    corrected_image = Im_min + (Im_max-Im_min).*corrected_image;
end

end


