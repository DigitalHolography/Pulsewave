function corrected_image = flat_field_correction(image, gaussianParams, borderAmount)
    % FLAT_FIELD_CORRECTION Corrects an image for uneven illumination using a fitted Gaussian.
    %
    % Inputs:
    %   image: Input image (2D array).
    %   gaussianParams: Gaussian parameters [A, mu, sigma, C].
    %   borderAmount: Fraction of the image border to exclude (optional, default = 0).
    %
    % Output:
    %   corrected_image: Flat-field corrected image.

    % Set default borderAmount if not provided
    if nargin < 3
        borderAmount = 0;
    end

    flag = 0;
    % Check if image is normalized between 0 and 1
    if min(image, [], 'all') < 0 || max(image, [], 'all') > 1
        Im_max = max(image, [], 'all');
        Im_min = min(image, [], 'all');
        image = (image - Im_min) ./ (Im_max - Im_min);
        flag = 1;
    end

    % Define the non-border region
    if borderAmount == 0
        a = 1;
        b = size(image, 1);
        c = 1;
        d = size(image, 2);
    else
        a = ceil(size(image, 1) * borderAmount);
        b = floor(size(image, 1) * (1 - borderAmount));
        c = ceil(size(image, 2) * borderAmount);
        d = floor(size(image, 2) * (1 - borderAmount));
    end

    % Compute the sum of intensities in the non-border region
    ms = sum(image(a:b, c:d), [1 2]);

    % Compute the Gaussian illumination pattern
    [numX, numY, ~] = size(image);
    [X, Y] = meshgrid(1:numY, 1:numX);
    R = sqrt((X - numX/2).^2 + (Y - numY/2).^2); % Distance from the center
    illumination_pattern = gaussianParams(1) * exp(-((R) / gaussianParams(3)).^2) + gaussianParams(4);

    % Normalize the image using the Gaussian illumination pattern
    image = image ./ illumination_pattern;

    % Rescale to maintain total intensity in the non-border region
    ms2 = sum(image(a:b, c:d), [1 2]);
    corrected_image = (ms / ms2) .* image;

    % Rescale back to original intensity range if normalization was applied
    if flag == 1
        corrected_image = Im_min + (Im_max - Im_min) .* corrected_image;
    end
end