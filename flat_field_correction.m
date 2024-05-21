function corrected_image = flat_field_correction(image, gw, borderAmount)

    % this function should only be used just before display,
    % with [0 1]-valued mat2gray(array) arg.
    % avoid to use it for zero-mean arrays

    % typical blur parameter : gw = 0.07*size(image,1)

    flag = 0;
    % Check if image is normalized between 0 and 1
    if min(image, [], 'all') < 0 || max(image, [], 'all') > 1
        Im_max = max(image, [], 'all');
        Im_min = min(image, [], 'all');
        image = (image - Im_min) ./ (Im_max - Im_min);
        flag = 1;
    end

    if (borderAmount == 0)
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

    ms = sum(image(a:b, c:d), [1 2]);
    image = image ./ imgaussfilt(image, gw);
    ms2 = sum(image(a:b, c:d), [1 2]);
    corrected_image = (ms / ms2) .* image;

    if flag == 1
        corrected_image = Im_min + (Im_max - Im_min) .* corrected_image;
    end

end
