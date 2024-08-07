function corrected_image = flat_field_correction(image, gw, borderAmount)

    % this function should only be used just before display,
    % with [0 1]-valued mat2gray(array) arg.
    % avoid to use it for zero-mean arrays

    % typical blur parameter : gw = 0.07*size(image,1)

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

end
