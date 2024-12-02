function [hue, sat, val, cmap] = API2HSVmap(API, Im, mask)

    tolVal = [0.02, 0.98];

    hue = zeros(size(Im)); %ARImap
    adjusted_image = mat2gray(imadjust(Im, stretchlim(Im, tolVal)));

    sat_max = 1.6;
    sat = API .* mask .* adjusted_image ./ sat_max;
    healthy_API = 0.75;
    sat = sat_max .* (sat - healthy_API) .* mask ./ (1 - healthy_API) ;
    sat(sat<0) = 0;
    val = adjusted_image;

    cmap = ones(256, 3);
    cmap(193:256, 2:3) = [linspace(1, 0, 64)' linspace(1, 0, 64)'];

end
