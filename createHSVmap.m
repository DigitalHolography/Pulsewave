function [hue, sat, val, cmap] = createHSVmap(Im, mask, hueValMin, hueValMax)

    %Im grayscale image with value between 0 and 1
    [numX, numY] = size(Im);
    Im = rescale(Im);
    Ones = ones(numX, numY);
    min_hue = min(Im .* mask + Ones .* (~mask), [], 'all');
    Im_hue = mat2gray(Im .* mask + Ones .* (~mask) .* min_hue);
    hue = (Im_hue * (hueValMax - hueValMin) + hueValMin) .* mask;
    %sat = (1.0-0.5*Im).*mask;
    sat = Ones .* mask;
    val = Im;
    tolVal = [0.01, 0.99];
    val = mat2gray(imadjust(val, stretchlim(val, tolVal)));

    % inflexion_point_hue = mean(Im.*mask,'all');
    % slope_hue = 3;

    %hue = ((hue_val_max-hue_val_min)./(1+exp(-(Im-inflexion_point_hue).*slope_hue)) + hue_val_min).*mask;

    ARI_x = linspace(0, 1, 256);
    ARI_h = (ARI_x * (hueValMax - hueValMin) + hueValMin);
    ARI_s = ones(1, 256);
    ARI_v = (1 - 0.7) * sigmoid(ARI_x, 0.3, 5) + 0.7;

    cmap = squeeze(hsv2rgb(ARI_h, ARI_s, ARI_v));

end

%sat = (1.0 - abs(0.5 * mat2gray(v))) .* double(or(maskArtery,maskVein)) ;
