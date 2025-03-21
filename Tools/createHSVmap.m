function [hue, sat, val, cmap] = createHSVmap(Im, mask, hueValMin, hueValMax)

%Im grayscale image with value between 0 and 1
[numX, numY] = size(Im);
Ones = ones(numX, numY);

Im_hue = mat2gray(Im .* mask + Ones .* (~mask) .* hueValMin);
hue = (Im_hue * (hueValMax - hueValMin) + hueValMin);

sat = Ones .* mask;

val = sigmoid(Im, 0.3, 5) .* mask + Im .* ~mask;

ARI_x = linspace(0, 1, 256);
ARI_h = (ARI_x * (hueValMax - hueValMin) + hueValMin);
ARI_s = ones(1, 256);
ARI_v = rescale(sigmoid(ARI_x, 0.3, 5));

cmap = hsv2rgb(ARI_h, ARI_s, ARI_v);
cmap = squeeze(cmap);
end
