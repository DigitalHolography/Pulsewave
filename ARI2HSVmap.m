function [hue, sat, val, cmap] = ARI2HSVmap(ARI, Im, mask, ToolBox)

    tolVal = [0.02, 0.98];

    hue = zeros(size(Im)); %ARImap
    adjusted_image = mat2gray(imadjust(Im, stretchlim(Im, tolVal)));

    sat = ARI .* mask .* adjusted_image;
    healthy_ARI = 0.75;
    sat = (sat - healthy_ARI) / (1 - healthy_ARI);
    sat(sat <= 0) = 0;

    val = adjusted_image;

    cmap = ones(256, 3);
    cmap(193:256, 2:3) = [linspace(1, 0, 64)' linspace(1, 0, 64)'];

end
