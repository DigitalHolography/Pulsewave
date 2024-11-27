function [hue, sat, val, cmap] = ARI2HSVmap(ARI, Im, mask, ToolBox)

    hue = zeros(size(Im)); %ARImap

    sat = ARI .* mask;
    healthy_ARI = 0.75;
    sat = (sat - healthy_ARI) / (1 - healthy_ARI);
    sat(sat <= 0) = 0;

    val = sigmoid(Im, 0.3, 5) .* mask + Im .* ~mask;

    cmap = ones(256, 3);
    cmap(193:256, 2:3) = [linspace(1, 0, 64)' linspace(1, 0, 64)'];

end
