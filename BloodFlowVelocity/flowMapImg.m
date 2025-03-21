function [combined_img] = flowMapImg(img, masks, cmaps, opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

arguments
    img
    masks
    cmaps
    opt.circles = []
    opt.background = []
end

N = size(masks, 2);

if N ~= size(cmaps, 2)
    error('Number of masks and colormaps must be the same')
end

if ~isempty(opt.circles)

    for i = 1:N
        masks{N + i} = masks{i} & opt.circles;
        masks{i} = masks{i} & ~opt.circles;
        cmaps{N + i} = flip(cmaps{i}, 1);
    end

    N = 2 * N;
end

[numX, numY, ~] = size(img);

layers = zeros(numX, numY, 3, N);

for i = 1:size(masks, 2)
    mask = masks{i};
    cmap = cmaps{i};
    layers(:, :, :, i) = setcmap(img, mask, cmap);
end

maskVessel = zeros(numX, numY);

for i = 1:N
    maskVessel = maskVessel + masks{i};
end

if isempty(opt.circles)
    opt.circles = zeros(numX, numY);
end

if isempty(opt.background)
    opt.background = zeros(numX, numY);
end

combined_img = sum(layers, 4) + (opt.background .* ~maskVessel .* ~opt.circles) + opt.circles .* ~maskVessel;

end
