function img = cropCircle(subImg)
subImgHW = min([size(subImg, 1) size(subImg, 2)]);
radius = round(subImgHW / 2);
%FIXME anamorph.
center = [round(size(subImg, 1) / 2) round(size(subImg, 2) / 2)];
[xx, yy] = meshgrid(1:size(subImg, 1), 1:size(subImg, 2));
circular_mask = false(size(subImg, 1), size(subImg, 2));
circular_mask = circular_mask | hypot(xx - center(1), yy - center(2)) <= radius;
img = double(circular_mask) .* subImg;
end
