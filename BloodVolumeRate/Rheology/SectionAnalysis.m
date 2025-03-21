function [avg_blood_rate, avg_blood_velocity, size_section, new_mask] = SectionAnalysis(mask, k)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

skel = bwskel(logical(mask));

[row, col] = find(skel);

x_center = round(mean(row));
y_center = round(mean(col));

size_rect = 2 * max((max(row) - min(row)), (max(col) - min(col)));
rect = [y_center - size_rect x_center - size_rect 2 * size_rect 2 * size_rect];
H_cropped = imcrop(mask, rect);

theta = linspace(0, 180, 181);
projx = zeros(size(H_cropped, 1), length(theta));
projy = zeros(size(H_cropped, 2), length(theta));

for tt = 1:length(theta)
    tmpImg = imrotate(H_cropped, theta(tt), 'bilinear', 'crop');
    projx(:, tt) = squeeze(sum(tmpImg, 1));
    projy(:, tt) = squeeze(sum(tmpImg, 2));
end

figure(3001)
imagesc(projx)

figure(3002)
imagesc(projy)

[~, tilt_angle] = max(sum((projx == 0), 1));
H_cropped = imrotate(H_cropped, tilt_angle, 'crop');
H_cropped = H_cropped .* (H_cropped > 0);

AVG_blood_rate = squeeze(sum(H_cropped, 1)) ./ squeeze(sum(H_cropped > 0, 1) + ~sum(H_cropped > 0, 1));

size_section = nnz(AVG_blood_rate > 0);
section_area = pi * (((size_section) / 2) * (0.7 * params.json.BloodVolumeRate.PixelSize / (2))) ^ 2;

avg_blood_velocity = sum(AVG_blood_rate) / nnz(AVG_blood_rate);
avg_blood_rate = avg_blood_velocity * section_area * 60;

x_center_bis = round(size(H_cropped, 1) / 2);
y_center_bis = round(size(H_cropped, 1) / 2);
size_rect = round(size(H_cropped, 1) / 4);
rect = [y_center_bis - size_rect x_center_bis - size_rect 2 * size_rect 2 * size_rect];
H_cropped = imcrop(H_cropped, rect);
H_cropped = imrotate(H_cropped, -tilt_angle);
size_rect = ceil(size(H_cropped, 1) / 2);
mask = zeros(size(mask, 1), size(mask, 2));
mask((x_center - size_rect + 1):(x_center - size_rect + size(H_cropped, 1)), (y_center - size_rect + 1):(y_center - size_rect + size(H_cropped, 2))) = H_cropped;
new_mask = mask > 0;
figure(k)
imagesc(H_cropped)
end
