function sectionImage(M0_ff_img, mask_r, name)

ToolBox = getGlobalToolBox;
[numX, numY] = size(M0_ff_img);
numCircles = size(mask_r, 3);

colors = lines(numCircles);
imgRGB = repmat(M0_ff_img, 1, 1, 3);

for circleIdx = 1:numCircles
    mask = mask_r(:, :, circleIdx);
    indxs = find(mask > 0);
    imgRGB(indxs) = colors(circleIdx, 1);
    imgRGB(numY * numX + indxs) = colors(circleIdx, 2);
    imgRGB(2 * numY * numX + indxs) = colors(circleIdx, 3);

    if circleIdx > 1 % intersections should be drawn in white
        previous_mask = mask_r(:, :, circleIdx - 1);
        indxs = find(mask > 0 & previous_mask > 0);
        imgRGB(indxs) = 1;
        imgRGB(numY * numX + indxs) = 1;
        imgRGB(2 * numY * numX + indxs) = 1;
    end

end

figure("Visible", "off")
imshow(imgRGB)
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_%s.png", ToolBox.main_foldername, sprintf('%s_sections%s', name))))

end
