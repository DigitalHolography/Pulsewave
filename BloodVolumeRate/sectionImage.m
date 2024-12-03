function sectionImage(numCircles, M0_ff_img, width_avg_r, name)

ToolBox = getGlobalToolBox;
[numX, numY] = size(M0_ff_img);

colors = lines(numCircles);
imgRGB = repmat(M0_ff_img, 1, 1, 3);

for circleIdx = 1:numCircles
    width_avg = width_avg_r{circleIdx};
    indxs = find(width_avg(:, :) > 0);
    imgRGB(indxs) = colors(circleIdx, 1);
    imgRGB(numY * numX + indxs) = colors(circleIdx, 2);
    imgRGB(2 * numY * numX + indxs) = colors(circleIdx, 3);

    if circleIdx > 1 % intersections should be drawn in white
        previous_width = width_avg_r{circleIdx - 1};
        indxs = find(width_avg(:, :) > 0 & previous_width(:, :) > 0);
        imgRGB(indxs) = 1;
        imgRGB(numY * numX + indxs) = 1;
        imgRGB(2 * numY * numX + indxs) = 1;
    end

end

figure("Visible","off")
imshow(imgRGB)
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('4_%s_sections.png', name))))

end