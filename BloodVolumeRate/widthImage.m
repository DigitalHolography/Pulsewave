function widthImage(subimage_size, sub_images_r, numSections, name)

ToolBox = getGlobalToolBox;

numCircles = size(sub_images_r, 2);
figure("Visible", "off")
numSectionMax = max(numSections);
% fill with zero images the zeros parts

sub_images_mat = zeros(subimage_size(1), subimage_size(2), numSectionMax * numCircles);

for circleIdx = 1:numCircles

    for sectionIdx = 1:numSectionMax

        if size(sub_images_r{1, circleIdx}, 2) < numSectionMax
            sub_images_r{1, circleIdx}(end + 1) = {zeros(subimage_size, 'single')};
        end

        sub_images_mat(:, :, ((sectionIdx - 1) * numCircles) + circleIdx) = sub_images_r{1, circleIdx}{1, sectionIdx};

    end

end

montage(sub_images_mat, "Size", [max(1, numSectionMax), max(1, numCircles)]);
exportgraphics(gca, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('5_all_%s_sections_with_increasing_radius.png', name))))

end
