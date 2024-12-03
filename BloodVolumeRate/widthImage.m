function widthImage(sub_images_r, numSections, name)

numCircles = size(sub_images_r, 1);
ToolBox = getGlobalToolBox;

figure("Visible","off")
numSectionMax = max(numSections);
% fill with zero images the zeros parts
subimage_cells = sub_images_r{1};
subimage_size = size(subimage_cells{1}, 1);

for circleIdx = 1:numCircles

    for j = 1:numSectionMax

        if isempty(sub_images_r{circleIdx}{j})
            sub_images_r{circleIdx}{j} = zeros(subimage_size, 'single');
        end

    end

end

montage(cell2mat(sub_images_r), "Size", [numSectionMax, numCircles])
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('5_all_%s_sections_with_increasing_radius.png', name))))

end