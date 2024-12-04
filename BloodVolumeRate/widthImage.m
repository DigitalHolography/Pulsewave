function widthImage(sub_images_r,  name)

numCircles = size(sub_images_r, 1);
ToolBox = getGlobalToolBox;

figure("Visible","off")
numSectionMax = size(sub_images_r, 1);
% fill with zero images the zeros parts
subimage_size = size(sub_images_r{1,1},1);

for circleIdx = 1:numCircles

    for sectionIdx = 1:numSectionMax

        if isempty(sub_images_r{circleIdx, sectionIdx})
            sub_images_r{circleIdx, sectionIdx} = zeros(subimage_size, 'single');
        end

    end

end

montage(sub_images_r(1:numCircles,1:numSectionMax),"Size",[numSectionMax, numCircles])
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('5_all_%s_sections_with_increasing_radius.png', name))))

end