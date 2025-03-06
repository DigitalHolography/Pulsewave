function widthHistogram(width, width_std, area, name)

TB = getGlobalToolBox;

figure("Visible","off")
histogram(2 * sqrt(area(area ~= 0) / pi) * 1000, 50, FaceColor = 'k');
set(gca, 'Linewidth', 2)

aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);
title(sprintf('Histogram of %s sections width (µm)', name));

exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', sprintf("%s_%s", TB.main_foldername, sprintf('histogram_of_%s_section_width.png', name))))

%csv output of the widths
T = table();
numR = length(width); % number of radii

for rIdx = 1:numR
    numSection = length(width{rIdx});
    for sectionIdx = 1:numSection
        T.(sprintf('Width_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(width{rIdx}(sectionIdx)));
        T.(sprintf('STD_Width_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(width_std{rIdx}(sectionIdx)));
    end
end

writetable(T,fullfile(TB.path_txt, strcat(TB.main_foldername, '_', 'WidthTable', '_', name, '.csv')));

end