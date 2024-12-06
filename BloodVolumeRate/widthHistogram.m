function widthHistogram(area, width_std, name)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);

figure("Visible","off")
histogram(2 * sqrt(area(area ~= 0) / pi) * 1000, 50, FaceColor = 'k');
set(gca, 'Linewidth', 2)

aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);
title(sprintf('Histogram of %s sections width (µm)', name));

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('histogram_of_%s_section_width.png', name))))

writematrix(2 * sqrt(area / pi) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s_section_widths.txt', name))));
writematrix(width_std * PW_params.cropSection_pixelSize / (2 ^ PW_params.k) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_standard_deviation_%s_section_width.txt", ToolBox.main_foldername, name)));

end