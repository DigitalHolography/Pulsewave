function crossSectionWidthImage(M0_ff_img, xy_barycenter, area_r, numSections, width_avg_r, locs, name)

ToolBox = getGlobalToolBox;
numCircles = size(area_r, 2);

section_width_plot = figure("Visible","off");

mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate'), 'sectionsWidth')
mkdir(fullfile(ToolBox.PW_path_eps, 'volumeRate'), 'sectionsWidth')

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

section_width_plot.Position = [200 200 600 600];
[~, ~, width, height] = gca.Position;
vesselWidthsVideo = zeros(width, height, 3, numCircles);

for circleIdx = 1:numCircles
    crossSectionWidthArtery = 2 * sqrt(area_r(1:numSections(circleIdx), :) / pi) * 1000;
    etiquettes_frame_values = append(string(round(crossSectionWidthArtery, 1)), "µm");
    graphMaskTags(section_width_plot, M0_ff_img, squeeze(width_avg_r{circleIdx}(:, :)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Fontsize = 12);
    title(sprintf("%s", 'Cross section width in arteries (µm)'));
    set(gca, 'FontSize', 14)
    vesselWidthsVideo(:, :, :, circleIdx) = frame2im(getframe(section_width_plot));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_crossSectionWidth%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_crossSectionWidth%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end

writeGifOnDisc(vesselWidthsVideo, 'sectionsWidth.gif', 0.1);
end