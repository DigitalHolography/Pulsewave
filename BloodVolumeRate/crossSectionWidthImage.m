function crossSectionWidthImage(M0_ff_img, xy_barycenter, area, mask_r, locs, name)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
numCircles = size(area, 1);
exportVideos = PW_params.exportVideos;

section_width_plot = figure("Visible","off");

mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate'), 'sectionsWidth')
mkdir(fullfile(ToolBox.PW_path_eps, 'volumeRate'), 'sectionsWidth')

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

section_width_plot.Position = [200 200 600 600];
vesselWidthsVideo = zeros(600, 600, 3, numCircles);

if strcmp(name, 'Artery') == 1
    color = [1 0 0];
    titleFig = "arteries";
else
    color = [0 0 1];
    titleFig = "veins";
end

ratio_etiquette = 1.2;

for circleIdx = 1:numCircles
    crossSectionWidthArtery = 2 * sqrt(squeeze(area(circleIdx, :)) / pi) * 1000;
    etiquettes_frame_values = append(string(round(crossSectionWidthArtery, 0)), "µm");

    image_RGB = repmat(M0_ff_img - M0_ff_img .* squeeze(mask_r(:, :, circleIdx)), 1, 1, 3) + reshape(color, 1, 1, 3) .* squeeze(mask_r(:, :, circleIdx)) .* M0_ff_img; % adding the Red value to the mask pixels
    imagesc(image_RGB);
    axis image
    axis off

    if ~isempty(locs{circleIdx})

        for etIdx = 1:length(locs{circleIdx})
            new_x = x_barycenter + ratio_etiquette * (locs{circleIdx}(etIdx, 2) - x_barycenter);
            new_y = y_barycenter + ratio_etiquette * (locs{circleIdx}(etIdx, 1) - y_barycenter);
            text(new_x, new_y, sprintf(string(etiquettes_frame_values(etIdx))), "FontWeight", "bold", "FontSize", 12, "Color", "white", "BackgroundColor", "black");
        end

    end

    title(sprintf("Cross section width in %s (µm)", titleFig));
    set(gca, 'FontSize', 14)
    vesselWidthsVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(section_width_plot)));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_crossSectionWidth%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_crossSectionWidth%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end
if exportVideos
    writeGifOnDisc(vesselWidthsVideo, sprintf('sectionsWidth%s', name), 0.1, 10);
end
end