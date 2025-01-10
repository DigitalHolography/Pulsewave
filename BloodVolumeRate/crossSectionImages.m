function crossSectionImages(M0_ff_img, xy_barycenter, area, bvr, velocity, mask_r, locs, name)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
numCircles = size(area, 1);
numVesselMax = size(area, 2);
exportVideos = PW_params.exportVideos;

section_width_plot = figure("Visible","off");

mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate'), 'sectionsImages')
mkdir(fullfile(ToolBox.PW_path_eps, 'volumeRate'), 'sectionsImages')

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

section_width_plot.Position = [200 200 600 600];
vesselWidthsVideo = zeros(600, 600, 3, numCircles);
vesselNumVideo = zeros(600, 600, 3, numCircles);
vesselBVRVideo = zeros(600, 600, 3, numCircles);
vesselMaxVelocityVideo = zeros(600, 600, 3, numCircles);

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

        for etIdx = 1:size(locs{circleIdx}, 1)
            new_x = x_barycenter + ratio_etiquette * (locs{circleIdx}(etIdx, 2) - x_barycenter);
            new_y = y_barycenter + ratio_etiquette * (locs{circleIdx}(etIdx, 1) - y_barycenter);
            text(new_x, new_y, sprintf(string(etiquettes_frame_values(etIdx))), "FontWeight", "bold", "FontSize", 12, "Color", "white", "BackgroundColor", "black");
        end

    end

    title(sprintf("Cross section width in %s (µm)", titleFig));
    set(gca, 'FontSize', 14)
    vesselWidthsVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(section_width_plot)));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_crossSectionWidth%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_crossSectionWidth%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end
NumVessel = (1:numVesselMax);
if name == "Artery"
    col = [1,0,0];
else
    col = [0,0,1];
end
figure(312)
title(sprintf("Numerotation in %s (µm)", titleFig));
set(gca, 'FontSize', 14)
for circleIdx = 1:numCircles
    etiquettes_frame_values = append(sprintf("%s°",name(1)),string(NumVessel(1:find(area(circleIdx, :),1,'last')))); % last non zero cross section area
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = col);
    vesselNumVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(312)));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_Numerotation%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_Numerotation%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end
figure(312)
title(sprintf("Mean Blood Volume Rate in %s (µm)", titleFig));
set(gca, 'FontSize', 14)
for circleIdx = 1:numCircles
    etiquettes_frame_values = round(mean(bvr(circleIdx,1:find(area(circleIdx, :),1,'last'),:),3),1); 
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = col);
    vesselBVRVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(312)));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_MeanBVR%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_MeanBVR%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end
figure(312)
title(sprintf("Max Velocity in %s (µm)", titleFig));
set(gca, 'FontSize', 14)
for circleIdx = 1:numCircles
    vel = velocity{circleIdx}';
    sz = size(vel{1});
    vel = reshape(cell2mat(vel),[],sz(1),sz(2));
    etiquettes_frame_values = round(mean(max(vel,[],2),3),1); 
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = col);
    vesselMaxVelocityVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(312)));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_MaxVelocity%sImage.png", ToolBox.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsImages', sprintf("%s_circle_%d_MaxVelocity%sImage.eps", ToolBox.main_foldername, circleIdx, name)))
end

widthtable = table;
widthtable.Label = (1:numVesselMax);
for circleIdx = 1:numCircles
    crossSectionWidthArtery = 2 * sqrt(squeeze(area(circleIdx, :)) / pi) * 1000;
    widthtable.(sprintf("circle %d",circleIdx)) = crossSectionWidthArtery;
end
tablefilename = fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_', 'Vessels_Widths_Table', '.txt'));
writetable(widthtable,tablefilename);
if exportVideos
    writeGifOnDisc(vesselWidthsVideo, sprintf('sectionsWidth%s', name), 0.15, 10);
    writeGifOnDisc(vesselNumVideo, sprintf('sectionsWidth%s', name), 0.15, 10);
    writeGifOnDisc(vesselBVRVideo, sprintf('sectionsWidth%s', name), 0.15, 10);
    writeGifOnDisc(vesselMaxVelocityVideo, sprintf('sectionsWidth%s', name), 0.15, 10);

end

close(312);
end