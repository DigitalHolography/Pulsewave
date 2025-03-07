function crossSectionImages(M0_ff_img, xy_barycenter, area, bvr, velocity, mask_r, locs, name)

TB = getGlobalToolBox;
params = TB.getParams;
numCircles = size(area, 1);
numVesselMax = size(area, 2);
exportVideos = params.exportVideos;

section_width_plot = figure("Visible","off");

mkdir(fullfile(TB.path_png, 'volumeRate'), 'sectionsImages')
mkdir(fullfile(TB.path_eps, 'volumeRate'), 'sectionsImages')
mkdir(fullfile(TB.path_png, 'volumeRate', 'sectionsImages'),  'widths')
mkdir(fullfile(TB.path_eps, 'volumeRate', 'sectionsImages'), 'widths')
mkdir(fullfile(TB.path_png, 'volumeRate', 'sectionsImages'), 'num')
mkdir(fullfile(TB.path_eps, 'volumeRate', 'sectionsImages'), 'num')
mkdir(fullfile(TB.path_png, 'volumeRate', 'sectionsImages'), 'bvr')
mkdir(fullfile(TB.path_eps, 'volumeRate', 'sectionsImages'), 'bvr')
mkdir(fullfile(TB.path_png, 'volumeRate', 'sectionsImages'), 'vel')
mkdir(fullfile(TB.path_eps, 'volumeRate', 'sectionsImages'), 'vel')

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

section_width_plot.Position = [200 200 600 600];
vesselWidthsVideo = zeros(465, 465, 3, numCircles);
vesselNumVideo = zeros(465, 465, 3, numCircles);
vesselBVRVideo = zeros(465, 465, 3, numCircles);
vesselMaxVelocityVideo = zeros(465, 465, 3, numCircles);

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

    image_RGB = repmat(M0_ff_img - M0_ff_img .* mask_r(:, :, circleIdx), 1, 1, 3) + reshape(color, 1, 1, 3) .* mask_r(:, :, circleIdx) .* M0_ff_img; % adding the Red value to the mask pixels
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

    title(sprintf("cross section width in %s (µm)", titleFig));
    set(gca, 'FontSize', 14)
    vesselWidthsVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(gca)));

    exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', 'sectionsImages', 'widths', sprintf("%s_circle_%d_crossSectionWidth%sImage.png", TB.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(TB.path_eps, 'volumeRate', 'sectionsImages', 'widths', sprintf("%s_circle_%d_crossSectionWidth%sImage.eps", TB.main_foldername, circleIdx, name)))
end
NumVessel = (1:numVesselMax);
figure(312)
parfor circleIdx = 1:numCircles
    str = name;
    etiquettes_frame_values = append(sprintf("%s°",str(1)),string(NumVessel(1:find(area(circleIdx, :),1,'last')))); % last non zero cross section area
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = name, Title = sprintf("Numerotation in %s", titleFig));
    vesselNumVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(gca)));

    exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', 'sectionsImages', 'num', sprintf("%s_circle_%d_Numerotation%sImage.png", TB.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(TB.path_eps, 'volumeRate', 'sectionsImages', 'num', sprintf("%s_circle_%d_Numerotation%sImage.eps", TB.main_foldername, circleIdx, name)))
end
figure(312)
parfor circleIdx = 1:numCircles
    etiquettes_frame_values = round(mean(bvr(circleIdx,1:find(area(circleIdx, :),1,'last'),:),3),1); 
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = name, Title = sprintf("Mean Blood Volume Rate in %s (µL/min)", titleFig));
    vesselBVRVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(gca)));

    exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', 'sectionsImages', 'bvr', sprintf("%s_circle_%d_MeanBVR%sImage.png", TB.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(TB.path_eps, 'volumeRate', 'sectionsImages', 'bvr', sprintf("%s_circle_%d_MeanBVR%sImage.eps", TB.main_foldername, circleIdx, name)))
end
figure(312)
parfor circleIdx = 1:numCircles
    vel = velocity{circleIdx}';
    if not(isempty(vel))
        sz = size(vel{1});
        vel = reshape(cell2mat(vel),[],sz(1),sz(2));
    else
        sz = [1 1] ;
        vel = zeros([],sz(1),sz(2));
    end
    
    etiquettes_frame_values = round(mean(max(vel,[],2),3),1); 
    graphMaskTags(312, M0_ff_img, squeeze(mask_r(:, :, circleIdx)), locs{circleIdx}, etiquettes_frame_values, x_barycenter, y_barycenter, Color = name, Title = sprintf("Max Velocity in %s (mm/s)", titleFig));
    vesselMaxVelocityVideo(:, :, :, circleIdx) = rescale(frame2im(getframe(gca)));

    exportgraphics(gca, fullfile(TB.path_png, 'volumeRate', 'sectionsImages', 'vel', sprintf("%s_circle_%d_MaxVelocity%sImage.png", TB.main_foldername, circleIdx, name)))
    exportgraphics(gca, fullfile(TB.path_eps, 'volumeRate', 'sectionsImages', 'vel', sprintf("%s_circle_%d_MaxVelocity%sImage.eps", TB.main_foldername, circleIdx, name)))
end

widthtable = table;
widthtable.Label = (1:numVesselMax);
for circleIdx = 1:numCircles
    crossSectionWidthArtery = 2 * sqrt(squeeze(area(circleIdx, :)) / pi) * 1000;
    widthtable.(sprintf("circle %d",circleIdx)) = crossSectionWidthArtery;
end
tablefilename = fullfile(TB.path_txt, strcat(TB.main_foldername, '_', 'Vessels_Widths_Table', '.txt'));
writetable(widthtable,tablefilename);
if exportVideos
    writeGifOnDisc(vesselWidthsVideo, sprintf('sectionsWidth%s', name), 0.15, 10);
    writeGifOnDisc(vesselNumVideo, sprintf('Numerotation%s', name), 0.15, 10);
    writeGifOnDisc(vesselBVRVideo, sprintf('MeanBVR%s', name), 0.15, 10);
    writeGifOnDisc(vesselMaxVelocityVideo, sprintf('MaxVelocity%s', name), 0.15, 10);
end

close(312);
end