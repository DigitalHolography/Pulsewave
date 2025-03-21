function crossSectionImages(M0_ff_img, xy_barycenter, area, Q, v, mask, locs, name)

% Get global toolbox and parameters
ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;
path_txt = ToolBox.path_txt;
params = ToolBox.getParams;
main_folder = ToolBox.main_foldername;
exportVideos = params.exportVideos;

% Determine color and title based on vessel type
if strcmp(name, 'Artery')
    color = [1, 0, 0]; % Red for arteries
    titleFig = "arteries";
else
    color = [0, 0, 1]; % Blue for veins
    titleFig = "veins";
end

% Extract barycenter coordinates
x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

% Initialize video arrays
numCircles = size(area, 1);
vesselD_video = zeros(465, 465, 3, numCircles);
vesselNum_video = zeros(465, 465, 3, numCircles);
Q_video = zeros(465, 465, 3, numCircles);
velocity_video = zeros(465, 465, 3, numCircles);

% Plot cross-section widths
for cIdx = 1:numCircles
    % Calculate cross-section widths in micrometers
    crossSectionWidth = 2 * sqrt(squeeze(area(cIdx, :)) / pi) * 1000;
    etiquettes_frame_values = append(string(round(crossSectionWidth, 0)), "µm");

    % Create RGB image with mask overlay
    image_RGB = repmat(M0_ff_img - M0_ff_img .* mask(:, :, cIdx), 1, 1, 3) + ...
        reshape(color, 1, 1, 3) .* mask(:, :, cIdx) .* M0_ff_img;

    % Plot the image
    figure("Visible", "off");
    imagesc(image_RGB);
    axis image;
    axis off;

    % Add labels for cross-section widths
    if ~isempty(locs{cIdx})
        ratio_etiquette = 1.2;

        for etIdx = 1:size(locs{cIdx}, 1)
            new_x = x_barycenter + ratio_etiquette * (locs{cIdx}(etIdx, 2) - x_barycenter);
            new_y = y_barycenter + ratio_etiquette * (locs{cIdx}(etIdx, 1) - y_barycenter);
            text(new_x, new_y, etiquettes_frame_values(etIdx), ...
                "FontWeight", "bold", "FontSize", 12, "Color", "white", "BackgroundColor", "black");
        end

    end

    % Add title and format plot
    title(sprintf("Cross-section width in %s (µm)", titleFig));
    set(gca, 'FontSize', 12);

    % Capture and resize frame for video
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    vesselD_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'sectionsImages', 'widths', ...
        sprintf("%s_circle_%d_crossSectionWidth%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'volumeRate', 'sectionsImages', 'widths', ...
        sprintf("%s_circle_%d_crossSectionWidth%sImage.eps", main_folder, cIdx, name)));
end

% Plot vessel numeration
numVesselMax = size(area, 2);
NumVessel = (1:numVesselMax);
Initial = name(1);
figure(312)

parfor cIdx = 1:numCircles
    etiquettes_frame_values = append(sprintf("%s°", Initial), string(NumVessel(1:find(area(cIdx, :), 1, 'last')))); % Last non-zero cross-section area
    graphMaskTags(312, M0_ff_img, squeeze(mask(:, :, cIdx)), locs{cIdx}, etiquettes_frame_values, ...
        x_barycenter, y_barycenter, Color = name, Title = sprintf("Numerotation in %s", titleFig));

    % Capture and resize frame for video
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    vesselNum_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'sectionsImages', 'num', ...
        sprintf("%s_circle_%d_Numerotation%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'volumeRate', 'sectionsImages', 'num', ...
        sprintf("%s_circle_%d_Numerotation%sImage.eps", main_folder, cIdx, name)));
end

% Plot mean blood volume rate (Q)
figure(312)

parfor cIdx = 1:numCircles
    etiquettes_frame_values = round(mean(Q(cIdx, 1:find(area(cIdx, :), 1, 'last'), :), 3), 1);
    graphMaskTags(312, M0_ff_img, squeeze(mask(:, :, cIdx)), locs{cIdx}, etiquettes_frame_values, ...
        x_barycenter, y_barycenter, Color = name, Title = sprintf("Average Blood Volume Rate in %s (µL/min)", titleFig));

    % Capture and resize frame for video
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    Q_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'sectionsImages', 'bvr', ...
        sprintf("%s_circle_%d_MeanBVR%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'volumeRate', 'sectionsImages', 'bvr', ...
        sprintf("%s_circle_%d_MeanBVR%sImage.eps", main_folder, cIdx, name)));
end

% Plot maximum velocity
figure(312)

parfor cIdx = 1:numCircles
    vel = v{cIdx};

    if ~isempty(vel)
        etiquettes_frame_values = round(mean(max(vel, [], 2), 3), 1);
    else
        etiquettes_frame_values = zeros(1, 0); % Handle empty velocity case
    end

    graphMaskTags(312, M0_ff_img, squeeze(mask(:, :, cIdx)), locs{cIdx}, etiquettes_frame_values, ...
        x_barycenter, y_barycenter, Color = name, Title = sprintf("Average Velocity in %s (mm/s)", titleFig));

    % Capture and resize frame for video
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    velocity_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'volumeRate', 'sectionsImages', 'vel', ...
        sprintf("%s_circle_%d_MaxVelocity%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'volumeRate', 'sectionsImages', 'vel', ...
        sprintf("%s_circle_%d_MaxVelocity%sImage.eps", main_folder, cIdx, name)));
end

% Create and save table of vessel widths
widthtable = table;
widthtable.Label = (1:numVesselMax)';

for cIdx = 1:numCircles
    crossSectionWidth = 2 * sqrt(squeeze(area(cIdx, :)) / pi) * 1000;
    widthtable.(sprintf("circle_%d", cIdx)) = crossSectionWidth';
end

writetable(widthtable, fullfile(path_txt, strcat(main_folder, '_Vessels_Widths_Table.txt')));

% Export videos if enabled
if exportVideos
    writeGifOnDisc(vesselD_video, sprintf('sectionsWidth%s', name), 0.15, 10);
    writeGifOnDisc(vesselNum_video, sprintf('Numerotation%s', name), 0.15, 10);
    writeGifOnDisc(Q_video, sprintf('MeanBVR%s', name), 0.15, 10);
    writeGifOnDisc(velocity_video, sprintf('MaxVelocity%s', name), 0.15, 10);
end

close(312);

end
