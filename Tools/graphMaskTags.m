function graphMaskTags(figId, Image, mask, etiquettes_locs, etiquettes_values, x_center, y_center, opt)
% Plots on an existing fig the image combination of a raw image and a mask displayed in red
arguments
    figId
    Image
    mask
    etiquettes_locs
    etiquettes_values
    x_center
    y_center
    opt.Fontsize double = 14
    opt.Color = 'none'
    opt.Title = []
    opt.Visible = false
    opt.circles = [];

end

ratio_etiquette = 1.2;

fig = figure(figId);
fig.Position = [200 200 600 600];
fig.Visible = opt.Visible;

if strcmp(opt.Color, 'Artery') == 1
    cmap = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
elseif strcmp(opt.Color, 'Vein') == 1
    cmap = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
else
    cmap = cmapLAB(256, [0 0 0], 0, [1 1 1], 1);
end

image_RGB = setcmap(Image, mask, cmap) + Image .* ~mask;

if ~isempty(opt.circles)
    image_RGB = image_RGB .* ~opt.circles + opt.circles .* 0.7;
end

% image_RGB = repmat(Image - Image .* mask, 1, 1, 3) + reshape(NameValueArgs.Color, 1, 1, 3) .* mask .* Image; % adding the Red value to the mask pixels
imagesc(image_RGB);
axis image
axis off

if ~isempty(opt.Title)
    title(opt.Title);
    set(gca, 'FontSize', 14);
end

% Get the image dimensions
[imgHeight, imgWidth, ~] = size(Image);

if ~isempty(etiquettes_locs)

    for etIdx = 1:size(etiquettes_locs, 1)

        try
            % Calculate the new position for the text
            new_x = x_center + ratio_etiquette * (etiquettes_locs(etIdx, 2) - x_center);
            new_y = y_center + ratio_etiquette * (etiquettes_locs(etIdx, 1) - y_center);

            % Ensure the text is within the image bounds
            new_x = max(1, min(new_x, imgWidth)); % Clamp x to image width
            new_y = max(1, min(new_y, imgHeight)); % Clamp y to image height

            % Add the text
            text(new_x, new_y, sprintf(string(etiquettes_values(etIdx))), ...
                "FontWeight", "bold", ...
                "FontSize", opt.Fontsize, ...
                "Color", "white", ...
                "BackgroundColor", "black");
        catch
        end

    end

end
