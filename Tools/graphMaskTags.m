function graphMaskTags(figId, Image, mask, etiquettes_locs, etiquettes_values, x_center, y_center, NameValueArgs)
% Plots on an existing fig the image combination of a raw image and a mask displayed in red
arguments
    figId
    Image
    mask
    etiquettes_locs
    etiquettes_values
    x_center
    y_center
    NameValueArgs.Fontsize double = 14
    NameValueArgs.Color (1, 3) double = [1 0 0]
    NameValueArgs.Title = []
    NameValueArgs.Visible = false
    
end

ratio_etiquette = 1.2;
fig = figure(figId);
fig.Position = [200 200 600 600];
fig.Visible = NameValueArgs.Visible;




image_RGB = repmat(Image - Image .* mask, 1, 1, 3) + reshape(NameValueArgs.Color, 1, 1, 3) .* mask .* Image; % adding the Red value to the mask pixels
imagesc(image_RGB);
axis image
axis off
if ~isempty(NameValueArgs.Title)
    title(NameValueArgs.Title);
    set(gca, 'FontSize', 14);
end
if ~isempty(etiquettes_locs)
    
    for etIdx = 1:size(etiquettes_locs,1)
        new_x = x_center + ratio_etiquette * (etiquettes_locs(etIdx, 2) - x_center);
        new_y = y_center + ratio_etiquette * (etiquettes_locs(etIdx, 1) - y_center);
        text(new_x, new_y, sprintf(string(etiquettes_values(etIdx))), "FontWeight", "bold", "FontSize", NameValueArgs.Fontsize, "Color", "white", "BackgroundColor", "black");
    end
    
end

end
