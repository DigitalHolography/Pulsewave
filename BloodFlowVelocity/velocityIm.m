function [f] = velocityIm(v_mean, mask, cmap, options)

arguments
    v_mean 
    mask 
    cmap 
    options.colorbarOn 
end

ToolBox = getGlobalToolBox;

f = figure("Visible","off");
colormap(cmap)
imagesc(v_mean .* mask)
if options.colorbarOn
    c = colorbar;
end
c.Label.String = 'mm/s';
axis image; axis off;
exportgraphics(gcf, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, "v_imagesc_Mean_arteries.png")), 'BackgroundColor', 'none', 'ContentType', 'vector', 'Resolution', 300);

end