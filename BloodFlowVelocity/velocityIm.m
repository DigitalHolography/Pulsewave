function [f] = velocityIm(v_mean, mask, cmap, name, options)

arguments
    v_mean 
    mask 
    cmap 
    name
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
exportgraphics(gcf, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_v_imagesc_Mean_%s.png", ToolBox.main_foldername, name)), 'BackgroundColor', 'none', 'ContentType', 'vector', 'Resolution', 300);

end