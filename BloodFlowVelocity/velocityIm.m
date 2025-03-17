function [f] = velocityIm(v_mean, mask, cmap, name, options)

arguments
    v_mean 
    mask 
    cmap 
    name
    options.colorbarOn 
end

TB = getGlobalToolBox;

f = figure("Visible","off");
colormap(cmap)
imagesc(v_mean .* mask)
if options.colorbarOn
    c = colorbar;
    c.Label.String = 'mm/s';
end
axis image; axis off;

% Export the figure as a PNG with a white background
exportgraphics(gcf, fullfile(TB.path_png, 'bloodFlowVelocity', sprintf("%s_v_mean_%s.png", TB.main_foldername, name)), 'Resolution', 300);

end