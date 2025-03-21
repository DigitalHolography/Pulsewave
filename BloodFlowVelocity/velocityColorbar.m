function [colorfig, hCB] = velocityColorbar(cmap, v_min, v_max, name)

arguments
    cmap
    v_min
    v_max
    name
end

ToolBox = getGlobalToolBox;

% Save colorbar
colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(cmap)
hCB = colorbar('north', 'Ticks', [0, 1], 'TickLabels', {string(round(v_min, 1)), string(round(v_max, 1))});
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
hCB.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");

exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_colorbarVelocity%s.png", ToolBox.main_foldername, name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, 'bloodFlowVelocity', sprintf("%s_colorbarVelocity%s.eps", ToolBox.main_foldername, name)))
end
