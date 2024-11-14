function graphSignal(ToolBox, filename, folder, x, y, style, color, opt)
% Graph Function for Signal display

arguments
    ToolBox
    filename {mustBeText}
    folder {mustBeText}
    x {mustBeNumeric}
end

arguments(Repeating)
    y {mustBeNumeric}
    style
    color
end

arguments
    opt.xlabel {mustBeText}  = 'Frame'
    opt.ylabel {mustBeText}  = 'Signal (a.u.)'
    opt.Title {mustBeText}  = 'Signal'
    opt.LineWidth = 2
    opt.Fontsize = 14
end

figure(Visible="off") 
hold on
for n = 1:length(y)
    plot(x, y{n}, style{n}, 'Color', color{n}, 'LineWidth', opt.LineWidth)
end

title(opt.Title)

fontsize(gca, opt.Fontsize, "points");
xlabel(opt.xlabel, 'FontSize', opt.Fontsize);
ylabel(opt.ylabel, 'FontSize', opt.Fontsize);

pbaspect([1.618 1 1]);
set(gca, 'LineWidth', opt.LineWidth);

axis padded
ax = axis;
axis([x(1), x(end), ax(3), ax(4)])
box on

exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s.png", ToolBox.main_foldername, filename)))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s.eps", ToolBox.main_foldername, filename)))
end
