function graphSignal(ToolBox, filename, folder, x, y, style, color, opt)
% Graph Function for Signal display

arguments
    ToolBox
    filename {mustBeText}
    folder {mustBeText}
end

arguments(Repeating)
    x {mustBeNumeric}
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
    opt.Legends
    opt.TxtName = {}
    opt.TxtFigX
    opt.TxtFigY
    opt.TxtFigString
end

figure(Visible="off")
hold on
for n = 1:length(y)
    plot(x{n}, y{n}, style{n}, 'Color', color{n}, 'LineWidth', opt.LineWidth)

end

if ~isempty(opt.TxtFigX)
    for n = 1:length(opt.TxtFigX)
        text(opt.TxtFigX{n}, opt.TxtFigY, opt.TxtFigString)
    end
end

title(opt.Title)

fontsize(gca, opt.Fontsize, "points");
xlabel(opt.xlabel, 'FontSize', opt.Fontsize);
ylabel(opt.ylabel, 'FontSize', opt.Fontsize);

pbaspect([1.618 1 1]);
set(gca, 'LineWidth', opt.LineWidth);

if ~isempty(opt.Legends)
    legends(opt.Legends)
end

if ~isempty(opt.TxtName)
    for n = 1:length(opt.TxtName)
        plot2txt(x{n}, y{n}, opt.TxtName{n}, ToolBox)
    end
end

% Axis are the minimum value of every X and the maximum value of every X
minX = min(x{1});
maxX = max(x{1});

for n = 2:length(x)
    minX = max(minX, x{n});
    maxX = max(maxX, x{n});
end

axis padded
ax = axis;
axis([minX, maxX, ax(3), ax(4)])
box on

exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s.png", ToolBox.main_foldername, filename)))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s.eps", ToolBox.main_foldername, filename)))

end
