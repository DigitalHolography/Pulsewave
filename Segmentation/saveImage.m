function saveImage(I, ToolBox, suffix, opt)

arguments
    I
    ToolBox
    suffix
    opt.cmap = []
    opt.isStep = false
end

main_folder = ToolBox.main_foldername;

if opt.isStep
    folderPath = fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, suffix));
else
    folderPath = fullfile(ToolBox.path_png, 'mask', sprintf("%s_%s", main_folder, suffix));
end

if isempty(opt.cmap)

    if size(I, 3) == 1 && isa(I, 'numeric')
        I = rescale(I);
    end

    imwrite(I, folderPath, 'png');
else

    if size(I, 3) == 1 && ~islogical(I), I = (size(opt.cmap, 1) - 1) * rescale(I);
    end

    imwrite(I, opt.cmap, folderPath, 'png');
end

end
