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
    folderPath = fullfile(ToolBox.PW_path_png, 'mask', 'steps', sprintf("%s_%s", main_folder, suffix));
else
    folderPath = fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", main_folder, suffix));
end


if isempty(opt.cmap)
    imwrite(I, folderPath, 'png');
else
    imwrite(I, opt.cmap, folderPath, 'png');
end
end