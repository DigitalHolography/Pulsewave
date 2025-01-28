function [mask] = clearMasks(ToolBox, mask, name, colormap)
    % clearMasks - Processes a binary mask to remove small objects, close gaps, and dilate.
    %
    % Inputs:
    %   ToolBox - Structure containing paths and parameters.
    %   mask    - Binary mask to be processed.
    %   name    - Name identifier for saving intermediate results.
    %
    % Output:
    %   maskDilated - The final processed mask after dilation.

    % Validate inputs
    if nargin < 3
        error('Not enough input arguments. Expected ToolBox, mask, and name.');
    end

    if ~islogical(mask)
        error('Mask must be a logical array.');
    end

    if ~ischar(name) && ~isstring(name)
        error('Name must be a character array or string.');
    end

    % Load parameters
    PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
    main_folder = ToolBox.main_foldername;

    % Ensure the output directory exists
    outputDir = fullfile(ToolBox.PW_path_png, 'mask', 'steps');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Step 1: Remove small objects using area opening
    mask = bwareaopen(mask, PW_params.masks_minSize, 4);
    imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_1_AreaOpened.png", main_folder, name)));

    % Step 2: Fill
    % mask = imfill(mask, "holes");
    % imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_2_Filled.png", main_folder, name)));

    % Step 3: Close small gaps using morphological closing
    imcloseSE = strel('disk', PW_params.masks_imclose_radius);
    mask = imclose(mask, imcloseSE);
    imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_3_Closed.png", main_folder, name)));

    % Step 4: Ensure minimum mask width using skeletonization and dilation
    minWidthSE = strel('disk', PW_params.masks_min_width);
    skel = bwskel(mask);
    mask = mask | imdilate(skel, minWidthSE);
    imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_4_MinWidth.png", main_folder, name)));

    % Step 5: Final dilation
    dilationSE = strel('disk', PW_params.masks_imdilateFinal);
    mask = imdilate(mask, dilationSE);
    imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_5_Dilated.png", main_folder, name)));

end