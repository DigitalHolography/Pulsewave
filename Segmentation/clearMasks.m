function [mask] = clearMasks(mask, name, colormap, ToolBox)
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
params = ToolBox.getParams;
main_folder = ToolBox.main_foldername;

minPixelSize = params.json.Mask.MinPixelSize;
imCloseRadius = params.json.Mask.ImcloseRadius;
minWidth = params.json.Mask.MinimumVesselWidth;
imDilateSize = params.json.Mask.FinalDilation;

% Ensure the output directory exists
outputDir = fullfile(ToolBox.path_png, 'mask', 'steps');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Step 1: Remove small objects using area opening
mask = bwareaopen(mask, minPixelSize, 4);
imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_1_AreaOpened.png", main_folder, name)));

% Step 2: Close small gaps using morphological closing
imcloseSE = strel('disk', imCloseRadius);
mask = imclose(mask, imcloseSE);
imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_2_Closed.png", main_folder, name)));

% Step 3: Ensure minimum mask width using skeletonization and dilation
minWidthSE = strel('disk', minWidth);
skel = bwskel(mask);
mask = mask | imdilate(skel, minWidthSE);
imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_3_MinWidth.png", main_folder, name)));

% Step 4: Final dilation
dilationSE = strel('disk', imDilateSize);
mask = imdilate(mask, dilationSE);
imwrite(mask, colormap, fullfile(outputDir, sprintf("%s_%s_4_Dilated.png", main_folder, name)));

end
