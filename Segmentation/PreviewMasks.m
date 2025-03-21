function PreviewMasks(app)

if ~app.flag_is_load
    app.LoadfolderButtonPushed();

    if ~app.flag_is_load
        error("Please select a correct HoloDoppler folder")
    end

end

params = Parameters_json(app.file.directory, app.file.param_name);

% Create a uifigure for mask preview with dark gray background
d = uifigure('Position', [300, 300, 1000, 750], ...
    'Color', [0.2, 0.2, 0.2], ... % Dark gray background for the figure
    'Name', 'Mask Preview Tool', ...
    'Resize', 'on', ...
    'WindowStyle', 'normal');

% Create a GridLayout for better control over resizing
layout = uigridlayout(d, [2, 3]); % Two rows, two columns
layout.BackgroundColor = [0.2, 0.2, 0.2];
layout.RowHeight = {60, '1x'}; % Set the first row for the button and second row to fill the remaining space
layout.ColumnWidth = {210, '1x', '1x'}; % The first column is smaller, second column is wider for the image and panel

% Create "Update Parameters" button in the top-right grid
updateParamsButton = uibutton(layout, ...
    'FontName', 'Helvetica', ...
    'BackgroundColor', [0.5, 0.5, 0.5], ...
    'FontColor', [0.9 0.9 0.9], ... % Use TextColor instead of ForegroundColor
    'FontWeight', 'bold', ...
    'Text', 'Update Parameters', ...
    'ButtonPushedFcn', @UpdateParameters);

% Position the "Update Parameters" button inside grid row 1, column 2
updateParamsButton.Layout.Row = 1;
updateParamsButton.Layout.Column = 1;

% Create "Create Mask" button using uibutton for layout compatibility
createMaskButton = uibutton(layout, ...
    'FontName', 'Helvetica', ...
    'BackgroundColor', [0.5, 0.5, 0.5], ...
    'FontColor', [0.9 0.9 0.9], ... % Use TextColor instead of ForegroundColor
    'FontWeight', 'bold', ...
    'Text', 'Create Mask', ...
    'ButtonPushedFcn', @createMask);

% Position the button inside grid row 1, column 1
createMaskButton.Layout.Row = 1;
createMaskButton.Layout.Column = 2;

% Create the panel for additional UI elements (numeric edit fields, etc.) with dark gray background
paramPanel = uipanel(layout, ...
    'ForegroundColor', [0.8 0.8 0.8], ...
    'Title', 'Parameters', ...
    'BackgroundColor', [0.2, 0.2, 0.2]); % Dark gray background for the panel

% Position the panel inside grid row 2, column 1
paramPanel.Layout.Row = 2;
paramPanel.Layout.Column = 1;

histoGrid = uigridlayout(layout, [3 1], ...
    'BackgroundColor', [0.2, 0.2, 0.2]);

histoGrid.Layout.Row = 2;
histoGrid.Layout.Column = 3;

% Create "Open Folder" button
openFolderButton = uibutton(layout, ...
    'FontName', 'Helvetica', ...
    'BackgroundColor', [0.5, 0.5, 0.5], ...
    'FontColor', [0.9 0.9 0.9], ...
    'FontWeight', 'bold', ...
    'Text', 'Open Folder', ...
    'ButtonPushedFcn', @openFolder);

% Position the button inside grid row 1, column 3
openFolderButton.Layout.Row = 1;
openFolderButton.Layout.Column = 3;

% Create numeric edit fields for each parameter

% Mask Thresholds
uicontrol(paramPanel, 'Style', 'text', 'String', 'Vascular Threshold:', 'Position', [10, 520, 120, 20]);
vascularThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.VascularThreshold, 'Position', [140, 520, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Arterial Threshold:', 'Position', [10, 400, 120, 20]);
arterialThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.ArterialThreshold, 'Position', [140, 400, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Venous Threshold:', 'Position', [10, 320, 120, 20]);
venousThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.VenousThreshold, 'Position', [140, 320, 30, 22]);

diasysCheckbox = uicheckbox(paramPanel, 'Text', 'Diastole-Systole Analysis', 'Position', [10, 440, 170, 20], ...
    'FontColor', [0.9 0.9 0.9], 'Value', params.json.Mask.DiaSysAnalysis);

% Geometric Parameters
uicontrol(paramPanel, 'Style', 'text', 'String', 'Diaphragm Radius:', 'Position', [10, 600, 120, 20]);
diaphragmRadiusEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.DiaphragmRadius, 'Position', [140, 600, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Coroid Radius:', 'Position', [10, 560, 120, 20]);
cropCoroidEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.CropChoroidRadius, 'Position', [140, 560, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Min Area Size:', 'Position', [10, 200, 120, 20]);
minSizeEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.MinPixelSize, 'Position', [140, 200, 30, 22]);

% Image Processing Parameters
uicontrol(paramPanel, 'Style', 'text', 'String', 'Gaussian Filter Size:', 'Position', [10, 240, 120, 20]);
gaussianFilterEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.Blur, 'Position', [140, 240, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Imclose Radius:', 'Position', [10, 160, 120, 20]);
imcloseRadiusEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.ImcloseRadius, 'Position', [140, 160, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Min Width:', 'Position', [10, 120, 120, 20]);
minWidthEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.MinimumVesselWidth, 'Position', [140, 120, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Imdilate Final:', 'Position', [10, 80, 120, 20]);
imdilateFinalEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.FinalDilation, 'Position', [140, 80, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Force Width:', 'Position', [10, 40, 120, 20]);
forceVesselWidthEdit = uieditfield(paramPanel, 'numeric', 'Value', params.json.Mask.ForceVesselWidth, 'Position', [140, 40, 30, 22]);

% Create a table for the VascularClasses, ArterialClasses, VenousClasses

vascularClasses = sprintf('%g,', params.json.Mask.VascularClasses);
arterialClasses = sprintf('%g,', params.json.Mask.ArterialClasses);
venousClasses = sprintf('%g,', params.json.Mask.VenousClasses);

% Create table for VascularClasses
uicontrol(paramPanel, 'Style', 'text', 'String', 'Vascular Classes:', 'Position', [10, 480, 120, 20]);
vascularClassesTable = uieditfield(paramPanel, 'Value', vascularClasses(1:end - 1), 'Position', [140, 480, 60, 22]);

% Create table for ArterialClasses
uicontrol(paramPanel, 'Style', 'text', 'String', 'Arterial Classes:', 'Position', [10, 360, 120, 20]);
arterialClassesTable = uieditfield(paramPanel, 'Value', arterialClasses(1:end - 1), 'Position', [140, 360, 60, 22]);

% Create table for VenousClasses
uicontrol(paramPanel, 'Style', 'text', 'String', 'Venous Classes:', 'Position', [10, 280, 120, 20]);
venousClassesTable = uieditfield(paramPanel, 'Value', venousClasses(1:end - 1), 'Position', [140, 280, 60, 22]);

% Create the axes in the dialog to display the image with dark gray background
ax = axes(layout, 'Color', [0.2, 0.2, 0.2]); % Set dark gray background for the axes

% Position the axes inside grid row 1, column 2
ax.Layout.Row = 2;
ax.Layout.Column = 2;

% Create the axes in the dialog to display the image with dark gray background
ax0 = axes(histoGrid, 'Color', [0.2, 0.2, 0.2]); % Set dark gray background for the axes

% Position the axes inside grid row 1, column 2
ax0.Layout.Row = 1;
ax0.Layout.Column = 1;

% Create the axes in the dialog to display the image with dark gray background
ax1 = axes(histoGrid, 'Color', [0.2, 0.2, 0.2]); % Set dark gray background for the axes

% Position the axes inside grid row 1, column 2
ax1.Layout.Row = 2;
ax1.Layout.Column = 1;

% Create the axes in the dialog to display the image with dark gray background
ax2 = axes(histoGrid, 'Color', [0.2, 0.2, 0.2]); % Set dark gray background for the axes

% Position the axes inside grid row 1, column 2
ax2.Layout.Row = 3;
ax2.Layout.Column = 1;

% Wait for user interaction with the dialog
uiwait(d);

% Nested function for mask creation
function createMask(~, ~)

    ToolBox = ToolBoxClass(app.file.directory, app.file.param_name, true);

    % Ensure necessary steps are done before proceeding
    if ~app.file.is_preprocessed
        app.file.preprocessData();
    end

    % Update parameters and get the image and masks
    M0_ff_img = squeeze(mean(app.file.M0_ff_video, 3)); % Mean of the M0_ff_video
    f_AVG_mean = squeeze(mean(app.file.f_AVG_video, 3)); % Mean of the f_AVG_video

    % Generate the masks
    [maskArtery, maskVein] = createMasks(app.file.M0_ff_video, f_AVG_mean);

    % Prepare RGB image to display
    cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
    cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
    cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

    M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
    M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
    M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

    M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);

    % Display the image with the masks on the axes in the dialog
    imshow(M0_RGB, 'Parent', ax);
    imshow(imread(fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_all_16_Histo.png", ToolBox.main_foldername))), 'Parent', ax0);

    if isfile(fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_artery_23_Histo.png", ToolBox.main_foldername)))
        imshow(imread(fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_artery_23_Histo.png", ToolBox.main_foldername))), 'Parent', ax1);
    end

    if isfile(fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_vein_23_Histo.png", ToolBox.main_foldername)))
        imshow(imread(fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_vein_23_Histo.png", ToolBox.main_foldername))), 'Parent', ax2);
    end

end

function UpdateParameters(~, ~)

    % Extract parameter values from the UI components
    params.json.Mask.DiaphragmRadius = diaphragmRadiusEdit.Value;
    params.json.Mask.Blur = gaussianFilterEdit.Value;
    params.json.Mask.CropChoroidRadius = cropCoroidEdit.Value;
    params.json.Mask.VascularCorrelationMapThreshold = vascularThresholdEdit.Value;
    params.json.Mask.ArterialCorrelationMapThreshold = arterialThresholdEdit.Value;
    params.json.Mask.VenousCorrelationMapThreshold = venousThresholdEdit.Value;
    params.json.Mask.MinPixelSize = minSizeEdit.Value;
    params.json.Mask.ImcloseRadius = imcloseRadiusEdit.Value;
    params.json.Mask.MinimumVesselWidth = minWidthEdit.Value;
    params.json.Mask.FinalDilation = imdilateFinalEdit.Value;
    params.json.Mask.DiaSysAnalysis = diasysCheckbox.Value;
    params.json.Mask.ForceVesselWidth = forceVesselWidthEdit.Value;

    % Get the data from tables
    params.json.Mask.VascularClasses = str2double(split(vascularClassesTable.Value, ','));
    params.json.Mask.ArterialClasses = str2double(split(arterialClassesTable.Value, ','));
    params.json.Mask.VenousClasses = str2double(split(venousClassesTable.Value, ','));

    params.WriteParametersToJson(fullfile(app.file.directory, 'eyeflow', 'json', app.file.param_name));
end

function openFolder(~, ~)

    ToolBox = getGlobalToolBox;

    if isempty(ToolBox)
        fprintf(2, "You must create the masks first\n")
        return
    end

    % Specify the folder to open
    folderPath = fullfile(ToolBox.path_png, 'mask', 'steps');

    % Check if the folder exists
    if isfolder(folderPath)
        winopen(folderPath); % Open folder in the system's file explorer
    else
        uialert(d, 'The folder does not exist!', 'Error', 'Icon', 'error');
    end

end

% Close the dialog after the user interaction
delete(d); % You might want to remove or comment this line if the user needs to manually close the dialog
end
