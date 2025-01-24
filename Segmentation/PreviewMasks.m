function PreviewMasks(app)

if ~app.flag_is_load
    app.LoadfolderButtonPushed();
    if ~app.flag_is_load
        error("Please select a correct HoloDoppler folder")
    end
end

ToolBox = getGlobalToolBox;
checkPulsewaveParamsFromJson(app.file.directory);
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);

% Create a uifigure for mask preview with dark gray background
d = uifigure('Position', [300, 300, 1000, 750],...
    'Color', [0.2, 0.2, 0.2],...   % Dark gray background for the figure
    'Name', 'Mask Preview Tool',...
    'Resize', 'on',...
    'WindowStyle', 'normal');

% Create a GridLayout for better control over resizing
layout = uigridlayout(d, [2, 2]); % Two rows, two columns
layout.BackgroundColor = [0.2, 0.2, 0.2];
layout.RowHeight = {60, '1x'};     % Set the first row for the button and second row to fill the remaining space
layout.ColumnWidth = {'1x', '1x' '3x'}; % The first column is smaller, second column is wider for the image and panel

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
updateParamsButton.Layout.Column = 2;

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
createMaskButton.Layout.Column = 3;

% Create the panel for additional UI elements (numeric edit fields, etc.) with dark gray background
paramPanel = uipanel(layout, ...
    'ForegroundColor', [0.8 0.8 0.8], ...
    'Title', 'Parameters', ...
    'BackgroundColor', [0.2, 0.2, 0.2]); % Dark gray background for the panel

% Position the panel inside grid row 2, column 1
paramPanel.Layout.Row = 2;
paramPanel.Layout.Column = 1;

% Create the panel for additional UI elements (numeric edit fields, etc.) with dark gray background
classPanel = uipanel(layout, ...
    'ForegroundColor', [0.8 0.8 0.8], ...
    'Title', 'Classes', ...
    'BackgroundColor', [0.2, 0.2, 0.2]); % Dark gray background for the panel

% Position the panel inside grid row 2, column 1
classPanel.Layout.Row = 2;
classPanel.Layout.Column = 2;

% Create numeric edit fields for each parameter

% Vesselness Parameters
uicontrol(paramPanel, 'Style', 'text', 'String', 'Vesselness Sigma:', 'Position', [10, 320, 120, 20]);
vesselnessSigmaEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_vesselness_sigma, 'Position', [140, 320, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Vesselness Beta:', 'Position', [10, 280, 120, 20]);
vesselnessBetaEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_vesselness_beta, 'Position', [140, 280, 30, 22]);

% Mask Thresholds
uicontrol(paramPanel, 'Style', 'text', 'String', 'Vascular Threshold:', 'Position', [10, 240, 120, 20]);
vascularThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_vascular_threshold, 'Position', [140, 240, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Arterial Threshold:', 'Position', [10, 200, 120, 20]);
arterialThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_arterial_threshold, 'Position', [140, 200, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Venous Threshold:', 'Position', [10, 160, 120, 20]);
venousThresholdEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_venous_threshold, 'Position', [140, 160, 30, 22]);

% Geometric Parameters
uicontrol(paramPanel, 'Style', 'text', 'String', 'Diaphragm Radius:', 'Position', [10, 120, 120, 20]);
diaphragmRadiusEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_diaphragmRadius, 'Position', [140, 120, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Coroid Radius:', 'Position', [10, 80, 120, 20]);
cropCoroidEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_crop_radius, 'Position', [140, 80, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Min Area Size:', 'Position', [10, 40, 120, 20]);
minSizeEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_minSize, 'Position', [140, 40, 30, 22]);

% Image Processing Parameters
uicontrol(paramPanel, 'Style', 'text', 'String', 'Gaussian Filter Size:', 'Position', [10, -40, 120, 20]);
gaussianFilterEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.gauss_filt_size_for_barycenter, 'Position', [140, -40, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Imclose Radius:', 'Position', [10, -80, 120, 20]);
imcloseRadiusEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_imclose_radius, 'Position', [140, -80, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Min Width:', 'Position', [10, -120, 120, 20]);
minWidthEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_min_width, 'Position', [140, -120, 30, 22]);

uicontrol(paramPanel, 'Style', 'text', 'String', 'Imdilate Final:', 'Position', [10, -160, 120, 20]);
imdilateFinalEdit = uieditfield(paramPanel, 'numeric', 'Value', PW_params.masks_imdilateFinal, 'Position', [140, -160, 30, 22]);

% Create a table for the VascularClasses, ArterialClasses, VenousClasses
% Create table for VascularClasses
uicontrol(classPanel, 'Style', 'text', 'String', 'Vascular Classes:', 'Position', [10, 340, 120, 20]);
vascularClassesTable = uitable(classPanel, 'Position', [10, 300, 150, 100], 'Data', PW_params.masks_vascular_classes, ColumnEditable=true);

% Create table for ArterialClasses
uicontrol(classPanel, 'Style', 'text', 'String', 'Arterial Classes:', 'Position', [10, 180, 120, 20]);
arterialClassesTable = uitable(classPanel, 'Position', [10, 140, 150, 100], 'Data', PW_params.masks_arterial_classes, ColumnEditable=true);

% Create table for VenousClasses
uicontrol(classPanel, 'Style', 'text', 'String', 'Venous Classes:', 'Position', [10, 60, 120, 20]);
venousClassesTable = uitable(classPanel, 'Position', [10, 20, 150, 100], 'Data', PW_params.masks_venous_classes, ColumnEditable=true);

% Create the axes in the dialog to display the image with dark gray background
ax = axes(layout, 'Color', [0.2, 0.2, 0.2]);  % Set dark gray background for the axes

% Position the axes inside grid row 1, column 2
ax.Layout.Row = 2;
ax.Layout.Column = 3;

% Wait for user interaction with the dialog
uiwait(d);

% Nested function for mask creation
    function createMask(~, ~)
        % Ensure necessary steps are done before proceeding
        if ~app.file.is_preprocessed
            app.PreProcessButtonPushed();
        end

        % Update parameters and get the image and masks
        M0_ff_img = squeeze(mean(app.file.M0_ff_video, 3)); % Mean of the M0_ff_video
        f_AVG_mean = squeeze(mean(app.file.f_AVG_video, 3)); % Mean of the f_AVG_video

        % Generate the masks
        [maskArtery, maskVein] = createMasks(app.file.M0_ff_video, f_AVG_mean);

        % Prepare RGB image to display
        RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArtery; % Red channel for artery mask
        RGBM0(:, :, 2) = rescale(M0_ff_img);              % Green channel for base image
        RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVein;   % Blue channel for vein mask

        % Display the image with the masks on the axes in the dialog
        imshow(RGBM0, 'Parent', ax);
    end

    function UpdateParameters(~, ~)

        % Extract parameter values from the UI components
        PW_params.params.CreationOfMasks.DiaphragmRadius = diaphragmRadiusEdit.Value;
        PW_params.params.CreationOfMasks.VesselnesParameterSigma = vesselnessSigmaEdit.Value;
        PW_params.params.CreationOfMasks.VesselnesParameterBeta = vesselnessBetaEdit.Value;
        PW_params.params.CreationOfMasks.GaussianFilterSizeForBarycenter = gaussianFilterEdit.Value;
        PW_params.params.CreationOfMasks.CropCoroidRadius = cropCoroidEdit.Value;
        PW_params.params.CreationOfMasks.VascularCorrelationMapThreshold = vascularThresholdEdit.Value;
        PW_params.params.CreationOfMasks.ArterialCorrelationMapThreshold = arterialThresholdEdit.Value;
        PW_params.params.CreationOfMasks.VenousCorrelationMapThreshold = venousThresholdEdit.Value;
        PW_params.params.CreationOfMasks.MinimumSeedAreaSize = minSizeEdit.Value;
        PW_params.params.CreationOfMasks.ImcloseRadius = imcloseRadiusEdit.Value;
        PW_params.params.CreationOfMasks.MinimumVesselWidth = minWidthEdit.Value;
        PW_params.params.CreationOfMasks.FinalDilation = imdilateFinalEdit.Value;

        % Get the data from tables
        PW_params.params.CreationOfMasks.VascularClasses = vascularClassesTable.Data;
        PW_params.params.CreationOfMasks.ArterialClasses = arterialClassesTable.Data;
        PW_params.params.CreationOfMasks.VenousClasses = venousClassesTable.Data;

        PW_params.WriteParametersToJson(fullfile(ToolBox.PW_path_main, 'json', 'InputPulsewaveParams.json'));
    end

% Close the dialog after the user interaction
delete(d); % You might want to remove or comment this line if the user needs to manually close the dialog
end
