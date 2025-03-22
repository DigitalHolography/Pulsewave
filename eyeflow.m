classdef eyeflow < matlab.apps.AppBase

% Properties that correspond to app components
properties (Access = public)
    EyeFlowUIFigure matlab.ui.Figure

    % Load
    LoadFolderButton matlab.ui.control.Button
    LoadHoloButton matlab.ui.control.Button
    ClearButton matlab.ui.control.Button
    FolderManagementButton matlab.ui.control.Button
    ReferenceDirectory matlab.ui.control.TextArea
    statusLamp matlab.ui.control.Lamp
    EditMasksButton matlab.ui.control.Button
    EditParametersButton matlab.ui.control.Button
    OpenDirectoryButton matlab.ui.control.Button
    MaskToolButton matlab.ui.control.Button
    PlayMomentsButton matlab.ui.control.Button

    % Checkboxes
    segmentationCheckBox matlab.ui.control.CheckBox
    bloodFlowAnalysisCheckBox matlab.ui.control.CheckBox
    bloodFlowVelocityFigCheckBox matlab.ui.control.CheckBox
    bloodVolumeRateCheckBox matlab.ui.control.CheckBox
    bloodVolumeRateFigCheckBox matlab.ui.control.CheckBox
    spectralAnalysisCheckBox matlab.ui.control.CheckBox

    % Execute
    NumberofWorkersSpinner matlab.ui.control.Spinner
    NumberofWorkersSpinnerLabel matlab.ui.control.Label
    OverWriteCheckBox matlab.ui.control.CheckBox
    ExecuteButton matlab.ui.control.Button
    ImageDisplay matlab.ui.control.Image

    % Files
    file
    drawer_list = {}
end

methods (Access = public)

    function Load(app, path)
        % Update lamp color to indicate loading
        app.statusLamp.Color = [1, 1/2, 0]; % Orange
        drawnow;

        if isfolder(path)
            path = strcat(path, '\');
        end

        totalLoadingTime = tic;

        try
            % Add file
            tic
            fprintf("\n----------------------------------\nVideo Loading\n----------------------------------\n")
            app.file = ExecutionClass(path);
            fprintf("- Video Loading took : %ds\n", round(toc))

            % Compute the mean of M0_data_video along the third dimension
            mean_M0 = mean(app.file.M0_data_video, 3);
            % Display the mean image in the uiimage component
            img = repmat(rescale(mean_M0), [1 1 3]);
            [numX, numY] = size(img);
            app.ImageDisplay.ImageSource = imresize(img, [max(numX, numY) max(numX, numY)]); % Rescale the image for display

            %% Enable buttons
            app.ExecuteButton.Enable = true;
            app.ClearButton.Enable = true;
            app.EditParametersButton.Enable = true;
            app.OverWriteCheckBox.Enable = true;
            app.EditMasksButton.Enable = true;
            app.PlayMomentsButton.Enable = true;
            app.OpenDirectoryButton.Enable = true;
            app.ReferenceDirectory.Value = path;

        catch ME
            MEdisp(ME, path)
            diary off
            % Update lamp color to indicate error
            app.statusLamp.Color = [1, 0, 0]; % Red
        end

        % Update lamp color to indicate success
        app.statusLamp.Color = [0, 1, 0]; % Green
        % Update checkbox states after loading
        app.CheckboxValueChanged();

        fprintf("----------------------------------\n")
        fprintf("- Total Load timing took : %ds\n", round(toc(totalLoadingTime)))
    end

end

% Callbacks that handle component events
methods (Access = public)

    % Code that executes after component creation
    function startupFcn(app)

        if exist("version.txt", 'file')
            v = readlines('version.txt');
            fprintf("==========================================\n " + ...
                "Welcome to EyeFlow %s\n" + ...
                "------------------------------------------\n" + ...
                "Developed by the DigitalHolographyFoundation\n" + ...
                "==========================================\n", v(1));
        end

        % Add necessary paths
        addpath("BloodFlowVelocity\", "BloodFlowVelocity\Elastography\", "BloodVolumeRate\", ...
            "BloodVolumeRate\Rheology\", "Loading\", "Parameters\", "Preprocessing\", ...
            "PulseAnalysis\", "Scripts\", "Segmentation\", "SHAnalysis\", "Tools\");

        % Set the UI title
        app.EyeFlowUIFigure.Name = ['EyeFlow ', char(v(1))];

        % Display splash screen
        displaySplashScreen();

        % Initialize checkbox states
        app.CheckboxValueChanged();
    end

    function LoadFromTxt(app)

        [selected_file, path] = uigetfile('*.txt');

        if (selected_file)
            files_lines = readlines(fullfile(path, selected_file));

            for nn = 1:length(files_lines)

                if ~isempty(files_lines(nn))
                    app.drawer_list{end + 1} = files_lines(nn);
                end

            end

        end

    end

    % Button pushed function: LoadfolderButton
    function LoadfolderButtonPushed(app, ~)

        if ~isempty(app.file)
            last_dir = app.file.directory;
        else
            last_dir = [];
        end

        % Clearing before loading
        ClearButtonPushed(app)

        % Store original WindowStyle
        originalWindowStyle = app.EyeFlowUIFigure.WindowStyle;
        app.EyeFlowUIFigure.WindowStyle = 'modal'; % Prevent minimizing

        selected_dir = uigetdir(last_dir);

        if selected_dir == 0
            fprintf(2, 'No folder selected\n');
        else
            app.Load(selected_dir);
        end

        app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore

    end

    % Button pushed function: LoadHoloButton
    function LoadHoloButtonPushed(app, ~)
        % Clearing before loading
        ClearButtonPushed(app)

        % Store original WindowStyle
        originalWindowStyle = app.EyeFlowUIFigure.WindowStyle;
        app.EyeFlowUIFigure.WindowStyle = 'modal'; % Prevent minimizing

        [selected_holo, path_holo] = uigetfile('*.holo');

        if selected_holo == 0
            fprintf(2, 'No file selected\n');
        else
            app.Load(fullfile(path_holo, selected_holo));
        end

        app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore

    end

    % Button pushed function: ExecuteButton
    function err = ExecuteButtonPushed(app, ~)

        err = [];

        if isempty(app.file)
            fprintf(2, "No input loaded.\n")
            return
        end

        % Update lamp color to indicate processing
        app.statusLamp.Color = [1, 1/2, 0]; % Orange
        drawnow;

        % Actualizes the input Parameters
        app.file.params_names = checkEyeFlowParamsFromJson(app.file.directory); % checks compatibility between found EF params and Default EF params of this version of EF.
        params = Parameters_json(app.file.directory, app.file.params_names{1});

        if params.json.NumberOfWorkers > 0 && params.json.NumberOfWorkers < app.NumberofWorkersSpinner.Limits(2)
            app.NumberofWorkersSpinner.Value = params.json.NumberOfWorkers;
        end

        parfor_arg = app.NumberofWorkersSpinner.Value;

        poolobj = gcp('nocreate'); % check if a pool already exist

        if isempty(poolobj)
            parpool(parfor_arg); % create a new pool
        elseif poolobj.NumWorkers ~= parfor_arg
            delete(poolobj); %close the current pool to create a new one with correct num of workers
            parpool(parfor_arg);
        end

        for i = 1:length(app.file.params_names)

            app.file.param_name = app.file.params_names{i};

            fprintf("==========================================\n")

            app.file.flag_segmentation = app.segmentationCheckBox.Value;
            app.file.flag_bloodFlowVelocity_analysis = app.bloodFlowAnalysisCheckBox.Value;
            app.file.flag_bloodFlowVelocity_figures = app.bloodFlowVelocityFigCheckBox.Value;
            app.file.flag_bloodVolumeRate_analysis = app.bloodVolumeRateCheckBox.Value;
            app.file.flag_bloodVolumeRate_figures = app.bloodVolumeRateFigCheckBox.Value;
            app.file.flag_spectral_analysis = app.spectralAnalysisCheckBox.Value;

            app.file.OverWrite = app.OverWriteCheckBox.Value;

            try
                app.file.ToolBoxMaster = ToolBoxClass(app.file.directory, app.file.param_name, app.file.OverWrite);

                if ~app.file.is_preprocessed
                    app.file = app.file.preprocessData();
                end

                app.file = app.file.analyzeData(app);

            catch ME
                err = ME;
                MEdisp(ME, app.file.directory)

                % Update lamp color to indicate warning
                app.statusLamp.Color = [1, 0, 0]; % Red
                diary off
            end

        end

        % Update lamp color to indicate success
        app.statusLamp.Color = [0, 1, 0]; % Green

        % Update checkbox states after execution
        app.CheckboxValueChanged();

    end

    function PlayMomentsButtonPushed(app, ~)

        try

            if app.file.is_preprocessed
                disp('inputs after preprocess.')
            else
                disp('inputs before preprocess.')
            end

            implay(rescale(app.file.M0_data_video));
            implay(rescale(app.file.M1_data_video));
            implay(rescale(app.file.M2_data_video));
        catch
            disp('Input not well loaded')
        end

    end

    function OverWriteCheckBoxChanged(app, ~)

        try
            app.file.OverWrite = app.OverWriteCheckBox.Value;
        catch
            disp('Couldnt force overwrite')
        end

    end

    % Button pushed function: ClearButton
    function ClearButtonPushed(app, ~)
        app.file = [];
        app.ReferenceDirectory.Value = "";

        app.ExecuteButton.Enable = false;
        app.ClearButton.Enable = false;
        app.EditParametersButton.Enable = false;
        app.OverWriteCheckBox.Enable = false;
        app.OpenDirectoryButton.Enable = false;
        app.EditMasksButton.Enable = false;
        app.PlayMomentsButton.Enable = false;

        % Update checkbox states
        app.CheckboxValueChanged();
    end

    % Callback function for Open Directory button
    function OpenDirectoryButtonPushed(app, ~)

        try
            % Open the directory in the system file explorer
            winopen(fullfile(app.file.directory, 'eyeflow')); % For Windows
            % Use `open(app.file.directory)` for macOS/Linux
        catch
            fprintf(2, "No valid directory loaded.\n");
        end

    end

    % Button pushed function: FolderManagementButton
    function FolderManagementButtonPushed(app, ~)
        d = dialog('Position', [300, 300, 750, 190 + length(app.drawer_list) * 14], ...
            'Color', [0.2, 0.2, 0.2], ...
            'Name', 'Folder management', ...
            'Resize', 'on', ...
            'WindowStyle', 'normal');

        txt = uicontrol('Parent', d, ...
            'Style', 'text', ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.2, 0.2, 0.2], ...
            'ForegroundColor', [0.8, 0.8, 0.8], ...
            'Position', [20, 70, 710, length(app.drawer_list) * 14], ...
            'HorizontalAlignment', 'left', ...
            'String', app.drawer_list);

        uicontrol('Parent', d, ...
            'Position', [20, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Select folder', ...
            'Callback', @select);

        uicontrol('Parent', d, ...
            'Position', [140, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Select entire folder', ...
            'Callback', @select_all);

        uicontrol('Parent', d, ...
            'Position', [260, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Clear list', ...
            'Callback', @clear_drawer);

        uicontrol('Parent', d, ...
            'Position', [380, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Load from text', ...
            'Callback', @load_from_txt);

        uicontrol('Parent', d, ...
            'Position', [500, 70, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Clear Parameters', ...
            'Callback', @clear_params);

        uicontrol('Parent', d, ...
            'Position', [500, 120, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Import Parameter', ...
            'Callback', @import_param);

        uicontrol('Parent', d, ...
            'Position', [500, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Render', ...
            'Callback', @render);

        uicontrol('Parent', d, ...
            'Position', [620, 20, 100, 25], ...
            'FontName', 'Helvetica', ...
            'BackgroundColor', [0.5, 0.5, 0.5], ...
            'ForegroundColor', [0.9 0.9 0.9], ...
            'FontWeight', 'bold', ...
            'String', 'Show Results', ...
            'Callback', @show_outputs);

        uiwait(d);

        function select_all(~, ~)
            %                 %% selection of one processed folder with uigetdir
            %                 selected_dir = uigetdir();
            %                 if (selected_dir)
            %                     app.drawer_list{end + 1} = selected_dir;
            %                 end
            %                 txt.String = app.drawer_list;
            %                 d.Position(4) = 100 + length(app.drawer_list) * 14;
            %                 txt.Position(4) = length(app.drawer_list) * 14;

            %% selection of the measurement folder with uigetdir to analyze all processed folders
            selected_dir = uigetdir();
            % List of Subfolders within the measurement folder
            tmp_dir = dir(selected_dir);
            % remove all files (isdir property is 0)
            subfoldersName = tmp_dir([tmp_dir(:).isdir]);
            % remove '.' and '..' folders
            subfoldersName = subfoldersName(~ismember({subfoldersName(:).name}, {'.', '..'}));
            subfoldersName = {subfoldersName.name};
            % remove of other folders (ex: 'config' subfolders)
            for ii = 1:length(subfoldersName)

                if contains(subfoldersName{ii}, '_HD_') || contains(subfoldersName{ii}, '_HW_')
                    app.drawer_list{end + 1} = fullfile(selected_dir, '\', subfoldersName{ii});
                    txt.String = app.drawer_list;
                    d.Position(4) = 100 + length(app.drawer_list) * 14;
                    txt.Position(4) = length(app.drawer_list) * 14;
                end

            end

        end

        function select(~, ~)
            %% selection of one processed folder with uigetdir
            selected_dir = uigetdir();

            if (selected_dir)
                app.drawer_list{end + 1} = selected_dir;
            end

            txt.String = app.drawer_list;
            d.Position(4) = 100 + length(app.drawer_list) * 14;
            txt.Position(4) = length(app.drawer_list) * 14;
        end

        function clear_drawer(~, ~)
            app.drawer_list = {};
            txt.String = app.drawer_list;
            d.Position(4) = 100 + length(app.drawer_list) * 14;
            txt.Position(4) = length(app.drawer_list) * 14;
        end

        function load_from_txt(~, ~)
            app.LoadFromTxt();
            txt.String = app.drawer_list;
            d.Position(4) = 100 + length(app.drawer_list) * 14;
            txt.Position(4) = length(app.drawer_list) * 14;
        end

        function clear_params(~, ~)
            tic
            ClearParams(app.drawer_list)
            toc
        end

        function import_param(app, ~)
            tic
            % Store the current WindowStyle of the main GUI
            originalWindowStyle = app.EyeFlowUIFigure.WindowStyle;

            % Temporarily set the WindowStyle to 'modal' to prevent minimizing
            app.EyeFlowUIFigure.WindowStyle = 'modal';

            % Open the file selection dialog
            [selected_json, json_path] = uigetfile('*.json');

            if selected_json == 0
                disp('No file selected');
                % Restore the original WindowStyle
                app.EyeFlowUIFigure.WindowStyle = originalWindowStyle;
                return;
            end

            % Restore the original WindowStyle
            app.EyeFlowUIFigure.WindowStyle = originalWindowStyle;

            % Process the selected file
            for ind = 1:length(app.drawer_list)
                path_json = fullfile(app.drawer_list{ind}, 'eyeflow', 'json');

                if ~isfolder(path_json)
                    mkdir(path_json);
                end

                copyfile(fullfile(json_path, selected_json), path_json);

                % Get idx for renaming
                idx = 0;
                list_dir = dir(path_json);

                for i = 1:numel(list_dir)
                    match = regexp(list_dir(i).name, '\d+$', 'match');

                    if ~isempty(match) && str2double(match{1}) >= idx
                        idx = str2double(match{1}); % Suffix
                    end

                end

                % Renaming
                copyfile(fullfile(path_json, selected_json), fullfile(path_json, sprintf('InputEyeFlowParams_%d.json', idx)));
                delete(fullfile(path_json, selected_json));
            end

            toc
        end

        function render(~, ~)

            num_drawers = length(app.drawer_list);
            error_list = cell(1, num_drawers);
            faulty_folders = cell(1, num_drawers);
            error_count = 0;

            for i = 1:num_drawers
                tic

                app.Load(app.drawer_list{i});
                ME = app.ExecuteButtonPushed();

                if ~isempty(ME)
                    error_count = error_count + 1;
                    error_list{error_count} = ME;
                    faulty_folders{error_count} = app.drawer_list{i};
                end

                app.ClearButtonPushed();
                toc
            end

            if error_count > 0
                disp('Errors in rendering:')

                for i = 1:error_count
                    exception = error_list{i};
                    path = faulty_folders{i};
                    fprintf(2, "==========================================\nERROR N°:%d\n==========================================\n", i)
                    fprintf(2, 'Folder : %s\n', path)
                    fprintf(2, "%s\n", exception.identifier)
                    fprintf(2, "%s\n", exception.message)

                    for stackIdx = 1:size(exception.stack, 1)
                        fprintf(2, "%s : %s, line : %d\n", exception.stack(stackIdx).file, exception.stack(stackIdx).name, exception.stack(stackIdx).line);
                    end

                end

            end

        end

        function show_outputs(~, ~)
            out_dir_path = fullfile(app.drawer_list{1}, 'Multiple_Results');
            mkdir(out_dir_path) % creates if it doesn't exists
            tic
            ShowOutputs(app.drawer_list, out_dir_path)
            toc
        end

        delete(d);
    end

    % Button pushed function: MaskToolButtonPushed
    function MaskToolButtonPushed(app, ~)

        parfor_arg = app.NumberofWorkersSpinner.Value;
        poolobj = gcp('nocreate'); % check if a pool already exist

        if isempty(poolobj)
            parpool(parfor_arg); % create a new pool
        elseif poolobj.NumWorkers ~= parfor_arg
            delete(poolobj); %close the current pool to create a new one with correct num of workers
            parpool(parfor_arg);
        end

        PreviewMasks(app);
    end

    % Button pushed function: EditParametersButton
    function EditParametersButtonPushed(app, ~)

        try
            main_path = fullfile(app.file.directory, 'eyeflow');

            if isfile(fullfile(main_path, 'json', app.file.param_name))
                disp(['opening : ', fullfile(main_path, 'json', app.file.param_name)])
                winopen(fullfile(main_path, 'json', app.file.param_name));
            else
                disp(['couldn''t open : ', fullfile(main_path, 'json', app.file.param_name)])
            end

        catch
            fprintf(2, "no input loaded\n")
        end

    end

    % Button pushed function: EditMasksButton
    function EditMasksButtonPushed(app, ~)
        ToolBox = getGlobalToolBox;

        if isempty(ToolBox) || ~strcmp(app.file.directory, ToolBox.EF_path)
            ToolBox = ToolBoxClass(app.file.directory, app.file.param_name, 1);
        end

        if isempty(app.file)

            if ~isfolder(fullfile(ToolBox.path_main, 'mask'))
                mkdir(fullfile(ToolBox.path_main, 'mask'))
            end

            try
                winopen(fullfile(ToolBox.path_main, 'mask'));
            catch
                disp("opening failed.")
            end

            try
                list_dir = dir(ToolBox.path_main);
                idx = 0;

                for i = 1:length(list_dir)

                    if contains(list_dir(i).name, ToolBox.EF_name)
                        match = regexp(list_dir(i).name, '\d+$', 'match');

                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end

                    end

                end

                path_dir = fullfile(ToolBox.path_main, ToolBox.folder_name);

                disp(['Copying from : ', fullfile(path_dir, 'png', 'mask')])
                copyfile(fullfile(path_dir, 'png', 'mask', sprintf("%s_maskArtery.png", ToolBox.main_foldername)), fullfile(ToolBox.path_main, 'mask', 'MaskArtery.png'));
                copyfile(fullfile(path_dir, 'png', 'mask', sprintf("%s_maskVein.png", ToolBox.main_foldername)), fullfile(ToolBox.path_main, 'mask', 'MaskVein.png'));
            catch
                disp("last auto mask copying failed.")
            end

            try

                copyfile(fullfile(ToolBox.EF_path, 'png', sprintf("%s_M0.png", ToolBox.main_foldername)), fullfile(ToolBox.path_main, 'mask', 'M0.png'));
                folder_name = strcat(ToolBox.main_foldername, '_EF');
                list_dir = dir(ToolBox.path_main);
                idx = 0;

                for i = 1:length(list_dir)

                    if contains(list_dir(i).name, folder_name)
                        match = regexp(list_dir(i).name, '\d+$', 'match');

                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end

                    end

                end

                folder_name = sprintf('%s_%d', folder_name, idx);
                copyfile(fullfile(path_dir, 'gif', sprintf("%s_M0.gif", folder_name)), fullfile(ToolBox.path_main, 'mask', 'M0.gif'));
            catch

                disp("last M0 png and gif copying failed")
            end

            try
                v = VideoReader(fullfile(ToolBox.EF_path, 'avi', sprintf("%s_M0.avi", ToolBox.main_foldername)));
                M0_video = read(v); clear v;
                M0_video = rescale(single(squeeze(mean(M0_video, 3))));
                sz = size(M0_video);
                [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_video, diskMask(sz(1), sz(2), 0.45));
                diasysArtery = M0_Systole_img - M0_Diastole_img;
                RGBdiasys = labDuoImage(mean(M0_video, 3), diasysArtery);
                imwrite(RGBdiasys, fullfile(ToolBox.path_main, 'mask', 'DiaSysRGB.png'), 'png');
            catch

                disp("Diasys png failed")

            end

            % try
            % %   Commented until further fixes MESSAGE TO ZACHARIE
            %     openmaskinpaintnet(fullfile(ToolBox.path_main,'mask','M0.png'), fullfile(ToolBox.path_main,'mask','DiaSysRGB.png'));
            % catch
            %     disp("paint.net macro failed")
            % end

        else

            fprintf(2, "no input loaded\n")

        end

    end

    function CheckboxValueChanged(app, ~)
        % Callback function triggered when any checkbox (except OverWriteCheckBox) is clicked.
        % This function enforces the rules for enabling/disabling checkboxes.

        % Segmentation Checkbox is always enabled
        app.segmentationCheckBox.Enable = true;

        % Enable/disable bloodFlowAnalysisCheckBox and spectralAnalysisCheckBox
        is_segmented = false;
        is_pulseAnalyzed = false;
        is_bloodVolumeRateAnalyzed = false;

        if ~isempty(app.file)
            is_segmented = app.file.is_segmented;
            is_pulseAnalyzed = app.file.is_pulseAnalyzed;
            is_bloodVolumeRateAnalyzed = app.file.is_bloodVolumeRateAnalyzed;
        end

        if app.segmentationCheckBox.Value || is_segmented
            app.bloodFlowAnalysisCheckBox.Enable = true;
            app.spectralAnalysisCheckBox.Enable = true;
        else
            app.bloodFlowAnalysisCheckBox.Enable = false;
            app.spectralAnalysisCheckBox.Enable = false;
            app.bloodFlowAnalysisCheckBox.Value = false; % Turn off if disabled
            app.spectralAnalysisCheckBox.Value = false; % Turn off if disabled
        end

        % Enable/disable bloodFlowVelocityFigCheckBox and bloodVolumeRateCheckBox
        if app.bloodFlowAnalysisCheckBox.Value || is_pulseAnalyzed
            app.bloodFlowVelocityFigCheckBox.Enable = true;
            app.bloodVolumeRateCheckBox.Enable = true;
        else
            app.bloodFlowVelocityFigCheckBox.Enable = false;
            app.bloodVolumeRateCheckBox.Enable = false;
            app.bloodFlowVelocityFigCheckBox.Value = false; % Turn off if disabled
            app.bloodVolumeRateCheckBox.Value = false; % Turn off if disabled
        end

        % Enable/disable bloodVolumeRateFigCheckBox
        if app.bloodVolumeRateCheckBox.Value || is_bloodVolumeRateAnalyzed
            app.bloodVolumeRateFigCheckBox.Enable = true;
        else
            app.bloodVolumeRateFigCheckBox.Enable = false;
            app.bloodVolumeRateFigCheckBox.Value = false; % Turn off if disabled
        end

    end

end

% Component initialization
methods (Access = private)

    % Create UIFigure and components
    function createComponents(app)
        % Create UIFigure and components with the specified layout.

        pathToMLAPP = fileparts(mfilename('fullpath'));

        % Create EyeFlowUIFigure and hide until all components are created
        app.EyeFlowUIFigure = uifigure('Visible', 'off');
        app.EyeFlowUIFigure.Color = [0.149 0.149 0.149];
        app.EyeFlowUIFigure.Position = [100 100 1050 421];
        app.EyeFlowUIFigure.Name = 'EyeFlow';
        app.EyeFlowUIFigure.Icon = fullfile(pathToMLAPP, 'eyeflow_logo.png');

        % Create a grid layout to manage resizing
        grid = uigridlayout(app.EyeFlowUIFigure);
        grid.RowHeight = {'fit', '1x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
        grid.ColumnWidth = {'1x', '1x', '1x', '1x', '3x'};
        grid.BackgroundColor = [0.149, 0.149, 0.149];

        % Top Row: Load Folder, Load Holo, Clear, Folder Management
        app.LoadFolderButton = uibutton(grid, 'push');
        app.LoadFolderButton.ButtonPushedFcn = createCallbackFcn(app, @LoadfolderButtonPushed, true);
        app.LoadFolderButton.BackgroundColor = [0.502 0.502 0.502];
        app.LoadFolderButton.FontSize = 16;
        app.LoadFolderButton.FontColor = [0.9412 0.9412 0.9412];
        app.LoadFolderButton.Layout.Row = 1;
        app.LoadFolderButton.Layout.Column = 1;
        app.LoadFolderButton.Text = 'Load Folder';

        app.LoadHoloButton = uibutton(grid, 'push');
        app.LoadHoloButton.ButtonPushedFcn = createCallbackFcn(app, @LoadHoloButtonPushed, true);
        app.LoadHoloButton.BackgroundColor = [0.502 0.502 0.502];
        app.LoadHoloButton.FontSize = 16;
        app.LoadHoloButton.FontColor = [0.9412 0.9412 0.9412];
        app.LoadHoloButton.Layout.Row = 1;
        app.LoadHoloButton.Layout.Column = 2;
        app.LoadHoloButton.Text = 'Load Holo';

        app.ClearButton = uibutton(grid, 'push');
        app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
        app.ClearButton.BackgroundColor = [0.502 0.502 0.502];
        app.ClearButton.FontSize = 16;
        app.ClearButton.FontColor = [0.9412 0.9412 0.9412];
        app.ClearButton.Enable = 'off';
        app.ClearButton.Layout.Row = 1;
        app.ClearButton.Layout.Column = 3;
        app.ClearButton.Text = 'Clear';

        app.FolderManagementButton = uibutton(grid, 'push');
        app.FolderManagementButton.ButtonPushedFcn = createCallbackFcn(app, @FolderManagementButtonPushed, true);
        app.FolderManagementButton.BackgroundColor = [0.502 0.502 0.502];
        app.FolderManagementButton.FontSize = 16;
        app.FolderManagementButton.FontColor = [0.9412 0.9412 0.9412];
        app.FolderManagementButton.Layout.Row = 1;
        app.FolderManagementButton.Layout.Column = 4;
        app.FolderManagementButton.Text = 'Folder Management';

        % Second Row: Directory

        dirgrid = uigridlayout(grid);
        dirgrid.Layout.Row = 2;
        dirgrid.Layout.Column = [1, 4];
        dirgrid.RowHeight = {'1x'};
        dirgrid.ColumnWidth = {'1x', 20};
        dirgrid.BackgroundColor = [0.149, 0.149, 0.149];

        app.ReferenceDirectory = uitextarea(dirgrid);
        app.ReferenceDirectory.Editable = 'off';
        app.ReferenceDirectory.FontSize = 16;
        app.ReferenceDirectory.FontColor = [0.9412 0.9412 0.9412];
        app.ReferenceDirectory.BackgroundColor = [0.149 0.149 0.149];
        app.ReferenceDirectory.Layout.Row = 1;
        app.ReferenceDirectory.Layout.Column = 1;

        % Create Lamp
        app.statusLamp = uilamp(dirgrid);
        app.statusLamp.Layout.Row = 1;
        app.statusLamp.Layout.Column = 2;
        app.statusLamp.Color = [0 1 0];

        % Third Row: Edit Parameters, Edit Masks, Play Inputs, Preview Masks
        app.EditParametersButton = uibutton(grid, 'push');
        app.EditParametersButton.ButtonPushedFcn = createCallbackFcn(app, @EditParametersButtonPushed, true);
        app.EditParametersButton.BackgroundColor = [0.502 0.502 0.502];
        app.EditParametersButton.FontSize = 16;
        app.EditParametersButton.FontColor = [0.9412 0.9412 0.9412];
        app.EditParametersButton.Layout.Row = 3;
        app.EditParametersButton.Layout.Column = 1;
        app.EditParametersButton.Enable = 'off';
        app.EditParametersButton.Text = 'Edit Parameters';
        app.EditParametersButton.Tooltip = 'Find the eyeflow parameters here.';

        app.EditMasksButton = uibutton(grid, 'push');
        app.EditMasksButton.ButtonPushedFcn = createCallbackFcn(app, @EditMasksButtonPushed, true);
        app.EditMasksButton.BackgroundColor = [0.502 0.502 0.502];
        app.EditMasksButton.FontSize = 16;
        app.EditMasksButton.FontColor = [0.9412 0.9412 0.9412];
        app.EditMasksButton.Layout.Row = 3;
        app.EditMasksButton.Layout.Column = 2;
        app.EditMasksButton.Enable = 'off';
        app.EditMasksButton.Text = 'Edit Masks';
        app.EditMasksButton.Tooltip = 'Open mask folder and use forceMaskArtery.png and forceMaskVein.png to force the segmentation';

        app.PlayMomentsButton = uibutton(grid, 'push');
        app.PlayMomentsButton.ButtonPushedFcn = createCallbackFcn(app, @PlayMomentsButtonPushed, true);
        app.PlayMomentsButton.BackgroundColor = [0.502 0.502 0.502];
        app.PlayMomentsButton.FontSize = 16;
        app.PlayMomentsButton.FontColor = [0.9412 0.9412 0.9412];
        app.PlayMomentsButton.Layout.Row = 3;
        app.PlayMomentsButton.Layout.Column = 3;
        app.PlayMomentsButton.Enable = 'off';
        app.PlayMomentsButton.Text = 'Play Moments';

        app.MaskToolButton = uibutton(grid, 'push');
        app.MaskToolButton.ButtonPushedFcn = createCallbackFcn(app, @MaskToolButtonPushed, true);
        app.MaskToolButton.BackgroundColor = [0.502 0.502 0.502];
        app.MaskToolButton.FontSize = 16;
        app.MaskToolButton.FontColor = [0.9412 0.9412 0.9412];
        app.MaskToolButton.Layout.Row = 3;
        app.MaskToolButton.Layout.Column = 4;
        app.MaskToolButton.Text = 'Mask Tool';

        % Checkboxes: Segmentation, Pulse analysis, Blood Flow Velocity, Blood Volume Rate, SH analysis
        app.segmentationCheckBox = uicheckbox(grid);
        app.segmentationCheckBox.Text = 'Segmentation';
        app.segmentationCheckBox.FontSize = 16;
        app.segmentationCheckBox.FontColor = [1 1 1];
        app.segmentationCheckBox.Layout.Row = 4;
        app.segmentationCheckBox.Layout.Column = [1, 4];
        app.segmentationCheckBox.Value = true;
        app.segmentationCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        app.bloodFlowAnalysisCheckBox = uicheckbox(grid);
        app.bloodFlowAnalysisCheckBox.Text = 'Blood Flow Analysis';
        app.bloodFlowAnalysisCheckBox.FontSize = 16;
        app.bloodFlowAnalysisCheckBox.FontColor = [1 1 1];
        app.bloodFlowAnalysisCheckBox.Layout.Row = 5;
        app.bloodFlowAnalysisCheckBox.Layout.Column = [1, 2];
        app.bloodFlowAnalysisCheckBox.Value = true;
        app.bloodFlowAnalysisCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        app.bloodFlowVelocityFigCheckBox = uicheckbox(grid);
        app.bloodFlowVelocityFigCheckBox.Text = 'Blood Flow Velocity Figures';
        app.bloodFlowVelocityFigCheckBox.FontSize = 16;
        app.bloodFlowVelocityFigCheckBox.FontColor = [1 1 1];
        app.bloodFlowVelocityFigCheckBox.Layout.Row = 5;
        app.bloodFlowVelocityFigCheckBox.Layout.Column = [3, 4];
        app.bloodFlowVelocityFigCheckBox.Value = true;
        app.bloodFlowVelocityFigCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        app.bloodVolumeRateCheckBox = uicheckbox(grid);
        app.bloodVolumeRateCheckBox.Text = 'Blood Volume Rate Analysis';
        app.bloodVolumeRateCheckBox.FontSize = 16;
        app.bloodVolumeRateCheckBox.FontColor = [1 1 1];
        app.bloodVolumeRateCheckBox.Layout.Row = 6;
        app.bloodVolumeRateCheckBox.Layout.Column = [1, 2];
        app.bloodVolumeRateCheckBox.Value = true;
        app.bloodVolumeRateCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        app.bloodVolumeRateFigCheckBox = uicheckbox(grid);
        app.bloodVolumeRateFigCheckBox.Text = 'Blood Volume Rate Figures';
        app.bloodVolumeRateFigCheckBox.FontSize = 16;
        app.bloodVolumeRateFigCheckBox.FontColor = [1 1 1];
        app.bloodVolumeRateFigCheckBox.Layout.Row = 6;
        app.bloodVolumeRateFigCheckBox.Layout.Column = [3, 4];
        app.bloodVolumeRateFigCheckBox.Value = true;
        app.bloodVolumeRateFigCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        app.spectralAnalysisCheckBox = uicheckbox(grid);
        app.spectralAnalysisCheckBox.Text = 'Spectral analysis';
        app.spectralAnalysisCheckBox.FontSize = 16;
        app.spectralAnalysisCheckBox.FontColor = [1 1 1];
        app.spectralAnalysisCheckBox.Layout.Row = 7;
        app.spectralAnalysisCheckBox.Layout.Column = [1, 4];
        app.spectralAnalysisCheckBox.Enable = true;
        app.spectralAnalysisCheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckboxValueChanged, true);

        % Bottom Left: Execute Button
        app.ExecuteButton = uibutton(grid, 'push');
        app.ExecuteButton.ButtonPushedFcn = createCallbackFcn(app, @ExecuteButtonPushed, true);
        app.ExecuteButton.BackgroundColor = [0.502 0.502 0.502];
        app.ExecuteButton.FontSize = 16;
        app.ExecuteButton.FontColor = [0.9412 0.9412 0.9412];
        app.ExecuteButton.Enable = 'off';
        app.ExecuteButton.Layout.Row = 8;
        app.ExecuteButton.Layout.Column = 1;
        app.ExecuteButton.Text = 'Execute';

        % Number of Workers Spinner and Label
        app.NumberofWorkersSpinnerLabel = uilabel(grid);
        app.NumberofWorkersSpinnerLabel.HorizontalAlignment = 'right';
        app.NumberofWorkersSpinnerLabel.FontColor = [0.902 0.902 0.902];
        app.NumberofWorkersSpinnerLabel.FontSize = 16;
        app.NumberofWorkersSpinnerLabel.Layout.Row = 8;
        app.NumberofWorkersSpinnerLabel.Layout.Column = 3;
        app.NumberofWorkersSpinnerLabel.Text = 'Number of Workers';

        app.NumberofWorkersSpinner = uispinner(grid);
        app.NumberofWorkersSpinner.FontSize = 16;
        app.NumberofWorkersSpinner.Layout.Row = 8;
        app.NumberofWorkersSpinner.Layout.Column = 4;
        maxWorkers = parcluster('local').NumWorkers;
        app.NumberofWorkersSpinner.Limits = [0 maxWorkers]; % Ensure valid range
        app.NumberofWorkersSpinner.Value = min(10, floor(maxWorkers / 2)); % Default to 10 or max available

        % Bottom Right: Overwrite Checkbox
        app.OverWriteCheckBox = uicheckbox(grid);
        app.OverWriteCheckBox.Text = 'Overwrite';
        app.OverWriteCheckBox.FontSize = 16;
        app.OverWriteCheckBox.FontColor = [0.8 0.8 0.8];
        app.OverWriteCheckBox.Layout.Row = 8;
        app.OverWriteCheckBox.Layout.Column = 2;
        app.OverWriteCheckBox.Value = false;
        app.OverWriteCheckBox.Enable = 'off';
        app.OverWriteCheckBox.ValueChangedFcn = createCallbackFcn(app, @OverWriteCheckBoxChanged, true);
        app.OverWriteCheckBox.Tooltip = 'Overwrite the new results in the last EF_ result folder (to save space\n Warning: it will supress previous figures.';

        % Add a new column for the image
        app.ImageDisplay = uiimage(grid);
        app.ImageDisplay.Layout.Row = [1, 8]; % Span all rows
        app.ImageDisplay.Layout.Column = 5; % Place in the new column
        app.ImageDisplay.ScaleMethod = 'fit'; % Adjust the image to fit the component

        % Add the new button under Preview Masks
        app.OpenDirectoryButton = uibutton(grid, 'push');
        app.OpenDirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @OpenDirectoryButtonPushed, true);
        app.OpenDirectoryButton.BackgroundColor = [0.502 0.502 0.502];
        app.OpenDirectoryButton.FontSize = 16;
        app.OpenDirectoryButton.FontColor = [0.9412 0.9412 0.9412];
        app.OpenDirectoryButton.Layout.Row = 4; % Adjust the row as needed
        app.OpenDirectoryButton.Layout.Column = 4; % Same column as Preview Masks
        app.OpenDirectoryButton.Text = 'Open Directory';
        app.OpenDirectoryButton.Enable = 'off'; % Disabled by default

        % Show the figure after all components are created
        app.EyeFlowUIFigure.Visible = 'on';
    end

end

% App creation and deletion
methods (Access = public)

    % Construct app
    function app = eyeflow

        % Create UIFigure and components
        createComponents(app)

        % Register the app with App Designer
        registerApp(app, app.EyeFlowUIFigure)

        % Execute the startup function
        runStartupFcn(app, @startupFcn)

        if nargout == 0
            clear app
        end

    end

    % Code that executes before app deletion
    function delete(app)

        % Delete UIFigure when app is deleted
        delete(app.EyeFlowUIFigure)
    end

end

end
