classdef eyeflow < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        EyeFlowUIFigure             matlab.ui.Figure
        PlayInputsButton              matlab.ui.control.Button
        EditParametersButton          matlab.ui.control.Button
        EditMasksButton               matlab.ui.control.Button
        NumberofWorkersSpinner        matlab.ui.control.Spinner
        NumberofWorkersSpinnerLabel   matlab.ui.control.Label
        SHanalysisCheckBox            matlab.ui.control.CheckBox
        FolderManagementButton        matlab.ui.control.Button
        PreviewMasksButton            matlab.ui.control.Button
        SegmentationCheckBox          matlab.ui.control.CheckBox
        PulseanalysisCheckBox         matlab.ui.control.CheckBox
        velocityCheckBox              matlab.ui.control.CheckBox
        bloodVolumeRateCheckBox       matlab.ui.control.CheckBox
        ReferenceDirectory            matlab.ui.control.TextArea
        OverWriteCheckBox             matlab.ui.control.CheckBox
        Lamp                          matlab.ui.control.Lamp
        ClearButton                   matlab.ui.control.Button
        LoadfolderButton              matlab.ui.control.Button
        LoadHoloButton                matlab.ui.control.Button
        ExecuteButton                 matlab.ui.control.Button
        file
        drawer_list = {}
        flag_is_load
    end

    methods (Access = private)
        function Load(app, path)
            app.Lamp.Color = [1, 0, 0];
            drawnow;
            if isfolder(path)
                path = strcat(path, '\');
            end

            totalLoadingTime = tic;

            try
                % add file
                tic
                fprintf("\n----------------------------------\nVideo Loading\n----------------------------------\n")
                app.file = ExecutionClass(path);
                fprintf("- Video Loading took : %ds\n", round(toc))

                % Compute the mean of M0_data_video along the third dimension
                mean_M0 = mean(app.file.M0_data_video, 3);
                % Display the mean image in the uiimage component
                app.ImageDisplay.ImageSource = rescale(mean_M0); % Rescale the image for display

                %% End
                app.LoadfolderButton.Enable = true ;
                app.ExecuteButton.Enable = true ;
                app.ClearButton.Enable = true ;
                app.EditParametersButton.Enable = true;
                % FIXME app.ReferenceDirectory.Value = fullfile(path,file)
                app.ReferenceDirectory.Value = path ;
                app.Lamp.Color = [0, 1, 0];
                app.flag_is_load = true;

            catch ME

                fprintf(2, "==========================================\nERROR\n==========================================\n")
                fprintf(2, 'Error while loading : %s\n', path)
                fprintf(2, "%s\n",ME.identifier)
                fprintf(2, "%s\n",ME.message)
                % for i = 1:size(exception.stack,1)
                %     stack = sprintf('%s : %s, line : %d \n', exception.stack(i).file, exception.stack(i).name, exception.stack(i).line);
                %     fprintf(stack);
                % end

                if ME.identifier == "MATLAB:audiovideo:VideoReader:FileNotFound"

                    fprintf(2, "No Raw File was found, please check 'save raw files' in HoloDoppler\n")

                else

                    for i = 1:numel(ME.stack)
                        fprintf(2, "%s", ME.stack(i))
                    end

                end

                fprintf(2, "==========================================\n")

                diary off
                app.Lamp.Color = [1, 1/2, 0];

            end

            fprintf("----------------------------------\n")
            fprintf("- Total Load timing took : %ds\n", round(toc(totalLoadingTime)))

        end
    end
    methods (Access = private)

        function LoadFromTxt(app)

            [selected_file,path] = uigetfile('*.txt');
            if (selected_file)
                files_lines = readlines(fullfile(path,selected_file));
                for nn = 1:length(files_lines)
                    if ~isempty(files_lines(nn))
                        app.drawer_list{end + 1} = files_lines(nn);
                    end
                end
            end

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


            % Initialize checkboxes and flags
            app.flag_is_load = false;

            % Display splash screen
            displaySplashScreen();
        end

        % Button pushed function: LoadfolderButton
        function LoadfolderButtonPushed(app, ~)
            % Clearing before loading
            if ~isempty(app.file)
                last_dir = app.file.directory;
            else
                last_dir = [];
            end

            app.file = [];
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;

            clear Parameters_json;

            if (app.flag_is_load)
                disp("Files already loaded");
                app.ClearButton.Enable = true;
                app.EditParametersButton.Enable = true;
            else
                % Store original WindowStyle
                originalWindowStyle = app.EyeFlowUIFigure.WindowStyle;
                app.EyeFlowUIFigure.WindowStyle = 'modal'; % Prevent minimizing

                selected_dir = uigetdir(last_dir);
                if selected_dir == 0
                    fprintf('No folder selected');
                    app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore
                    return;
                end

                app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore
                app.flag_is_load = true;
                app.Load(selected_dir);
            end
        end

        function LoadHoloButtonPushed(app, ~)
            % Clearing before loading
            app.file = [];
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;

            clear Parameters_json;

            if (app.flag_is_load)
                disp("Files already loaded");
                app.ClearButton.Enable = true;
                app.EditParametersButton.Enable = true;
            else
                % Store original WindowStyle
                originalWindowStyle = app.EyeFlowUIFigure.WindowStyle;
                app.EyeFlowUIFigure.WindowStyle = 'modal'; % Prevent minimizing

                [selected_holo, path_holo] = uigetfile('*.holo');
                if selected_holo == 0
                    disp('No file selected');
                    app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore
                    return;
                end

                app.EyeFlowUIFigure.WindowStyle = originalWindowStyle; % Restore
                app.flag_is_load = true;
                app.Load(fullfile(path_holo, selected_holo));
            end

            app.Lamp.Color = [0, 1, 0];
        end

        % Button pushed function: ExecuteButton
        function ExecuteButtonPushed(app, ~)
            if ~app.flag_is_load
                fprintf(2, "no input loaded\n")
                return
            end

            warning('off');
            parfor_arg = app.NumberofWorkersSpinner.Value ;

            poolobj = gcp('nocreate'); % check if a pool already exist
            if isempty(poolobj)
                parpool(parfor_arg); % create a new pool
            elseif poolobj.NumWorkers ~= parfor_arg
                delete(poolobj); %close the current pool to create a new one with correct num of workers
                parpool(parfor_arg);
            end

            clear Parameters_json
            app.Lamp.Color = [1, 0, 0];
            drawnow;

            % Actualizes the input Parameters
            app.file.params_names = checkEyeFlowParamsFromJson(app.file.directory); % checks compatibility between found EF params and Default EF params of this version of EF.


            for i = 1:length(app.file.params_names)

                app.file.param_name = app.file.params_names{i};

                fprintf("==========================================\n")
                app.file.flag_Segmentation = app.SegmentationCheckBox.Value;
                app.file.flag_SH_analysis = app.SHanalysisCheckBox.Value;
                app.file.flag_Pulse_analysis = app.PulseanalysisCheckBox.Value;
                app.file.flag_velocity_analysis = app.velocityCheckBox.Value;
                app.file.flag_bloodVolumeRate_analysis = app.bloodVolumeRateCheckBox.Value;

                app.file.OverWrite = app.OverWriteCheckBox.Value;

                try
                    if ~app.file.is_preprocessed
                        app.file = app.file.preprocessData();
                    end
                    app.file = app.file.analyzeData();

                catch ME

                    diary off

                    fprintf(2,"==========================================\nERROR\n==========================================\n");

                    fprintf(2, 'Error with file : %s\n%s\n%s', app.file.directory, ME.identifier, ME.message);

                    for stackIdx = 1:size(ME.stack, 1)
                        fprintf(2,"%s : %s, line : %d\n", ME.stack(stackIdx).file, ME.stack(stackIdx).name, ME.stack(stackIdx).line);
                    end

                    fprintf(2,"==========================================\n");

                    % Update lamp color to indicate warning
                    app.Lamp.Color = [1, 0.5, 0]; % Orange
                end
            end
            app.Lamp.Color = [0, 1, 0];
        end

        function PlayInputsButtonPushed(app, ~)
            if ~app.flag_is_load
                fprintf(2, "no input loaded\n")
                return
            end
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
            if ~app.flag_is_load
                fprintf(2, "no input loaded\n")
                return
            end
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
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;

            clear Parameters_json

        end

        % Button pushed function: FolderManagementButton
        function FolderManagementButtonPushed(app, ~)
            d = dialog('Position', [300, 300, 750, 190 + length(app.drawer_list) * 14],...
                'Color', [0.2, 0.2, 0.2],...
                'Name', 'Folder management',...
                'Resize', 'on',...
                'WindowStyle', 'normal');

            txt = uicontrol('Parent', d,...
                'Style', 'text',...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.2, 0.2, 0.2],...
                'ForegroundColor', [0.8, 0.8, 0.8],...
                'Position', [20, 70, 710, length(app.drawer_list) * 14],...
                'HorizontalAlignment', 'left',...
                'String', app.drawer_list);

            uicontrol('Parent', d,...
                'Position', [20, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Select folder',...
                'Callback', @select);

            uicontrol('Parent', d,...
                'Position', [140, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Select entire folder',...
                'Callback', @select_all);

            uicontrol('Parent', d,...
                'Position', [260, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Clear list',...
                'Callback', @clear_drawer);

            uicontrol('Parent', d,...
                'Position', [380, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Load from text',...
                'Callback', @load_from_txt);

            uicontrol('Parent', d,...
                'Position', [500, 70, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Clear Parameters',...
                'Callback', @clear_params);

            uicontrol('Parent', d,...
                'Position', [500, 120, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Import Parameter',...
                'Callback', @import_param);


            uicontrol('Parent', d,...
                'Position', [500, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Render',...
                'Callback', @render);

            uicontrol('Parent', d,...
                'Position', [620, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Show Results',...
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
                subfoldersName = subfoldersName(~ismember({subfoldersName(:).name},{'.','..'}));
                subfoldersName = {subfoldersName.name};
                % remove of other folders (ex: 'config' subfolders)
                for ii=1:length(subfoldersName)
                    if contains(subfoldersName{ii}, '_HD_') || contains(subfoldersName{ii}, '_HW_')
                        app.drawer_list{end + 1} = fullfile(selected_dir,'\',subfoldersName{ii});
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
                for i = 1:length(app.drawer_list)
                    tic
                    try
                    app.Load(app.drawer_list{i});
                    app.ExecuteButtonPushed();
                    catch
                        fprintf(2, 'Error in rendering: %s', ME.message)
                    end
                    app.ClearButtonPushed();
                    toc
                end
                %                 clear app.drawer_list
            end

            function show_outputs(~, ~)
                out_dir_path = fullfile(app.drawer_list{1},'Multiple_Results');
                mkdir(out_dir_path) % creates if it doesn't exists
                tic
                ShowOutputs(app.drawer_list,out_dir_path)
                toc
            end
            delete(d);
        end

        % Button pushed function: PreviewMasksButtonPushed
        function PreviewMasksButtonPushed(app, ~)

            parfor_arg = app.NumberofWorkersSpinner.Value ;
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
            if (app.flag_is_load)
                main_path = fullfile(app.file.directory, 'eyeflow');
                if isfile(fullfile(main_path,'json',app.file.param_name))
                    disp(['opening : ', fullfile(main_path,'json',app.file.param_name)])
                    winopen(fullfile(main_path,'json',app.file.param_name));
                else
                    disp(['couldn''t open : ',fullfile(main_path,'json',app.file.param_name)])
                end
            else
                fprintf(2, "no input loaded\n")
            end
        end

        % Button pushed function: EditMasksButton
        function EditMasksButtonPushed(app, ~)
            TB = getGlobalToolBox;

            if isempty(TB)
                TB = ToolBoxClass(app.file.directory, app.file.param_name, 1);
            end

            if (app.flag_is_load)
                if ~isfolder(fullfile(TB.path_main,'mask'))
                    mkdir(fullfile(TB.path_main,'mask'))
                end
                try
                    winopen(fullfile(TB.path_main,'mask'));
                catch
                    disp("opening failed.")
                end
                try
                    list_dir = dir(TB.path_main);
                    idx=0;
                    for i=1:length(list_dir)
                        if contains(list_dir(i).name, TB.EF_name)
                            match = regexp(list_dir(i).name, '\d+$', 'match');
                            if ~isempty(match) && str2double(match{1}) >= idx
                                idx = str2double(match{1}) ; %suffix
                            end
                        end
                    end
                    path_dir = fullfile(TB.path_main, TB.folder_name);

                    disp(['Copying from : ',fullfile(path_dir,'png','mask')])
                    copyfile(fullfile(path_dir,'png', 'mask', sprintf("%s_maskArtery.png", TB.main_foldername)), fullfile(TB.path_main, 'mask', 'MaskArtery.png'));
                    copyfile(fullfile(path_dir,'png', 'mask', sprintf("%s_maskVein.png", TB.main_foldername)), fullfile(TB.path_main, 'mask', 'MaskVein.png'));
                catch
                    disp("last auto mask copying failed.")
                end
                try

                    copyfile(fullfile(TB.path,'png',sprintf("%s_M0.png",TB.main_foldername)),fullfile(TB.path_main,'mask','M0.png'));
                    folder_name = strcat(TB.main_foldername, '_EF');
                    list_dir = dir(TB.path_main);
                    idx=0;
                    for i=1:length(list_dir)
                        if contains(list_dir(i).name, folder_name)
                            match = regexp(list_dir(i).name, '\d+$', 'match');
                            if ~isempty(match) && str2double(match{1}) >= idx
                                idx = str2double(match{1}) ; %suffix
                            end
                        end
                    end
                    folder_name = sprintf('%s_%d', folder_name, idx);
                    copyfile(fullfile(path_dir,'gif',sprintf("%s_M0.gif",folder_name)),fullfile(TB.path_main,'mask','M0.gif'));
                catch

                    disp("last M0 png and gif copying failed")
                end

                try
                    v = VideoReader(fullfile(TB.path,'avi',sprintf("%s_M0.avi",TB.main_foldername)));
                    M0_video = read(v); clear v;
                    M0_video = rescale(single(squeeze(mean(M0_video,3))));
                    sz = size(M0_video);
                    [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_video,  diskMask(sz(1),sz(2),0.9));
                    diasysArtery = M0_Systole_img - M0_Diastole_img;
                    RGBdiasys = labDuoImage(mean(M0_video,3), diasysArtery);
                    imwrite(RGBdiasys, fullfile(TB.path_main,'mask','DiaSysRGB.png'), 'png');
                catch

                    disp("Diasys png failed")


                end

                % try
                %     Commented until further fixes
                %     openmaskinpaintnet(fullfile(TB.path_main,'mask','M0.png'), fullfile(TB.path_main,'mask','DiaSysRGB.png'));
                % catch
                %     disp("paint.net macro failed")
                % end

            else

                fprintf(2, "no input loaded\n")

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
            app.EyeFlowUIFigure.Position = [100 100 640 421];
            app.EyeFlowUIFigure.Name = 'EyeFlow';
            app.EyeFlowUIFigure.Icon = fullfile(pathToMLAPP, 'eyeflow_logo.png');

            % Create a grid layout to manage resizing
            grid = uigridlayout(app.EyeFlowUIFigure);
            grid.RowHeight = {'fit', '1x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            grid.ColumnWidth = {'1x', '1x', '1x', '1x', '3x'};
            grid.BackgroundColor = [0.149, 0.149, 0.149];

            % Top Row: Load Folder, Load Holo, Clear, Folder Management
            app.LoadfolderButton = uibutton(grid, 'push');
            app.LoadfolderButton.ButtonPushedFcn = createCallbackFcn(app, @LoadfolderButtonPushed, true);
            app.LoadfolderButton.BackgroundColor = [0.502 0.502 0.502];
            app.LoadfolderButton.FontSize = 16;
            app.LoadfolderButton.FontColor = [0.9412 0.9412 0.9412];
            app.LoadfolderButton.Layout.Row = 1;
            app.LoadfolderButton.Layout.Column = 1;
            app.LoadfolderButton.Text = 'Load Folder';

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
            app.ReferenceDirectory.FontSize = 16;
            app.ReferenceDirectory.FontColor = [0.9412 0.9412 0.9412];
            app.ReferenceDirectory.BackgroundColor = [0.149 0.149 0.149];
            app.ReferenceDirectory.Layout.Row = 1;
            app.ReferenceDirectory.Layout.Column = 1;

            % Create Lamp
            app.Lamp = uilamp(dirgrid);
            app.Lamp.Layout.Row = 1;
            app.Lamp.Layout.Column = 2;

            % Third Row: Edit Parameters, Edit Masks, Play Inputs, Preview Masks
            app.EditParametersButton = uibutton(grid, 'push');
            app.EditParametersButton.ButtonPushedFcn = createCallbackFcn(app, @EditParametersButtonPushed, true);
            app.EditParametersButton.BackgroundColor = [0.502 0.502 0.502];
            app.EditParametersButton.FontSize = 16;
            app.EditParametersButton.FontColor = [0.9412 0.9412 0.9412];
            app.EditParametersButton.Layout.Row = 3;
            app.EditParametersButton.Layout.Column = 1;
            app.EditParametersButton.Text = 'Edit Parameters';
            app.EditParametersButton.Tooltip = 'Find the eyeflow parameters here.';

            app.EditMasksButton = uibutton(grid, 'push');
            app.EditMasksButton.ButtonPushedFcn = createCallbackFcn(app, @EditMasksButtonPushed, true);
            app.EditMasksButton.BackgroundColor = [0.502 0.502 0.502];
            app.EditMasksButton.FontSize = 16;
            app.EditMasksButton.FontColor = [0.9412 0.9412 0.9412];
            app.EditMasksButton.Layout.Row = 3;
            app.EditMasksButton.Layout.Column = 2;
            app.EditMasksButton.Text = 'Edit Masks';
            app.EditMasksButton.Tooltip = 'Open mask folder and use forceMaskArtery.png and forceMaskVein.png to force the segmentation';

            app.PlayInputsButton = uibutton(grid, 'push');
            app.PlayInputsButton.ButtonPushedFcn = createCallbackFcn(app, @PlayInputsButtonPushed, true);
            app.PlayInputsButton.BackgroundColor = [0.502 0.502 0.502];
            app.PlayInputsButton.FontSize = 16;
            app.PlayInputsButton.FontColor = [0.9412 0.9412 0.9412];
            app.PlayInputsButton.Layout.Row = 3;
            app.PlayInputsButton.Layout.Column = 3;
            app.PlayInputsButton.Text = 'Play Inputs';

            app.PreviewMasksButton = uibutton(grid, 'push');
            app.PreviewMasksButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewMasksButtonPushed, true);
            app.PreviewMasksButton.BackgroundColor = [0.502 0.502 0.502];
            app.PreviewMasksButton.FontSize = 16;
            app.PreviewMasksButton.FontColor = [0.9412 0.9412 0.9412];
            app.PreviewMasksButton.Layout.Row = 3;
            app.PreviewMasksButton.Layout.Column = 4;
            app.PreviewMasksButton.Text = 'Preview Masks';

            % Checkboxes: Segmentation, Pulse analysis, Blood Flow Velocity, Blood Volume Rate, SH analysis
            app.SegmentationCheckBox = uicheckbox(grid);
            app.SegmentationCheckBox.Text = 'Segmentation';
            app.SegmentationCheckBox.FontSize = 16;
            app.SegmentationCheckBox.FontColor = [1 1 1];
            app.SegmentationCheckBox.Layout.Row = 4;
            app.SegmentationCheckBox.Layout.Column = [1, 4];
            app.SegmentationCheckBox.Value = true;

            app.PulseanalysisCheckBox = uicheckbox(grid);
            app.PulseanalysisCheckBox.Text = 'Pulse analysis';
            app.PulseanalysisCheckBox.FontSize = 16;
            app.PulseanalysisCheckBox.FontColor = [1 1 1];
            app.PulseanalysisCheckBox.Layout.Row = 5;
            app.PulseanalysisCheckBox.Layout.Column = [1, 4];
            app.PulseanalysisCheckBox.Value = true;

            app.velocityCheckBox = uicheckbox(grid);
            app.velocityCheckBox.Text = 'Blood Flow Velocity';
            app.velocityCheckBox.FontSize = 16;
            app.velocityCheckBox.FontColor = [1 1 1];
            app.velocityCheckBox.Layout.Row = 6;
            app.velocityCheckBox.Layout.Column = [1, 4];
            app.velocityCheckBox.Value = true;

            app.bloodVolumeRateCheckBox = uicheckbox(grid);
            app.bloodVolumeRateCheckBox.Text = 'Blood Volume Rate';
            app.bloodVolumeRateCheckBox.FontSize = 16;
            app.bloodVolumeRateCheckBox.FontColor = [1 1 1];
            app.bloodVolumeRateCheckBox.Layout.Row = 7;
            app.bloodVolumeRateCheckBox.Layout.Column = [1, 4];
            app.bloodVolumeRateCheckBox.Value = true;

            app.SHanalysisCheckBox = uicheckbox(grid);
            app.SHanalysisCheckBox.Text = 'SH analysis';
            app.SHanalysisCheckBox.FontSize = 16;
            app.SHanalysisCheckBox.FontColor = [1 1 1];
            app.SHanalysisCheckBox.Layout.Row = 8;
            app.SHanalysisCheckBox.Layout.Column = [1, 4];
            app.SHanalysisCheckBox.Enable = true;

            % Bottom Left: Execute Button
            app.ExecuteButton = uibutton(grid, 'push');
            app.ExecuteButton.ButtonPushedFcn = createCallbackFcn(app, @ExecuteButtonPushed, true);
            app.ExecuteButton.BackgroundColor = [0.502 0.502 0.502];
            app.ExecuteButton.FontSize = 16;
            app.ExecuteButton.FontColor = [0.9412 0.9412 0.9412];
            app.ExecuteButton.Enable = 'off';
            app.ExecuteButton.Layout.Row = 9;
            app.ExecuteButton.Layout.Column = 1;
            app.ExecuteButton.Text = 'Execute';

            % Number of Workers Spinner and Label
            app.NumberofWorkersSpinnerLabel = uilabel(grid);
            app.NumberofWorkersSpinnerLabel.HorizontalAlignment = 'right';
            app.NumberofWorkersSpinnerLabel.FontColor = [0.902 0.902 0.902];
            app.NumberofWorkersSpinnerLabel.FontSize = 16;
            app.NumberofWorkersSpinnerLabel.Layout.Row = 9;
            app.NumberofWorkersSpinnerLabel.Layout.Column = 3;
            app.NumberofWorkersSpinnerLabel.Text = 'Number of Workers';

            app.NumberofWorkersSpinner = uispinner(grid);
            app.NumberofWorkersSpinner.FontSize = 16;
            app.NumberofWorkersSpinner.Layout.Row = 9;
            app.NumberofWorkersSpinner.Layout.Column = 4;
            maxWorkers = parcluster('local').NumWorkers;
            app.NumberofWorkersSpinner.Limits = [0 maxWorkers];  % Ensure valid range
            app.NumberofWorkersSpinner.Value = min(10, floor(maxWorkers/2)); % Default to 10 or max available

            % Bottom Right: Overwrite Checkbox
            app.OverWriteCheckBox = uicheckbox(grid);
            app.OverWriteCheckBox.Text = 'Overwrite';
            app.OverWriteCheckBox.FontSize = 16;
            app.OverWriteCheckBox.FontColor = [0.8 0.8 0.8];
            app.OverWriteCheckBox.Layout.Row = 9;
            app.OverWriteCheckBox.Layout.Column = 2;
            app.OverWriteCheckBox.Value = false;
            app.OverWriteCheckBox.ValueChangedFcn = createCallbackFcn(app, @OverWriteCheckBoxChanged, true);
            app.OverWriteCheckBox.Tooltip = 'Overwrite the new results in the last EF_ result folder (to save space)';

            % Add a new column for the image
            app.ImageDisplay = uiimage(grid);
            app.ImageDisplay.Layout.Row = [1, 9]; % Span all rows
            app.ImageDisplay.Layout.Column = 5;   % Place in the new column
            app.ImageDisplay.ScaleMethod = 'stretch'; % Adjust the image to fit the component

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