classdef pulse < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PulsewaveUIFigure             matlab.ui.Figure
        PreProcessButton              matlab.ui.control.Button
        PlayInputsButton              matlab.ui.control.Button
        EditParametersButton          matlab.ui.control.Button
        EditMasksButton               matlab.ui.control.Button
        NumberofWorkersSpinner        matlab.ui.control.Spinner
        NumberofWorkersSpinnerLabel   matlab.ui.control.Label
        SHanalysisCheckBox            matlab.ui.control.CheckBox
        FolderManagementButton        matlab.ui.control.Button
        PreviewMasksButton            matlab.ui.control.Button
        SegmentationCheckBox          matlab.ui.control.CheckBox
        PulsewaveanalysisCheckBox     matlab.ui.control.CheckBox
        ExtendedPulsewaveCheckBox     matlab.ui.control.CheckBox
        velocityCheckBox              matlab.ui.control.CheckBox
        bloodVolumeRateCheckBox       matlab.ui.control.CheckBox
        bloodVelocityProfileCheckBox  matlab.ui.control.CheckBox
        ReferenceDirectory            matlab.ui.control.TextArea
        OverWriteCheckBox             matlab.ui.control.CheckBox
        ErrorLabel                    matlab.ui.control.Label
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
                fprintf("\n----------------------------------\n")
                fprintf("Video Loading\n")
                fprintf("----------------------------------\n")
                app.file = OneCycleClass(path);
                fprintf("- Video Loading took : %ds\n", round(toc))

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
                    "Welcome to Pulsewave %s\n" + ...
                    "------------------------------------------\n" + ...
                    "Developed by the DigitalHolographyFoundation\n" + ...
                    "==========================================\n", v(1));
            end

            % Add necessary paths
            addpath("BloodFlowVelocity\", "BloodFlowVelocity\Elastography\", "BloodVolumeRate\", ...
                "BloodVolumeRate\Rheology\", "Loading\", "Parameters\", "Preprocessing\", ...
                "PulseAnalysis\", "Scripts\", "Segmentation\", "SHAnalysis\", "Tools\");

            % Set the UI title
            app.PulsewaveUIFigure.Name = ['Pulsewave ', char(v(1))];

            % Auto-detect the maximum available workers and update the spinner
            maxWorkers = parcluster('local').NumWorkers;
            app.NumberofWorkersSpinner.Limits = [0 maxWorkers];  % Ensure valid range
            app.NumberofWorkersSpinner.Value = min(10, floor(maxWorkers/2)); % Default to 10 or max available

            % Initialize checkboxes and flags
            app.updateCheckboxes();
            app.flag_is_load = false;

            % Display splash screen
            displaySplashScreen();
        end

        % Button pushed function: LoadfolderButton
        function LoadfolderButtonPushed(app, ~)
            % clearing before loading
            if ~isempty(app.file)
                last_dir = app.file.directory;
            else
                last_dir = [];
            end

            app.file = [];
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;

            clear Parameters_json
            if (app.flag_is_load)
                disp("Files already loaded")
                app.ClearButton.Enable = true ;
                app.EditParametersButton.Enable = true;
            else
                f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
                selected_dir = uigetdir(last_dir);
                if selected_dir == 0
                    disp('No folder selected')
                    return
                end
                delete(f); %delete the dummy figure
                app.flag_is_load = true;
                app.Load(selected_dir);
            end

        end
        function LoadHoloButtonPushed(app, ~)
            % clearing before loading

            app.file =  [];
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;

            clear Parameters_json
            if (app.flag_is_load)
                disp("Files already loaded")
                app.ClearButton.Enable = true ;
                app.EditParametersButton.Enable = true;
            else
                f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
                [selected_holo,path_holo] = uigetfile('*.holo');
                if selected_holo == 0
                    disp('No file selected')
                    return
                end
                delete(f); %delete the dummy figure
                app.flag_is_load = true;
                app.Load(fullfile(path_holo,selected_holo));

            end

            app.ErrorLabel.Text = "" ;
            app.Lamp.Color = [0, 1, 0];
        end

        % Button pushed function: ExecuteButton
        function ExecuteButtonPushed(app, ~)
            if ~app.flag_is_load
                fprintf(2, "no input loaded\n")
                return
            end

            if ~app.file.is_preprocessed
                app.PreProcessButtonPushed();
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
            app.ErrorLabel.Text = "" ;
            drawnow;

            % Actualizes the input Parameters
            app.file.PW_params_names = checkPulsewaveParamsFromJson(app.file.directory); % checks compatibility between found PW params and Default PW params of this version of PW.


            for i = 1:length(app.file.PW_params_names)

                app.file.PW_param_name = app.file.PW_params_names{i};

                fprintf("==========================================\n")
                app.file.flag_Segmentation = app.SegmentationCheckBox.Value;
                app.file.flag_SH_analysis = app.SHanalysisCheckBox.Value;
                app.file.flag_PulseWave_analysis = app.PulsewaveanalysisCheckBox.Value;
                app.file.flag_velocity_analysis = app.velocityCheckBox.Value;
                app.file.flag_ExtendedPulseWave_analysis = app.ExtendedPulsewaveCheckBox.Value;
                app.file.flag_bloodVolumeRate_analysis = app.bloodVolumeRateCheckBox.Value;
                app.file.flag_bloodVelocityProfile_analysis = app.bloodVelocityProfileCheckBox.Value;

                app.file.OverWrite = app.OverWriteCheckBox.Value;

                try
                    app.file = app.file.onePulse();

                catch ME

                    diary off

                    fprintf(2,"==========================================\n");
                    fprintf(2,"ERROR\n");
                    fprintf(2,"==========================================\n");

                    fprintf(2, 'Error with file : %s\n%s\n%s', app.file.directory, ME.identifier, ME.message);

                    for stackIdx = 1:size(ME.stack, 1)
                        fprintf(2,"%s : %s, line : %d\n", ME.stack(stackIdx).file, ME.stack(stackIdx).name, ME.stack(stackIdx).line);
                    end

                    fprintf(2,"==========================================\n");
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

        function PreProcessButtonPushed(app, ~)
            if ~app.flag_is_load
                fprintf(2, "no input loaded\n")
                return
            end
            if app.file.is_preprocessed
                disp('input already preprocessed. reload if needed.')
                return
            end
            app.Lamp.Color = [1, 0, 0];
            drawnow;

            parfor_arg = app.NumberofWorkersSpinner.Value ;
            poolobj = gcp('nocreate'); % check if a pool already exist
            if isempty(poolobj)
                parpool(parfor_arg); % create a new pool
            elseif poolobj.NumWorkers ~= parfor_arg
                delete(poolobj); %close the current pool to create a new one with correct num of workers
                parpool(parfor_arg);
            end

            totalPreProcessTime = tic;

            try
                fprintf("\n----------------------------------\n")
                fprintf("Video PreProcessing\n")
                fprintf("----------------------------------\n")
                app.file = app.file.preprocessData();
                app.Lamp.Color = [0, 1, 0];
            catch ME

                fprintf("==========================================\nERROR\n==========================================\n")
                if ~isempty(app.file)
                    fprintf('Error while preprocessing : %s\n', app.file.directory)
                else
                    fprintf('Error while preprocessing : %s\n', 'xx')
                end
                fprintf("%s\n",ME.identifier)
                fprintf("%s\n",ME.message)

                for i = 1:size(ME.stack,1)
                    fprintf('%s : %s, line : %d \n', ME.stack(i).file, ME.stack(i).name, ME.stack(i).line);
                end

                fprintf("==========================================\n")
                app.Lamp.Color = [1, 1/2, 0];

            end

            fprintf("----------------------------------\n")
            fprintf("- Total PreProcess timing took : %ds\n", round(toc(totalPreProcessTime)))

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


            function import_param(~, ~)
                tic
                f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
                [selected_json,path_json] = uigetfile('*.json');
                if selected_json == 0
                    disp('No file selected')
                    return
                end
                delete(f); %delete the dummy figure

                for ind = 1:length(app.drawer_list)
                    pw_path_json = fullfile(app.drawer_list{ind},'pulsewave','json');
                    if ~isfolder(pw_path_json)
                        mkdir(pw_path_json);
                    end
                    copyfile(fullfile(path_json,selected_json),pw_path_json);

                    %get idx for renaming
                    idx = 0;
                    list_dir = dir(pw_path_json);
                    for i = 1:numel(list_dir)
                        match = regexp(list_dir(i).name, '\d+$', 'match');
                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end
                    end

                    %renaming
                    copyfile(fullfile(pw_path_json,selected_json),fullfile(pw_path_json,sprintf('InputPulseWaveParams_%d.json',idx)));
                    delete(fullfile(pw_path_json,selected_json));
                end
                toc
            end

            function render(~, ~)
                for i = 1:length(app.drawer_list)
                    tic
                    app.Load(app.drawer_list{i});
                    app.PreProcessButtonPushed();
                    app.ExecuteButtonPushed();
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
            PreviewMasks(app);
        end


        % Checkbox update function:
        function updateCheckboxes(app, ~)
            % if not(isempty(app.files)) && not(isempty(app.files{end}.maskArtery)) % if segmentation masks exists
            %     app.PulsewaveanalysisCheckBox.Enable = true;
            %     if not(isempty(app.files)) && not(isempty(app.files{end}.vRMS)) % if velocity estimate exists
            %
            %         app.ExtendedPulsewaveCheckBox.Enable = true;
            %         app.velocityCheckBox.Enable = true;
            %         app.bloodVolumeRateCheckBox.Enable = true;
            %         app.bloodVelocityProfileCheckBox.Enable = true;
            %     else
            %         app.ExtendedPulsewaveCheckBox.Enable = false;
            %         app.velocityCheckBox.Enable = false;
            %         app.bloodVolumeRateCheckBox.Enable = false;
            %         app.bloodVelocityProfileCheckBox.Enable = false;
            %     end
            %
            % else
            %     app.PulsewaveanalysisCheckBox.Enable = false;
            %     app.ExtendedPulsewaveCheckBox.Enable = false;
            %     app.velocityCheckBox.Enable = false;
            %     app.bloodVolumeRateCheckBox.Enable = false;
            %     app.bloodVelocityProfileCheckBox.Enable = false;
            % end
        end

        % Button pushed function: EditParametersButton
        function EditParametersButtonPushed(app, ~)
            if (app.flag_is_load)
                main_path = fullfile(app.file.directory, 'pulsewave');
                if isfile(fullfile(main_path,'json',app.file.PW_param_name))
                    disp(['opening : ', fullfile(main_path,'json',app.file.PW_param_name)])
                    winopen(fullfile(main_path,'json',app.file.PW_param_name));
                else
                    disp(['couldn''t open : ',fullfile(main_path,'json',app.file.PW_param_name)])
                end
            else
                fprintf(2, "no input loaded\n")
            end
        end

        % Button pushed function: EditMasksButton
        function EditMasksButtonPushed(app, ~)
            ToolBox = getGlobalToolBox;

            if (app.flag_is_load)
                if ~isfolder(fullfile(ToolBox.PW_path_main,'mask'))
                    mkdir(fullfile(ToolBox.PW_path_main,'mask'))
                end
                try
                    winopen(fullfile(ToolBox.PW_path_main,'mask'));
                catch
                    disp("opening failed.")
                end
                try
                    list_dir = dir(ToolBox.PW_path_main);
                    idx=0;
                    for i=1:length(list_dir)
                        if contains(list_dir(i).name, ToolBox.PW_name)
                            match = regexp(list_dir(i).name, '\d+$', 'match');
                            if ~isempty(match) && str2double(match{1}) >= idx
                                idx = str2double(match{1}) ; %suffix
                            end
                        end
                    end
                    PW_path_dir = fullfile(ToolBox.PW_path_main, ToolBox.PW_folder_name);

                    disp(['Copying from : ',fullfile(PW_path_dir,'png','mask')])
                    copyfile(fullfile(PW_path_dir,'png', 'mask', sprintf("%s_maskArtery.png", ToolBox.main_foldername)), fullfile(ToolBox.PW_path_main, 'mask', 'MaskArtery.png'));
                    copyfile(fullfile(PW_path_dir,'png', 'mask', sprintf("%s_maskVein.png", ToolBox.main_foldername)), fullfile(ToolBox.PW_path_main, 'mask', 'MaskVein.png'));
                catch
                    disp("last auto mask copying failed.")
                end
                try

                    copyfile(fullfile(ToolBox.PW_path,'png',sprintf("%s_M0.png",ToolBox.main_foldername)),fullfile(ToolBox.PW_path_main,'mask','M0.png'));
                    PW_folder_name = strcat(ToolBox.main_foldername, '_PW');
                    list_dir = dir(ToolBox.PW_path_main);
                    idx=0;
                    for i=1:length(list_dir)
                        if contains(list_dir(i).name, PW_folder_name)
                            match = regexp(list_dir(i).name, '\d+$', 'match');
                            if ~isempty(match) && str2double(match{1}) >= idx
                                idx = str2double(match{1}) ; %suffix
                            end
                        end
                    end
                    PW_folder_name = sprintf('%s_%d', PW_folder_name, idx);
                    copyfile(fullfile(PW_path_dir,'gif',sprintf("%s_M0.gif",PW_folder_name)),fullfile(ToolBox.PW_path_main,'mask','M0.gif'));
                catch

                    disp("last M0 png and gif copying failed")
                end
            else

                fprintf(2, "no input loaded\n")

            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create PulsewaveUIFigure and hide until all components are created
            app.PulsewaveUIFigure = uifigure('Visible', 'off');
            app.PulsewaveUIFigure.Color = [0.149 0.149 0.149];
            app.PulsewaveUIFigure.Position = [100 100 640 421];
            app.PulsewaveUIFigure.Name = 'Pulsewave';
            app.PulsewaveUIFigure.Icon = fullfile(pathToMLAPP, 'pulsewave_logo_temp.png');

            % Create ExecuteButton
            app.ExecuteButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.ExecuteButton.ButtonPushedFcn = createCallbackFcn(app, @ExecuteButtonPushed, true);
            app.ExecuteButton.BackgroundColor = [0.502 0.502 0.502];
            app.ExecuteButton.FontSize = 16;
            app.ExecuteButton.FontColor = [0.9412 0.9412 0.9412];
            app.ExecuteButton.Enable = 'off';
            app.ExecuteButton.Position = [61 24 100 27];
            app.ExecuteButton.Text = 'Execute';

            % Create LoadfolderButton
            app.LoadfolderButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.LoadfolderButton.ButtonPushedFcn = createCallbackFcn(app, @LoadfolderButtonPushed, true);
            app.LoadfolderButton.BackgroundColor = [0.502 0.502 0.502];
            app.LoadfolderButton.FontSize = 16;
            app.LoadfolderButton.FontColor = [0.9412 0.9412 0.9412];
            app.LoadfolderButton.Position = [61 322 123 28];
            app.LoadfolderButton.Text = 'Load folder';

            % Create LoadHoloButton
            app.LoadHoloButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.LoadHoloButton.ButtonPushedFcn = createCallbackFcn(app, @LoadHoloButtonPushed, true);
            app.LoadHoloButton.BackgroundColor = [0.502 0.502 0.502];
            app.LoadHoloButton.FontSize = 16;
            app.LoadHoloButton.FontColor = [0.9412 0.9412 0.9412];
            app.LoadHoloButton.Position = [61 362 123 28];
            app.LoadHoloButton.Text = 'Load holo';

            % Create ClearButton
            app.ClearButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.BackgroundColor = [0.502 0.502 0.502];
            app.ClearButton.FontSize = 16;
            app.ClearButton.FontColor = [0.9412 0.9412 0.9412];
            app.ClearButton.Enable = 'off';
            app.ClearButton.Position = [63 240 100 27];
            app.ClearButton.Text = 'Clear';

            % Create Lamp
            app.Lamp = uilamp(app.PulsewaveUIFigure);
            app.Lamp.Position = [569 286 20 20];

            % Create ErrorLabel
            app.ErrorLabel = uilabel(app.PulsewaveUIFigure);
            app.ErrorLabel.HorizontalAlignment = 'center';
            app.ErrorLabel.FontSize = 16;
            app.ErrorLabel.FontColor = [1 0 0];
            app.ErrorLabel.Position = [61 372 542 23];
            app.ErrorLabel.Text = '';

            % Create ReferenceDirectory
            app.ReferenceDirectory = uitextarea(app.PulsewaveUIFigure);
            app.ReferenceDirectory.FontSize = 16;
            app.ReferenceDirectory.FontColor = [0.9412 0.9412 0.9412];
            app.ReferenceDirectory.BackgroundColor = [0.149 0.149 0.149];
            app.ReferenceDirectory.Position = [61 284 485 24];

            % Create SegmentationCheckBox
            app.SegmentationCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.SegmentationCheckBox.Text = 'Segmentation';
            app.SegmentationCheckBox.FontSize = 16;
            app.SegmentationCheckBox.FontColor = [1 1 1];
            app.SegmentationCheckBox.Position = [63 198 250 24];
            app.SegmentationCheckBox.Value = true;
            app.SegmentationCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);

            % Create PulsewaveanalysisCheckBox
            app.PulsewaveanalysisCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.PulsewaveanalysisCheckBox.Text = 'Pulse wave analysis';
            app.PulsewaveanalysisCheckBox.FontSize = 16;
            app.PulsewaveanalysisCheckBox.FontColor = [1 1 1];
            app.PulsewaveanalysisCheckBox.Position = [63 164 250 24];
            app.PulsewaveanalysisCheckBox.Value = true;
            app.PulsewaveanalysisCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);

            % Create velocityCheckBox
            app.velocityCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.velocityCheckBox.Text = 'Blood Flow Velocity';
            app.velocityCheckBox.FontSize = 16;
            app.velocityCheckBox.FontColor = [1 1 1];
            app.velocityCheckBox.Position = [63 130 250 24];
            app.velocityCheckBox.Value = true;
            app.velocityCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);

            % Create ExtendedPulsewaveCheckBox
            app.ExtendedPulsewaveCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.ExtendedPulsewaveCheckBox.Text = 'Extended Pulse Analysis';
            app.ExtendedPulsewaveCheckBox.FontSize = 16;
            app.ExtendedPulsewaveCheckBox.FontColor = [1 1 1];
            app.ExtendedPulsewaveCheckBox.Position = [250 164 250 24];
            app.ExtendedPulsewaveCheckBox.Value = false;
            app.ExtendedPulsewaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);

            % Create bloodVolumeRateCheckBox
            app.bloodVolumeRateCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.bloodVolumeRateCheckBox.Text = 'Blood Volume Rate';
            app.bloodVolumeRateCheckBox.FontSize = 16;
            app.bloodVolumeRateCheckBox.FontColor = [1 1 1];
            app.bloodVolumeRateCheckBox.Position = [63 96 250 24];
            app.bloodVolumeRateCheckBox.Value = true;
            app.bloodVolumeRateCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);

            % Create bloodVelocityProfileCheckBox
            app.bloodVelocityProfileCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.bloodVelocityProfileCheckBox.Text = 'Blood Velocity Profile';
            app.bloodVelocityProfileCheckBox.FontSize = 16;
            app.bloodVelocityProfileCheckBox.FontColor = [1 1 1];
            app.bloodVelocityProfileCheckBox.Position = [250 96 250 24];
            app.bloodVelocityProfileCheckBox.Value = false;
            app.bloodVolumeRateCheckBox.ValueChangedFcn = createCallbackFcn(app, @updateCheckboxes, true);


            % Create OverWriteCheckBox
            app.OverWriteCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.OverWriteCheckBox.Text = 'over write';
            app.OverWriteCheckBox.FontSize = 16;
            app.OverWriteCheckBox.FontColor = [0.8 0.8 0.8];
            app.OverWriteCheckBox.Position = [500 24 130 28];
            app.OverWriteCheckBox.Value = false;
            app.OverWriteCheckBox.ValueChangedFcn = createCallbackFcn(app, @OverWriteCheckBoxChanged, true);
            app.OverWriteCheckBox.Tooltip = 'OverWrite the new results in the last PW_ result folder (to save space)';


            % Create PreviewMasksButton
            app.PreviewMasksButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.PreviewMasksButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewMasksButtonPushed, true);
            app.PreviewMasksButton.BackgroundColor = [0.502 0.502 0.502];
            app.PreviewMasksButton.FontSize = 16;
            app.PreviewMasksButton.FontColor = [0.9412 0.9412 0.9412];
            app.PreviewMasksButton.Position = [400 200 130 28];
            app.PreviewMasksButton.Text = 'Preview Masks';


            % Create FolderManagementButton
            app.FolderManagementButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.FolderManagementButton.ButtonPushedFcn = createCallbackFcn(app, @FolderManagementButtonPushed, true);
            app.FolderManagementButton.BackgroundColor = [0.502 0.502 0.502];
            app.FolderManagementButton.FontSize = 16;
            app.FolderManagementButton.FontColor = [0.9412 0.9412 0.9412];
            app.FolderManagementButton.Position = [431 238 158 28];
            app.FolderManagementButton.Text = 'Folder Management';

            % Create SHanalysisCheckBox
            app.SHanalysisCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.SHanalysisCheckBox.Text = 'SH analysis';
            app.SHanalysisCheckBox.FontSize = 16;
            app.SHanalysisCheckBox.FontColor = [1 1 1];
            app.SHanalysisCheckBox.Position = [63 62 250 24];
            app.SHanalysisCheckBox.Enable = true;

            % Create NumberofWorkersSpinnerLabel
            app.NumberofWorkersSpinnerLabel = uilabel(app.PulsewaveUIFigure);
            app.NumberofWorkersSpinnerLabel.HorizontalAlignment = 'right';
            app.NumberofWorkersSpinnerLabel.FontColor = [0.902 0.902 0.902];
            app.NumberofWorkersSpinnerLabel.Position = [169 26 109 22];
            app.NumberofWorkersSpinnerLabel.Text = 'Number of Workers';

            % Create NumberofWorkersSpinner
            app.NumberofWorkersSpinner = uispinner(app.PulsewaveUIFigure);
            app.NumberofWorkersSpinner.Limits = [-1 32];
            app.NumberofWorkersSpinner.Position = [289 26 51 22];
            app.NumberofWorkersSpinner.Value = 8;

            % Create EditParametersButton
            app.EditParametersButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.EditParametersButton.ButtonPushedFcn = createCallbackFcn(app, @EditParametersButtonPushed, true);
            app.EditParametersButton.BackgroundColor = [0.502 0.502 0.502];
            app.EditParametersButton.FontSize = 16;
            app.EditParametersButton.FontColor = [0.9412 0.9412 0.9412];
            app.EditParametersButton.Position = [276 238 130 28];
            app.EditParametersButton.Text = 'Edit Parameters';
            app.EditParametersButton.Enable = 'on';
            app.EditParametersButton.Tooltip = 'Find the pulse wave parameters here.';

            % Create PreProcessButton
            app.PreProcessButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.PreProcessButton.ButtonPushedFcn = createCallbackFcn(app, @PreProcessButtonPushed, true);
            app.PreProcessButton.BackgroundColor = [0.502 0.502 0.502];
            app.PreProcessButton.FontSize = 16;
            app.PreProcessButton.FontColor = [0.9412 0.9412 0.9412];
            app.PreProcessButton.Position = [191 322 130 28];
            app.PreProcessButton.Text = 'Pre Process';
            app.PreProcessButton.Enable = 'on';

            % Create PlayInputsButton
            app.PlayInputsButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.PlayInputsButton.ButtonPushedFcn = createCallbackFcn(app, @PlayInputsButtonPushed, true);
            app.PlayInputsButton.BackgroundColor = [0.502 0.502 0.502];
            app.PlayInputsButton.FontSize = 16;
            app.PlayInputsButton.FontColor = [0.9412 0.9412 0.9412];
            app.PlayInputsButton.Position = [351 322 130 28];
            app.PlayInputsButton.Text = 'Play Inputs';
            app.PlayInputsButton.Enable = 'on';

            % Create EditMasksButton
            app.EditMasksButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.EditMasksButton.ButtonPushedFcn = createCallbackFcn(app, @EditMasksButtonPushed, true);
            app.EditMasksButton.BackgroundColor = [0.502 0.502 0.502];
            app.EditMasksButton.FontSize = 16;
            app.EditMasksButton.FontColor = [0.9412 0.9412 0.9412];
            app.EditMasksButton.Position = [250 200 130 28];
            app.EditMasksButton.Text = 'Edit Masks';
            app.EditMasksButton.Enable = 'on';
            app.EditMasksButton.Tooltip = 'Open mask folder and use forceMaskArtery.png and forceMaskVein.png to force the segmentation';

            % Show the figure after all components are created
            app.PulsewaveUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = pulse

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PulsewaveUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PulsewaveUIFigure)
        end
    end
end