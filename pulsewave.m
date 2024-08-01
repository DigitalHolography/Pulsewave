classdef pulsewave < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PulsewaveUIFigure             matlab.ui.Figure
        EditParametersButton          matlab.ui.control.Button
        NumberofWorkersSpinner        matlab.ui.control.Spinner
        NumberofWorkersSpinnerLabel   matlab.ui.control.Label
        ClearListButton               matlab.ui.control.Button
        LoadTextButton                matlab.ui.control.Button
        ExecutefromtextButton         matlab.ui.control.Button
        SHanalysisCheckBox            matlab.ui.control.CheckBox
        FolderManagementButton        matlab.ui.control.Button
        PulsewaveanalysisCheckBox     matlab.ui.control.CheckBox
        ReferenceDirectory            matlab.ui.control.TextArea
        ErrorLabel                    matlab.ui.control.Label
        Lamp                          matlab.ui.control.Lamp
        ApplyFlatFieldCheckBox        matlab.ui.control.CheckBox
        ClearButton                   matlab.ui.control.Button
        LoadfolderButton              matlab.ui.control.Button
        ExecuteButton                 matlab.ui.control.Button
        NumberofframesEditField       matlab.ui.control.NumericEditField
        NumberofframesEditFieldLabel  matlab.ui.control.Label
    end

    properties (Access = private)
        files = {}
        drawer_list = {}
        flag_is_load
    end

    methods (Access = private)
        function Load(app, path)
       
                for n = 1:length(app.files)
                    if length(app.files)~=app.files{n}.nbFiles 
                        app.ErrorLabel.Text = "Error: number of files not correct";
                        return ; 
                    end 
                end
                app.Lamp.Color = [1, 0, 0];
                drawnow;
                path = strcat(path, '\');
                
                parfor_arg = app.NumberofWorkersSpinner.Value ;

                poolobj = gcp('nocreate'); % check if a pool already exist
                if isempty(poolobj)
                    parpool(parfor_arg); % create a new pool
                elseif poolobj.NumWorkers ~= parfor_arg
                    delete(poolobj); %close the current pool to create a new one with correct num of workers
                    parpool(parfor_arg);
                end

                tic

                try
                    %% add files to the drawer list
                    app.files{end + 1} = OneCycleClass(path);

                    %% crop videos
                    for n = 1:length(app.files)
                        app.files{n} = app.files{n}.cropAllVideo();
                    end

                    %% moment normalize
                    for n = 1:length(app.files)
                        app.files{n} = app.files{n}.MomentNormalize();
                    end

                    %% Video resize (preprocess interpolation interpolate)
                    app.files{1} = app.files{1}.VideoResize();

                    %% interpolate
                    app.files{1} = app.files{1}.Interpolate();

                    %% End
                    app.LoadfolderButton.Enable = true ;
                    app.ExecuteButton.Enable = true ;
                    app.ClearButton.Enable = true ;
                    % FIXME app.ReferenceDirectory.Value = fullfile(path,file)
                    app.ReferenceDirectory.Value = path ;

                catch exception

                    fprintf("==============================\nERROR\n==============================\n")
                    fprintf('Error while loading : %s\n', path)
                    fprintf("%s\n",exception.identifier)
                    fprintf("%s\n",exception.message)
                    % for i = 1:size(exception.stack,1)
                    %     stack = sprintf('%s : %s, line : %d \n', exception.stack(i).file, exception.stack(i).name, exception.stack(i).line);
                    %     fprintf(stack);
                    % end

                    if exception.identifier == "MATLAB:audiovideo:VideoReader:FileNotFound"

                        fprintf("No Raw File was found, please check 'save raw files' in HoloDoppler")

                    else

                        for i = 1:numel(exception.stack)
                            disp(exception.stack(i))
                        end

                    end

                    fprintf("==============================\n")


                    app.Lamp.Color = [1, 1/2, 0];

                end

                disp('Load timing :')
                toc
                
        end
    end
methods (Access = private)

    function LoadFromTxt(app)

        [selected_file,path] = uigetfile('*.txt');
        if (selected_file)
            fileID = fopen(fullfile(path,selected_file),'r');
            % filecontent = fscanf(fileID,'%s');
            % files_cell = strsplit(filecontent,'\n');

            filecontent = fileread(fullfile(path,selected_file));
            files_cell = strsplit(filecontent,'\n');
            for nn = 1:length(files_cell)
                if ~isempty(files_cell{1,nn})
                    app.drawer_list{end + 1} = files_cell{1,nn};
                end
            end
            fclose(fileID);
        end

    end
end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            displaySplashScreen();
        end

        % Button pushed function: LoadfolderButton
        function LoadfolderButtonPushed(app, event)
            % clearing before loading
            app.files = {};
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;
            
            clear Parameters_json
            if (app.flag_is_load)
                disp("Files already loaded")
                app.ClearButton.Enable = true ;
            else
                f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
                selected_dir = uigetdir();
                delete(f); %delete the dummy figure
                app.flag_is_load = true;
                app.Load(selected_dir);
                
            end

            app.ErrorLabel.Text = "" ;
            app.Lamp.Color = [0, 1, 0];
        end

        % Button pushed function: ExecuteButton
        function ExecuteButtonPushed(app, event)
           
            %             delete(gcp('nocreate')); % closeparallel CPU pool
            %             parpool;% launch parallel CPU pool

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

            for n = 1:length(app.files)
              
                app.files{n}.flag_FlatField = app.ApplyFlatFieldCheckBox.Value ;
                app.files{n}.flag_SH_analysis = app.SHanalysisCheckBox.Value;
                app.files{n}.flag_PulseWave_analysis = app.PulsewaveanalysisCheckBox.Value;
                try
                    app.files{n}.onePulse(app.NumberofframesEditField.Value);
                catch exception
                    fprintf("==============================\nERROR\n==============================\n")
                    disp(['Error with file : ', app.files{n}.directory])
                    disp(exception.identifier)
                    disp(exception.message)
                    for i = 1:size(exception.stack,1)
                        fprintf('%s : %s, line : %d \n',exception.stack(i).file, exception.stack(i).name, exception.stack(i).line);
                    end
                    fprintf("==============================\n")
                end
            end
            app.Lamp.Color = [0, 1, 0];
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            app.files = {};
            app.ReferenceDirectory.Value = "";
            app.LoadfolderButton.Enable = true;
            app.flag_is_load = false;
            
            clear Parameters_json
            
        end

        % Button pushed function: FolderManagementButton
        function FolderManagementButtonPushed(app, event)
            d = dialog('Position', [300, 300, 750, 90 + length(app.drawer_list) * 14],...
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
                'String', 'Export to text',...
                'Callback', @export_text);


            uicontrol('Parent', d,...
                'Position', [500, 20, 100, 25],...
                'FontName', 'Helvetica',...
                'BackgroundColor', [0.5, 0.5, 0.5],...
                'ForegroundColor', [0.9 0.9 0.9],...
                'FontWeight', 'bold',...
                'String', 'Render',...
                'Callback', @render);


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
                [~,folder_name,~] = fileparts(selected_dir); 
                % List of Subfolders within the measurement folder
                tmp_dir = dir(selected_dir);
                % remove all files (isdir property is 0)
                subfoldersName = tmp_dir([tmp_dir(:).isdir]);
                % remove '.' and '..' folders
                subfoldersName = subfoldersName(~ismember({subfoldersName(:).name},{'.','..'}));
                subfoldersName = {subfoldersName.name};
                % remove of other folders (ex: 'config' subfolders)
                for ii=1:length(subfoldersName)
                    if contains(subfoldersName{ii},folder_name)
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

            function export_text(~, ~)
                selected_dir = uigetdir('C:\', 'Select file location');
                filename = fullfile(selected_dir, 'folderManagement.txt');
                fileID = fopen(filename,'a+') ;
                if ~isempty(app.drawer_list)
                    for i = 1:length(app.drawer_list)
                        fprintf(fileID,'%s \n',app.drawer_list{i});
                    end
                end
                
                fclose(fileID);
            end

            function render(~, ~)
                for i = 1:length(app.drawer_list)
                    tic
                    %app.PulsewaveanalysisCheckBox.Value = true;
                    app.Load(app.drawer_list{i});
                    app.ExecuteButtonPushed();
                    app.ClearButtonPushed();
                    toc
                end
%                 clear app.drawer_list
            end
            delete(d);
        end

        % Button pushed function: ExecutefromtextButton
        function ExecutefromtextButtonPushed(app, event)
            app.ClearButtonPushed();
            path_to_copy = uigetdir();
            if ~isempty(app.drawer_list) && ~isempty(path_to_copy)
                for i = 1:length(app.drawer_list)
                    app.Load(app.drawer_list{i});
                    app.ExecuteButtonPushed();
                    mkdir(fullfile(path_to_copy,app.files{1}.ToolBoxmaster.PW_folder_name))
                    copyfile(app.files{1}.ToolBoxmaster.PW_path_dir,fullfile(path_to_copy,app.files{1}.ToolBoxmaster.PW_folder_name))
                    app.ClearButtonPushed();
                end
            end

            % clear app.drawer_list
        end

        % Button pushed function: LoadTextButton
        function LoadTextButtonPushed(app, event)
            app.ClearButtonPushed();
            app.ClearListButtonPushed()
            app.LoadFromTxt()
            app.ClearButton.Enable = true;
        end

        % Button pushed function: ClearListButton
        function ClearListButtonPushed(app, event)
            app.drawer_list = {};
        end

        % Button pushed function: EditParametersButton
        function EditParametersButtonPushed(app, event)
            if (app.flag_is_load)
                winopen(fullfile(app.files{1}.ToolBoxmaster.PW_path_main,'json','InputPulsewaveParams.json'));
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PulsewaveUIFigure and hide until all components are created
            app.PulsewaveUIFigure = uifigure('Visible', 'off');
            app.PulsewaveUIFigure.Color = [0.149 0.149 0.149];
            app.PulsewaveUIFigure.Position = [100 100 640 319];
            app.PulsewaveUIFigure.Name = 'Pulsewave';

            % Create NumberofframesEditFieldLabel
            app.NumberofframesEditFieldLabel = uilabel(app.PulsewaveUIFigure);
            app.NumberofframesEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofframesEditFieldLabel.FontSize = 16;
            app.NumberofframesEditFieldLabel.FontColor = [0.9412 0.9412 0.9412];
            app.NumberofframesEditFieldLabel.Position = [406 221 138 22];
            app.NumberofframesEditFieldLabel.Text = 'Number of frames ';

            % Create NumberofframesEditField
            app.NumberofframesEditField = uieditfield(app.PulsewaveUIFigure, 'numeric');
            app.NumberofframesEditField.Limits = [1 Inf];
            app.NumberofframesEditField.FontSize = 16;
            app.NumberofframesEditField.FontColor = [0.149 0.149 0.149];
            app.NumberofframesEditField.Position = [545 221 44 22];
            app.NumberofframesEditField.Value = 256;

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
            app.LoadfolderButton.Position = [61 218 123 28];
            app.LoadfolderButton.Text = 'Load folder';

            % Create ClearButton
            app.ClearButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.BackgroundColor = [0.502 0.502 0.502];
            app.ClearButton.FontSize = 16;
            app.ClearButton.FontColor = [0.9412 0.9412 0.9412];
            app.ClearButton.Enable = 'off';
            app.ClearButton.Position = [63 136 100 27];
            app.ClearButton.Text = 'Clear';

            % Create ApplyFlatFieldCheckBox
            app.ApplyFlatFieldCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.ApplyFlatFieldCheckBox.Text = 'Apply Flat Field';
            app.ApplyFlatFieldCheckBox.FontSize = 16;
            app.ApplyFlatFieldCheckBox.FontColor = [0.9412 0.9412 0.9412];
            app.ApplyFlatFieldCheckBox.Position = [208 221 132 22];

            % Create Lamp
            app.Lamp = uilamp(app.PulsewaveUIFigure);
            app.Lamp.Position = [569 182 20 20];

            % Create ErrorLabel
            app.ErrorLabel = uilabel(app.PulsewaveUIFigure);
            app.ErrorLabel.HorizontalAlignment = 'center';
            app.ErrorLabel.FontSize = 16;
            app.ErrorLabel.FontColor = [1 0 0];
            app.ErrorLabel.Position = [61 268 542 23];
            app.ErrorLabel.Text = '';

            % Create ReferenceDirectory
            app.ReferenceDirectory = uitextarea(app.PulsewaveUIFigure);
            app.ReferenceDirectory.FontSize = 16;
            app.ReferenceDirectory.FontColor = [0.9412 0.9412 0.9412];
            app.ReferenceDirectory.BackgroundColor = [0.149 0.149 0.149];
            app.ReferenceDirectory.Position = [61 180 485 24];

            % Create PulsewaveanalysisCheckBox
            app.PulsewaveanalysisCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.PulsewaveanalysisCheckBox.Text = 'Pulse wave analysis';
            app.PulsewaveanalysisCheckBox.FontSize = 16;
            app.PulsewaveanalysisCheckBox.FontColor = [1 1 1];
            app.PulsewaveanalysisCheckBox.Position = [63 96 199 24];

            % Create FolderManagementButton
            app.FolderManagementButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.FolderManagementButton.ButtonPushedFcn = createCallbackFcn(app, @FolderManagementButtonPushed, true);
            app.FolderManagementButton.BackgroundColor = [0.502 0.502 0.502];
            app.FolderManagementButton.FontSize = 16;
            app.FolderManagementButton.FontColor = [0.9412 0.9412 0.9412];
            app.FolderManagementButton.Position = [431 134 158 28];
            app.FolderManagementButton.Text = 'Folder Management';

            % Create SHanalysisCheckBox
            app.SHanalysisCheckBox = uicheckbox(app.PulsewaveUIFigure);
            app.SHanalysisCheckBox.Text = 'SH analysis';
            app.SHanalysisCheckBox.FontSize = 16;
            app.SHanalysisCheckBox.FontColor = [1 1 1];
            app.SHanalysisCheckBox.Position = [63 62 146 25];

            % Create ExecutefromtextButton
            app.ExecutefromtextButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.ExecutefromtextButton.ButtonPushedFcn = createCallbackFcn(app, @ExecutefromtextButtonPushed, true);
            app.ExecutefromtextButton.BackgroundColor = [0.502 0.502 0.502];
            app.ExecutefromtextButton.FontSize = 16;
            app.ExecutefromtextButton.FontColor = [1 1 1];
            app.ExecutefromtextButton.Position = [431 24 163 32];
            app.ExecutefromtextButton.Text = 'Execute from text';

            % Create LoadTextButton
            app.LoadTextButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.LoadTextButton.ButtonPushedFcn = createCallbackFcn(app, @LoadTextButtonPushed, true);
            app.LoadTextButton.BackgroundColor = [0.502 0.502 0.502];
            app.LoadTextButton.FontSize = 16;
            app.LoadTextButton.FontColor = [1 1 1];
            app.LoadTextButton.Position = [347 70 136 28];
            app.LoadTextButton.Text = 'Load Text';

            % Create ClearListButton
            app.ClearListButton = uibutton(app.PulsewaveUIFigure, 'push');
            app.ClearListButton.ButtonPushedFcn = createCallbackFcn(app, @ClearListButtonPushed, true);
            app.ClearListButton.BackgroundColor = [0.502 0.502 0.502];
            app.ClearListButton.FontSize = 16;
            app.ClearListButton.FontColor = [1 1 1];
            app.ClearListButton.Position = [502 70 87 28];
            app.ClearListButton.Text = 'Clear List';

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
            app.EditParametersButton.Position = [276 134 130 28];
            app.EditParametersButton.Text = 'Edit Parameters';

            % Show the figure after all components are created
            app.PulsewaveUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = pulsewave

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