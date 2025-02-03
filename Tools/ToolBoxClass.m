classdef ToolBoxClass < handle
    % ToolBoxClass holds useful variables for pulsewave processing.

    properties
        % Paths
        PW_path char
        PW_path_main char
        PW_path_dir char
        PW_path_png char
        PW_path_eps char
        PW_path_gif char
        PW_path_txt char
        PW_path_avi char
        PW_path_mp4 char
        PW_path_json char
        PW_path_log char
        main_foldername char
        PW_param_name char
        PW_folder_name char
        % Parameters
        stride double
        fs double
        f1 double
        f2 double
        ScalingFactorVelocityInPlane double
    end

    methods

        function obj = ToolBoxClass(path, PW_param_name, OverWrite)
            % Constructor for ToolBoxClass: Initializes paths, parameters, and calculates scaling factors.

            % Store paths and parameters
            obj.PW_path = path;
            obj.PW_param_name = PW_param_name;
            obj.main_foldername = obj.extractFolderName(path);

            % Initialize Pulsewave-related paths
            obj.initializePaths(OverWrite);

            % Load parameters from cache or fall back to defaults
            obj.loadParameters(path);

            % Calculate Scaling Factors for velocity
            obj.calculateScalingFactors(path);

            % Set up logging (diary)
            obj.setupLogging();

            % Copy input parameters to result folder
            obj.copyInputParameters();

        end

        function mainFolder = extractFolderName(obj, path)
            % Helper function to extract the folder name
            split_path = strsplit(path, filesep);
            mainFolder = split_path{end - 1};
        end

        function initializePaths(obj, OverWrite)
            % Helper function to initialize paths for storing pulsewave-related data

            % Define main and subdirectories for storing data
            obj.PW_path_main = fullfile(obj.PW_path, 'pulsewave');
            PW_folder_name = strcat(obj.main_foldername, '_PW');
            
            % Create or identify a unique folder for the current run
            idx = obj.getUniqueFolderIndex(PW_folder_name, OverWrite);

            % Set the folder name and paths for various data types
            obj.PW_folder_name = sprintf('%s_%d', PW_folder_name, idx);
            obj.PW_path_dir = fullfile(obj.PW_path_main, obj.PW_folder_name);
            obj.PW_path_png = fullfile(obj.PW_path_dir, 'png');
            obj.PW_path_eps = fullfile(obj.PW_path_dir, 'eps');
            obj.PW_path_txt = fullfile(obj.PW_path_dir, 'txt');
            obj.PW_path_avi = fullfile(obj.PW_path_dir, 'avi');
            obj.PW_path_gif = fullfile(obj.PW_path_dir, 'gif');
            obj.PW_path_mp4 = fullfile(obj.PW_path_dir, 'mp4');
            obj.PW_path_json = fullfile(obj.PW_path_dir, 'json');
            obj.PW_path_log = fullfile(obj.PW_path_dir, 'log');

            % Create directories if they don't exist
            obj.createDirectories();
        end

        function idx = getUniqueFolderIndex(obj, folderBaseName, OverWrite)
            % Helper function to determine the unique folder index based on existing directories

            idx = 0;
            list_dir = dir(obj.PW_path_main);
            for i = 1:length(list_dir)
                if contains(list_dir(i).name, folderBaseName)
                    match = regexp(list_dir(i).name, '\d+$', 'match');
                    if ~isempty(match) && str2double(match{1}) >= idx
                        if isempty(OverWrite) || ~OverWrite
                            idx = str2double(match{1}) + 1;  % Use the next index
                        else
                            idx = str2double(match{1});
                        end
                    end
                end
            end
        end

        function createDirectories(obj)
            % Helper function to create necessary directories if they don't exist

            dirs = {obj.PW_path_dir, obj.PW_path_png, obj.PW_path_eps, obj.PW_path_gif, ...
                    obj.PW_path_txt, obj.PW_path_avi, obj.PW_path_mp4, obj.PW_path_json, ...
                    obj.PW_path_log};

            for i = 1:length(dirs)
                if ~isfolder(dirs{i})
                    mkdir(dirs{i});
                end
            end
        end

        function loadParameters(obj, path)
            % Load or fall back to default parameters from cache or config files

            % Try loading parameters from existing .mat or .json files
            if isfile(fullfile(path, 'mat', [obj.main_foldername, '.mat']))
                disp('Reading cache parameters from .mat');
                load(fullfile(path, 'mat', [obj.main_foldername, '.mat']), 'cache');
                obj.stride = cache.batch_stride;
                obj.fs = cache.Fs / 1000;  % Convert Hz to kHz
                obj.f1 = cache.time_transform.f1;
                obj.f2 = cache.time_transform.f2;
                disp('Done.');
            elseif isfile(fullfile(path, 'Holovibes_rendering_parameters.json'))
                disp('Reading cache parameters from Holovibes');
                json_txt = fileread(fullfile(path, 'Holovibes_rendering_parameters.json'));
                footer_parsed = jsondecode(json_txt);
                obj.stride = footer_parsed.compute_settings.image_rendering.time_transformation_stride;
                obj.fs = footer_parsed.info.camera_fps / 1000;  % Convert FPS to kHz
                obj.f1 = footer_parsed.compute_settings.view.z.start / footer_parsed.compute_settings.image_rendering.time_transformation_size * obj.fs;
                obj.f2 = obj.fs / 2;
                disp('Done.');
            else
                % Default values if no parameters are found
                disp('WARNING: No rendering parameters file found. Using default values.');
                obj.stride = 500;
                obj.fs = 34;  % Default value in kHz
                obj.f1 = 6;
                obj.f2 = 15;
            end
        end

        function calculateScalingFactors(obj, path)
            % Calculate scaling factors based on pulsewave parameters

            PW_params = Parameters_json(obj.PW_path, obj.PW_param_name);
            obj.ScalingFactorVelocityInPlane = 1000 * 1000 * 2 * PW_params.lambda / sin(PW_params.phi);  % 6.9 mm/s / kHz
        end

        function setupLogging(obj)
            % Set up logging (diary) for the current session

            diary off  % Turn off logging first to avoid logging this script
            diary_filename = fullfile(obj.PW_path_log, sprintf('%s_log.txt', obj.main_foldername));
            set(0, 'DiaryFile', diary_filename);
            diary on  % Turn on logging
            fprintf("==========================================\n");
            fprintf("Current Folder Path: %s\n", obj.PW_path);
            fprintf("Current File: %s\n", obj.PW_folder_name);
            fprintf("Start Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
            fprintf("==========================================\n");
        end

        function copyInputParameters(obj)
            % Copy the input parameters to the result folder

            path_dir_json = fullfile(obj.PW_path, 'pulsewave', 'json');
            path_file_json_params = fullfile(path_dir_json, obj.PW_param_name);
            copyfile(path_file_json_params, obj.PW_path_json);
        end

        function Params = getParams(obj)
            Params = Parameters_json(obj.PW_path,obj.PW_param_name);
        end

    end
end
