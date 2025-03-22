classdef ToolBoxClass < handle
% ToolBoxClass holds useful variables for eyeflow processing.

properties
    % Paths
    EF_path char
    EF_name char
    path_main char
    path_dir char
    path_png char
    path_eps char
    path_pdf char
    path_gif char
    path_txt char
    path_avi char
    path_mp4 char
    path_json char
    path_log char
    main_foldername char
    param_name char
    folder_name char
    % Parameters
    stride double
    fs double
    f1 double
    f2 double
end

methods

    function obj = ToolBoxClass(path, EF_param_name, OverWrite)
        % Constructor for ToolBoxClass: Initializes paths, parameters, and calculates scaling factors.

        % Store paths and parameters
        obj.EF_path = path;
        obj.param_name = EF_param_name;
        obj.main_foldername = obj.extractFolderName(path);

        % Initialize EyeFlow-related paths
        obj.initializePaths(OverWrite);

        % Load parameters from cache or fall back to defaults
        obj.loadParameters(path);

        % Set up logging (diary)
        obj.setupLogging();

        % Copy input parameters to result folder
        obj.copyInputParameters();

        obj.setGlobalToolBox;

    end

    function mainFolder = extractFolderName(~, path)
        % Helper function to extract the folder name
        split_path = strsplit(path, filesep);
        mainFolder = split_path{end - 1};
    end

    function setGlobalToolBox(obj)
        global ToolBoxGlobal
        ToolBoxGlobal = obj;
    end

    function initializePaths(obj, OverWrite)
        % Helper function to initialize paths for storing eyeflow-related data

        % Define main and subdirectories for storing data
        obj.path_main = fullfile(obj.EF_path, 'eyeflow');
        foldername_EF = strcat(obj.main_foldername, '_EF');

        % Create or identify a unique folder for the current run
        idx = obj.getUniqueFolderIndex(foldername_EF, OverWrite);

        % Set the folder name and paths for various data types
        obj.EF_name = foldername_EF;
        obj.folder_name = sprintf('%s_%d', foldername_EF, idx);
        obj.path_dir = fullfile(obj.path_main, obj.folder_name);
        obj.path_png = fullfile(obj.path_dir, 'png');
        obj.path_eps = fullfile(obj.path_dir, 'eps');
        obj.path_pdf = fullfile(obj.path_dir, 'pdf');
        obj.path_txt = fullfile(obj.path_dir, 'txt');
        obj.path_avi = fullfile(obj.path_dir, 'avi');
        obj.path_gif = fullfile(obj.path_dir, 'gif');
        obj.path_mp4 = fullfile(obj.path_dir, 'mp4');
        obj.path_json = fullfile(obj.path_dir, 'json');
        obj.path_log = fullfile(obj.path_dir, 'log');

        % Create directories if they don't exist
        obj.createDirectories();
    end

    function idx = getUniqueFolderIndex(obj, folderBaseName, OverWrite)
        % Helper function to determine the unique folder index based on existing directories

        idx = 0;
        list_dir = dir(obj.path_main);

        for i = 1:length(list_dir)

            if contains(list_dir(i).name, folderBaseName)
                match = regexp(list_dir(i).name, '\d+$', 'match');

                if ~isempty(match) && str2double(match{1}) >= idx

                    if isempty(OverWrite) || ~OverWrite
                        idx = str2double(match{1}) + 1; % Use the next index
                    else
                        idx = str2double(match{1});
                    end

                end

            end

        end

    end

    function createDirectories(obj)
        % Helper function to create necessary directories if they don't exist

        dirs = {obj.path_dir, obj.path_png, obj.path_eps, obj.path_gif, ...
                    obj.path_txt, obj.path_avi, obj.path_mp4, obj.path_json, ...
                    obj.path_log, obj.path_pdf};

        for i = 1:length(dirs)

            if ~isfolder(dirs{i})
                mkdir(dirs{i});
            end

        end

    end

    function loadParameters(obj, path)
        % Load or fall back to default parameters from cache or config files

        % Try loading parameters from existing .mat or .json files
        if ~isempty(dir(fullfile(path, ['*', 'RenderingParameters', '*']))) % since HD 2.0
            disp('Reading cache parameters from .json');
            fpath = fullfile(path, dir(fullfile(path, ['*', 'RenderingParameters', '*'])).name);
            decoded_data = jsondecode(fileread(fpath));
            obj.stride = decoded_data.batch_stride;
            obj.fs = decoded_data.fs; % Convert kHz to kHz
            obj.f1 = decoded_data.time_range(1);
            obj.f2 = decoded_data.time_range(2);
            disp('Done.');
        elseif isfile(fullfile(path, 'mat', [obj.main_foldername, '.mat']))
            disp('Reading cache parameters from .mat');
            load(fullfile(path, 'mat', [obj.main_foldername, '.mat']), 'cache');
            obj.stride = cache.batch_stride;
            obj.fs = cache.Fs / 1000; % Convert Hz to kHz
            obj.f1 = cache.time_transform.f1;
            obj.f2 = cache.time_transform.f2;
        elseif isfile(fullfile(path, 'Holovibes_rendering_parameters.json'))
            json_txt = fileread(fullfile(path, 'Holovibes_rendering_parameters.json'));
            footer_parsed = jsondecode(json_txt);
            obj.stride = footer_parsed.compute_settings.image_rendering.time_transformation_stride;
            obj.fs = footer_parsed.info.camera_fps / 1000; % Convert FPS to kHz
            obj.f1 = footer_parsed.compute_settings.view.z.start / footer_parsed.compute_settings.image_rendering.time_transformation_size * obj.fs;
            obj.f2 = obj.fs / 2;
        else
            % Default values if no parameters are found
            disp('WARNING: No rendering parameters file found. Using default values.');
            obj.stride = 500;
            obj.fs = 34; % Default value in kHz
            obj.f1 = 6;
            obj.f2 = 15;
        end

    end

    function setupLogging(obj)
        % Set up logging (diary) for the current session

        diary off % Turn off logging first to avoid logging this script
        diary_filename = fullfile(obj.path_log, sprintf('%s_log.txt', obj.main_foldername));
        set(0, 'DiaryFile', diary_filename);
        diary on % Turn on logging
        fprintf("==========================================\n");
        fprintf("Current Folder Path: %s\n", obj.EF_path);
        fprintf("Current File: %s\n", obj.folder_name);
        fprintf("Start Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
        fprintf("==========================================\n");
    end

    function copyInputParameters(obj)
        % Copy the input parameters to the result folder

        path_dir_json = fullfile(obj.EF_path, 'eyeflow', 'json');
        path_file_json_params = fullfile(path_dir_json, obj.param_name);
        copyfile(path_file_json_params, obj.path_json);
    end

    function Params = getParams(obj)
        Params = Parameters_json(obj.EF_path, obj.param_name);
    end

end

end
