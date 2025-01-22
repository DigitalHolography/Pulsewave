classdef ToolBoxClass < handle

    % Holds useful variables calculated ones and used in the rest of the
    % script

    properties
        %Path of the PW dir and the output dir inside
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
        PW_path_pulswave char
        main_foldername char
        PW_param_name char
        PW_folder_name char
        stride double
        fs double
        f1 double
        f2 double
        ScalingFactorVelocityInPlane double

    end

    methods

        function obj = ToolBoxClass(path, PW_param_name, OverWrite)

            obj.PW_path = path;
            obj.PW_param_name = PW_param_name;
            PW_params = Parameters_json(obj.PW_path,obj.PW_param_name);

            %% Creating paths
            idx = 0;
            split_path = strsplit(path, '\');
            obj.main_foldername = split_path{end - 1};
           
            obj.PW_path_main = fullfile(path, 'pulsewave');

            if ~exist(obj.PW_path_main, 'dir')
                mkdir(obj.PW_path_main);
            end

            PW_folder_name = strcat(obj.main_foldername, '_PW');

            list_dir = dir(obj.PW_path_main);
            for i=1:length(list_dir)
                if contains(list_dir(i).name, PW_folder_name)
                    match = regexp(list_dir(i).name, '\d+$', 'match');
                    if ~isempty(match) && str2double(match{1}) >= idx
                        if ~(OverWrite) || isempty(OverWrite)
                            idx = str2double(match{1}) + 1; %suffix
                        else
                            idx = str2double(match{1});
                        end
                    end
                end
            end

            % for naming with the minimum possible suffix
            %while (exist(fullfile(obj.PW_path_main, sprintf('%s_%d', PW_folder_name, idx)), 'dir'))
            %    idx = idx + 1;
            %end

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

            mkdir(obj.PW_path_dir);
            mkdir(obj.PW_path_png);
            mkdir(obj.PW_path_eps);
            mkdir(obj.PW_path_gif);
            mkdir(obj.PW_path_txt);
            mkdir(obj.PW_path_avi);
            mkdir(obj.PW_path_mp4);
            mkdir(obj.PW_path_json);
            mkdir(obj.PW_path_log);

            %% Reading Cache Parameters from .mat
            dir_path_mat = fullfile(path, 'mat');
            file_path_mat = fullfile(dir_path_mat, [obj.main_foldername, '.mat']);
            %[~,filename_mat,~] = fileparts(file_path_mat);

            if isfile(file_path_mat) % .mat with cache from holowaves is present, timeline can be computed
                disp('reading cache parameters');
                load(file_path_mat, 'cache');
                obj.stride = cache.batch_stride;
                obj.fs = (cache.Fs) / 1000;
                obj.f1 = cache.time_transform.f1;
                obj.f2 = cache.time_transform.f2;
                disp('done.')

            elseif isfile(fullfile(path, [obj.main_foldername, '.mat']))
                disp('reading cache parameters');
                load(fullfile(path, [obj.main_foldername, '.mat']), 'cache');
                obj.stride = cache.batch_stride;
                obj.fs = (cache.Fs) / 1000; %conversion in kHz
                obj.f1 = cache.time_transform.f1;
                obj.f2 = cache.time_transform.f2;
                disp('done.')
            elseif isfile(fullfile(path, 'Holovibes_rendering_parameters.json'))
                disp('reading cache parameters from holovibes');
                json_txt = fileread(fullfile(path, 'Holovibes_rendering_parameters.json'));
                footer_parsed = jsondecode(json_txt);
                obj.stride = footer_parsed.compute_settings.image_rendering.time_transformation_stride;
                obj.fs = footer_parsed.info.camera_fps/1000; %conversion in kHz
                time_transform_size = footer_parsed.compute_settings.image_rendering.time_transformation_size;
                obj.f1 = footer_parsed.compute_settings.view.z.start / time_transform_size * obj.fs; % on suppose que les fréquences choisies
                % sont symmétriques (la plage va jusqu'à end - start)
                obj.f2 = obj.fs / 2;
                disp('done.')
            else
                disp('WARNING : no rendering parameters file found');
                obj.stride = 500;
                obj.fs = 34; % default values faked
                obj.f1 = 6;
                obj.f2 = 15;
            end

            if isfile(fullfile(path, 'log', 'RenderingParameters.json')) % copies a simple log version for readability
                copyfile(fullfile(path, 'log', 'RenderingParameters.json'), obj.PW_path_log)
            end
            %% Calculation of the Sacling Factors

            % obj.ScalingFactorVelocityInPlane = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (3/PW_params.theta)^(1/2); % 1000 for kHz -> Hz and 1000 for m -> mm
            % obj.ScalingFactorVelocityInPlane = 30;
            obj.ScalingFactorVelocityInPlane = 1000 * 1000 * 2 * PW_params.lambda / sin(PW_params.phi); % ~~6.9 mm/s / kHz
            % obj.ScalingFactorVelocityInPlane = 1000 * 1000 * 2 * 2.4 * PW_params.lambda / 1.33;
            % obj.ScalingFactorVelocityCRA_AVG = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (PW_params.theta / 2); % 1000 for kHz -> Hz and 1000 for m -> mm
            % obj.ScalingFactorVelocityCRA_RMS = 1000 * 1000 * PW_params.lambda / PW_params.opticalIndex * (1 / (2 + 2 * (PW_params.theta ^ 3) / 3)) ^ (1/2); % 1000 for kHz -> Hz and 1000 for m -> mm

            % Turn On Diary Logging
            diary off
            % first turn off diary, so as not to log this script
            diary_filename = fullfile(obj.PW_path_log, sprintf('%s_log.txt', obj.main_foldername));
            % setup temp variable with filename + timestamp, echo off
            set(0, 'DiaryFile', diary_filename)
            % set the objectproperty DiaryFile of hObject 0 to the temp variable filename
            clear diary_filename
            % clean up temp variable
            diary on
            % turn on diary logging
            fprintf("==========================================\n")
            fprintf("Current Folder Path: %s\n", obj.PW_path)
            fprintf("Current File: %s\n", obj.PW_folder_name)
            fprintf("Start Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'))
            fprintf("==========================================\n")

            fprintf("Loading Input Parameters\n")

            % copying the input parameters to the result folder
            path_dir_json = fullfile(obj.PW_path, 'pulsewave', 'json');
            path_file_json_params = fullfile(path_dir_json, obj.PW_param_name);
            copyfile(path_file_json_params, obj.PW_path_json);

        end

    end

end
