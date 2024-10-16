classdef OneCycleClass

    properties

        M0_data_video % M0 raw
        M1_data_video % M1 raw
        M2_data_video % M2 raw
        M0_disp_video % M0 AVI

        f_RMS_video % RMS M2/M0
        f_AVG_video % AVG M1/M0

        SH_data_hypervideo % SH raw

        directory char
        filenames
        k double %interpolaton parameter
        load_logs char
        % For sectioning
        video_loaded
        pictureSection
        videoSection
        flag_SH_analysis
        flag_PulseWave_analysis
        flag_velocity_analysis
        flag_ARI_analysis
        flag_bloodVolumeRate_analysis
        flag_bloodVelocityProfile_analysis

        %% FIXME : relancer a chaque rendering
        ToolBoxmaster ToolBoxClass

    end

    methods

        function obj = OneCycleClass(path)

            arguments
                path
            end

            disp(path)

            if ~isfolder(path) % if holo file with moments inside is the input
                [filepath, name, ~] = fileparts(path);
                mkdir(fullfile(filepath, name));
                obj.directory = fullfile(filepath, name);
                obj.filenames = name;
            else
                obj.directory = path;
                tmp_idx = regexp(path, '\');
                obj.filenames = path(tmp_idx(end - 1) + 1:end - 1);
            end

            obj.k = 0;

            %% AVEC TXT

            %             try
            %                 checkPulsewaveParamsFromTxt(obj.directory);
            %                 PW_params = Parameters(obj.directory);
            %             catch
            %                 dir_path_txt = fullfile(path,'txt');
            %                 delete(fullfile(dir_path_txt,'InputPulsewaveParams.txt'));
            %                 checkPulsewaveParamsFromTxt(obj.directory);
            %                 PW_params = Parameters(obj.directory);
            %             end

            %% AVEC JSON

            checkPulsewaveParamsFromJson(obj.directory);
            PW_params = Parameters_json(obj.directory);

            obj.k = PW_params.k;
            obj.ToolBoxmaster = ToolBoxClass(obj.directory);

            obj.load_logs = '\n=== LOADING : \r\n';
            logs = obj.load_logs;

            %% Display Video loading

            if ~isfolder(path) % if holo file with moments inside is the input
                disp(['reading moments in : ', strcat(obj.directory, '.holo')]);
                [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
                readMomentsFooter(obj.directory);
                obj.M0_disp_video = ff_correction(videoM0, 30);
                obj.M0_data_video = videoM0;
                obj.M1_data_video = videoM1;
                obj.M2_data_video = videoM2;
            else
                %% Ref loading
                dir_path_avi = fullfile(obj.directory, 'avi');
                NameRefAviFile = strcat(obj.filenames, '_M0');
                RefAviFilePath = fullfile(dir_path_avi, NameRefAviFile);
                ext = '.avi';

                disp(['reading : ', RefAviFilePath]);
                RefAviFilePath(strfind(RefAviFilePath, '\')) = '/';
                str_tosave = sprintf('reading : %s', RefAviFilePath);
                logs = strcat(logs, '\r', str_tosave);

                V = VideoReader(fullfile(dir_path_avi, [NameRefAviFile, ext]));
                M0_disp_video = zeros(V.Height, V.Width, V.NumFrames);

                for n = 1:V.NumFrames
                    M0_disp_video(:, :, n) = rgb2gray(read(V, n));
                end

                obj.M0_disp_video = M0_disp_video;

                clear V M0_disp_video

                %% Data Videos loading

                refvideosize = size(obj.M0_disp_video);
                dir_path_raw = fullfile(obj.directory, 'raw');
                ext = '.raw';

                try
                    % Import Moment 0

                    NameRawFile = strcat(obj.filenames, '_moment0');
                    disp(['reading : ', fullfile(dir_path_raw, [NameRawFile, ext])]);
                    FilePathUnix = fullfile(dir_path_raw, [NameRawFile, ext]);
                    FilePathUnix(strfind(FilePathUnix, '\')) = '/';
                    str_tosave = sprintf('reading : %s', FilePathUnix);
                    logs = strcat(logs, '\r', str_tosave);

                    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
                    M0_data_video = fread(fileID, 'float32');
                    fclose(fileID);
                    obj.M0_data_video = reshape(M0_data_video, refvideosize);

                    % Import Moment 1

                    NameRawFile = strcat(obj.filenames, '_moment1');
                    disp(['reading : ', fullfile(dir_path_raw, [NameRawFile, ext])]);
                    FilePathUnix = fullfile(dir_path_raw, [NameRawFile, ext]);
                    FilePathUnix(strfind(FilePathUnix, '\')) = '/';
                    str_tosave = sprintf('reading : %s', FilePathUnix);
                    logs = strcat(logs, '\r', str_tosave);

                    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
                    M1_data_video = fread(fileID, 'float32');
                    fclose(fileID);
                    obj.M1_data_video = reshape(M1_data_video, refvideosize);

                    % Import Moment 2

                    NameRawFile = strcat(obj.filenames, '_moment2');
                    disp(['reading : ', fullfile(dir_path_raw, [NameRawFile, ext])]);
                    FilePathUnix = fullfile(dir_path_raw, [NameRawFile, ext]);
                    FilePathUnix(strfind(FilePathUnix, '\')) = '/';
                    str_tosave = sprintf('reading : %s', FilePathUnix);
                    logs = strcat(logs, '\r', str_tosave);

                    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
                    M2_data_video = fread(fileID, 'float32');
                    fclose(fileID);
                    obj.M2_data_video = reshape(M2_data_video, refvideosize);

                catch ME
                    disp(['ID: ' ME.identifier])
                    % fprintf("No Raw Files, please select a Holowave folder with raw files exported")
                    rethrow(ME)

                end

                obj.load_logs = logs;

            end

        end

        function obj = MomentNormalize(obj)

            PW_params = Parameters_json(obj.directory);
            M0_data_mean = mean(obj.M0_data_video, [1 2]);

            obj.f_RMS_video = sqrt(obj.M2_data_video ./ M0_data_mean);
            obj.f_AVG_video = obj.M1_data_video ./ M0_data_mean;
            obj.M0_disp_video = flat_field_correction(obj.M0_disp_video, ceil(PW_params.flatField_gwRatio * size(obj.M0_disp_video, 1)), PW_params.flatField_border);

        end

        function obj = registerVideo(obj)
            tic
            % Registers the video using intensity based registration
            PW_params = Parameters_json(obj.directory);

            if ~PW_params.registerVideoFlag
                return % do nothing if not required
            end

            video = obj.M0_disp_video;
            numX = size(video, 1);
            numY = size(video, 2);
            x = linspace(-numX / 2, numX / 2, numX);
            y = linspace(-numY / 2, numY / 2, numY);
            [X, Y] = meshgrid(x, y);

            disc_ratio = 0.7; % parametrize this coef if needed
            disc = X .^ 2 + Y .^ 2 < (disc_ratio * min(numX, numY) / 2) ^ 2;
            video_reg = video .* disc - disc .* sum(video .* disc, [1, 2]) / nnz(disc); % minus the mean in the disc of each frame
            video_reg = reshape(video_reg, size(video, 1), size(video, 2), 1, size(video, 3)); % insert a dimension to match reegistration functions

            video_reg = video_reg ./ (max(abs(video_reg), [], [1, 2])); % rescaling each frame but keeps mean at zero

            image_ref = mean(video_reg(:, :, PW_params.refAvgStart:PW_params.refAvgEnd), 3); % ref image is from 10 to 20
            [~, shifts] = register_video_from_reference(video_reg, image_ref);

            obj.M0_disp_video = register_video_from_shifts(video, shifts);

            obj.M0_data_video = register_video_from_shifts(obj.M0_data_video, shifts);
            obj.M1_data_video = register_video_from_shifts(obj.M1_data_video, shifts);
            obj.M2_data_video = register_video_from_shifts(obj.M2_data_video, shifts);

            toc

        end

        function obj = cropAllVideo(obj)
            %Crop a video (matrix dim 3)
            PW_params = Parameters_json(obj.directory);
            firstFrame = PW_params.videoStartFrameIndex;
            lastFrame = PW_params.videoEndFrameIndex;
            [~, ~, numFrames] = size(obj.M0_disp_video);

            logs = obj.load_logs;

            if firstFrame > 0 && firstFrame < numFrames && lastFrame > firstFrame && lastFrame <= numFrames

                obj.M0_disp_video = obj.M0_disp_video(:, :, firstFrame:lastFrame);
                obj.M0_data_video = obj.M0_data_video(:, :, firstFrame:lastFrame);
                obj.M1_data_video = obj.M1_data_video(:, :, firstFrame:lastFrame);
                obj.M2_data_video = obj.M2_data_video(:, :, firstFrame:lastFrame);

                disp(['Data cube frame: ', num2str(firstFrame), '/', num2str(numFrames), ' to ', num2str(lastFrame), '/', num2str(numFrames)])

                str_tosave = sprintf('Data cube frame: %s/%s to %s/%s', num2str(firstFrame), num2str(numFrames), num2str(lastFrame), num2str(numFrames));
                logs = strcat(logs, '\r', str_tosave);
            else
                disp('Wrong value for the first frame. Set as 1.')
                disp('Wrong value for the last frame. Set as the end.')
                disp(['Data cube frame: 1/', num2str(numFrames), ' to ', num2str(numFrames), '/', num2str(numFrames)])

                str_tosave = sprintf('Wrong value for the first frame. Set as 1. \rWrong value for the last frame. Set as the end. \rData cube frame: 1/%s to %s/%s', num2str(numFrames), num2str(numFrames), num2str(numFrames));
                logs = strcat(logs, '\r\n\n', str_tosave, '\n');
            end

            obj.load_logs = logs;
        end

        function obj = VideoResize(obj)
            tic
            PW_params = Parameters_json(obj.directory);
            out_height = PW_params.frameHeight;
            out_width = PW_params.frameHeight;
            out_numFrames = PW_params.videoLength;

            if out_numFrames < 0
                return % do nothing if not required
            end

            if out_numFrames < numFrames
                % we average the input images to get out_numFrames frames

                % batch = floor(numFrames / out_numFrames);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.M0_disp_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.M0_disp_video = single(tmp_ref);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.M0_data_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.M0_data_video = single(tmp_ref);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.M1_data_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.M1_data_video = single(tmp_ref);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.M2_data_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.M2_data_video = single(tmp_ref);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.f_AVG_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.f_AVG_video = single(tmp_ref);

                tmp_ref = zeros([numX, numY, out_numFrames]);
                tmp_calc = obj.f_RMS_video;

                for i = 1:out_numFrames
                    tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
                end

                obj.f_RMS_video = single(tmp_ref);

            end

            if out_height < 0 && out_width < 0
                return % do nothing if not required
            end

            if out_height < 0
                out_height = numX;
            end

            if out_width < 0
                out_width = numY;
            end

            if out_numFrames < 0
                out_numFrames = numFrames;
            end

            [Xq, Yq, Zq] = meshgrid(linspace(1, numY, out_width), linspace(1, numX, out_height), linspace(1, numFrames, out_numFrames));
            % tmp_ref = zeros(numX, numY, numFrames);

            tmp_calc_ref = obj.M0_disp_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.M0_disp_video = single(tmp_ref);

            tmp_calc_ref = obj.M0_data_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.M0_data_video = single(tmp_ref);

            tmp_calc_ref = obj.M1_data_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.M1_data_video = single(tmp_ref);

            tmp_calc_ref = obj.M2_data_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.M2_data_video = single(tmp_ref);

            tmp_calc_ref = obj.f_AVG_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.f_AVG_video = single(tmp_ref);

            tmp_calc_ref = obj.f_RMS_video;
            tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
            obj.f_RMS_video = single(tmp_ref);

            disp(['Resized data cube : ', num2str(out_width), 'x', num2str(out_height), 'x', num2str(out_numFrames)])
            % logs = obj.load_logs;
            % str_tosave = sprintf("Resized data cube : %s x %s x %s", num2str(out_width), num2str(out_height), num2str(out_numFrames));
            % logs = strcat(logs, '\r\n\n', str_tosave, '\n');
            toc
        end

        function obj = RemoveOutliers(obj)
            %% Outlier Cleaning
            [numX, numY, numFrames] = size(obj.f_RMS_video);
            window_size = ceil(numFrames / 50);

            tmp_f_RMS_cleaned = zeros(numX, numY, numFrames);
            tmp_f_AVG_cleaned = zeros(numX, numY, numFrames);
            tmp_f_RMS = obj.f_RMS_video;
            tmp_f_AVG = obj.f_AVG_video;

            parfor xx = 1:numX

                for yy = 1:numY
                    tmp_f_RMS_cleaned(xx, yy, :) = filloutliers(tmp_f_RMS(xx, yy, :), 'linear', 'movmedian', window_size);
                end

            end

            obj.f_RMS_video = tmp_f_RMS_cleaned;

            clear tmp_f_RMS_cleaned tmp_f_RMS

            parfor xx = 1:numX

                for yy = 1:numY
                    tmp_f_AVG_cleaned(xx, yy, :) = filloutliers(tmp_f_AVG(xx, yy, :), 'linear', 'movmedian', window_size);
                end

            end

            obj.f_AVG_video = tmp_f_AVG_cleaned;

            clear tmp_f_AVG_cleaned tmp_f_AVG
        end

        function obj = Interpolate(obj) %ref = TRUE indicates the object is the reference
            [numX, numY, numFrames] = size(obj.M0_disp_video);
            kInterp = obj.k;
            numX = (numX - 1) * (2 ^ kInterp - 1) + numX;
            numY = (numY - 1) * (2 ^ kInterp - 1) + numY;

            % Reference M0
            tmpReferenceM0 = zeros(numX, numY, numFrames);
            tmpCalcRef = obj.M0_disp_video;

            parfor frameIdx = 1:numFrames
                tmpReferenceM0(:, :, frameIdx) = interp2(tmpCalcRef(:, :, frameIdx), kInterp);
            end

            obj.M0_disp_video = tmpReferenceM0;
            clear tmpReferenceM0 tmpCalcRef

            % M0
            tmpM0 = zeros(numX, numY, numFrames);
            tmpCalcM0 = obj.M0_data_video;

            parfor frameIdx = 1:numFrames % loop over frames
                tmpM0(:, :, frameIdx) = interp2(tmpCalcM0(:, :, frameIdx), kInterp);
            end

            obj.M0_data_video = tmpM0;
            clear tmpM0 tmpCalcM0

            % M1
            tmpM1 = zeros(numX, numY, numFrames);
            tmpCalcM1 = obj.M1_data_video;

            parfor frameIdx = 1:numFrames % loop over frames
                tmpM1(:, :, frameIdx) = interp2(tmpCalcM1(:, :, frameIdx), kInterp);
            end

            obj.M1_data_video = tmpM1;
            clear tmpM1 tmpCalcM1

            % M2
            tmpM2 = zeros(numX, numY, numFrames);
            tmpCalcM2 = obj.M2_data_video;

            parfor frameIdx = 1:numFrames % loop over frames
                tmpM2(:, :, frameIdx) = interp2(tmpCalcM2(:, :, frameIdx), kInterp);
            end

            obj.M2_data_video = tmpM2;
            clear tmpM2 tmpCalcM2

            % M1M0
            tmpM1M0 = zeros(numX, numY, numFrames);
            tmpCalcM1M0 = obj.f_AVG_video;

            parfor frameIdx = 1:numFrames
                tmpM1M0(:, :, frameIdx) = interp2(tmpCalcM1M0(:, :, frameIdx), kInterp);
            end

            obj.f_AVG_video = tmpM1M0;
            clear tmpM1M0 tmpCalcM1M0

            % M2M0
            tmpM2M0 = zeros(numX, numY, numFrames);
            tmpCalcM2M0 = obj.f_RMS_video;

            parfor frameIdx = 1:numFrames
                tmpM2M0(:, :, frameIdx) = interp2(tmpCalcM2M0(:, :, frameIdx), kInterp);
            end

            obj.f_RMS_video = single(tmpM2M0);
            clear tmpM2M0 tmpCalcM2M0

        end

        function onePulse(obj, Ninterp)
            %  ------- This is the app main routine. --------

            % PW_params = Parameters_json(obj.directory);
            % ToolBox = obj.ToolBoxmaster;

            % progress_bar = waitbar(0,'');


            checkPulsewaveParamsFromJson(obj.directory);
            PW_params = Parameters_json(obj.directory);
            totalTime = tic;

            obj.k = PW_params.k;
            obj.ToolBoxmaster = ToolBoxClass(obj.directory);
            ToolBox = obj.ToolBoxmaster;

            fprintf("==================================\n")
            fprintf("File: %s\n", ToolBox.PW_folder_name)
            fprintf("==================================\n")
            fprintf("Loading Input Parameters\n")

            mkdir(ToolBox.PW_path_dir);
            mkdir(ToolBox.PW_path_png);
            mkdir(ToolBox.PW_path_eps);
            mkdir(ToolBox.PW_path_gif);
            mkdir(ToolBox.PW_path_txt);
            mkdir(ToolBox.PW_path_avi);
            mkdir(ToolBox.PW_path_mp4);
            mkdir(ToolBox.PW_path_json);
            mkdir(ToolBox.PW_path_log);

            try
                copyfile(fullfile(obj.directory, 'log', 'RenderingParameters.json'), ToolBox.PW_path_log)
            catch
                warning('no rendering parameters were found')
            end

            %             path_dir_txt = fullfile(obj.directory,'txt');
            %             path_file_txt_params = fullfile(path_dir_txt,'InputPulsewaveParams.txt');
            %             copyfile(path_file_txt_params,ToolBox.PW_path_txt );

            path_dir_json = fullfile(obj.directory, 'pulsewave', 'json');
            path_file_json_params = fullfile(path_dir_json, 'InputPulsewaveParams.json');
            copyfile(path_file_json_params, ToolBox.PW_path_json);

            %saving times
            path_file_txt_exe_times = fullfile(ToolBox.PW_path_log, sprintf('%s_execution_times.txt', ToolBox.PW_folder_name));
            fileID = fopen(path_file_txt_exe_times, 'w');
            fprintf(fileID, 'EXECUTION TIMES : \r\n==================\n\r\n');
            fclose(fileID);

            % SAVING GIT VERSION
            % In the txt file in the folder : "log"

            name_log = strcat(ToolBox.PW_folder_name, '_log.txt');
            path_file_log = fullfile(ToolBox.PW_path_log, name_log);

            gitBranchCommand = 'git symbolic-ref --short HEAD';
            [statusBranch, resultBranch] = system(gitBranchCommand);

            if statusBranch == 0
                resultBranch = strtrim(resultBranch);
                MessBranch = 'Current branch : %s \r';
            else
                MessBranch = 'Error getting current branch name.\r Git command output: %s \r';
            end

            gitHashCommand = 'git rev-parse HEAD';
            [statusHash, resultHash] = system(gitHashCommand);

            if statusHash == 0 %hash command was successful
                resultHash = strtrim(resultHash);
                MessHash = 'Latest Commit Hash : %s \r';
            else
                MessHash = 'Error getting latest commit hash.\r Git command output: %s \r';
            end

            fileID = fopen(path_file_log, 'w');

            fprintf(fileID, '==================\rGIT VERSION :\r');
            fprintf(fileID, MessBranch, resultBranch);
            fprintf(fileID, MessHash, resultHash);
            fprintf(fileID, '==================\r\n ');

            fprintf(fileID, obj.load_logs);

            fprintf(fileID, '\r\n=== EXECUTION \r\n\n');

            fclose(fileID);

            %% Creating Masks
            % waitbar(0.1,progress_bar,"Creating Masks");

            fileID = fopen(path_file_txt_exe_times, 'a+');
            fprintf(fileID, '\r\n');
            fclose(fileID);

            createMasksTiming = tic;

            fprintf("\n----------------------------------\n")
            fprintf("Mask Creation\n")
            fprintf("----------------------------------\n")

            f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));

            try
                [maskArtery, maskVein, ~, maskBackground, ~, ~, ~] = forceCreateMasks(obj.M0_disp_video, f_AVG_mean, obj.directory, ToolBox);
            catch
                [maskArtery, maskVein, ~, maskBackground, ~, ~, ~] = createMasks(obj.M0_disp_video, obj.f_RMS_video, f_AVG_mean, obj.directory, ToolBox);
            end

            time_create_masks = toc(createMasksTiming);
            fprintf("- Mask Creation took : %ds\n", round(time_create_masks))
            save_time(path_file_txt_exe_times, 'CreateMasks', time_create_masks)

            fileID = fopen(path_file_txt_exe_times, 'a+');
            fprintf(fileID, '\r\n----------\r\n');
            fclose(fileID);

            %% PulseWave Analysis
            % waitbar(0.25,progress_bar,"PulseWave analysis");

            if obj.flag_PulseWave_analysis
                close all

                fileID = fopen(path_file_txt_exe_times, 'a+');
                fprintf(fileID, 'PULSEWAVE ANALYSIS : \r\n\n');
                fclose(fileID);

                findSystoleTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Find Systole\n")
                fprintf("----------------------------------\n")

                [sysIdxList, ~] = find_systole_index(obj.M0_disp_video, obj.directory, maskArtery);

                time_sys_idx = toc(findSystoleTimer);
                fprintf("- FindSystoleIndex took : %ds\n", round(time_sys_idx))
                save_time(path_file_txt_exe_times, 'Find Systole Index', time_sys_idx)

                pulseAnalysisTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Pulse Analysis\n")
                fprintf("----------------------------------\n")

                [vOneCycle, vRMS, exec_times] = pulseAnalysis(Ninterp, obj.f_RMS_video, f_AVG_mean, obj.M2_data_video, obj.M0_data_video, sysIdxList, maskArtery, maskVein, maskBackground, ToolBox, obj.directory);

                time_pulseanalysis = toc(pulseAnalysisTimer);
                fprintf("- Pulse Analysis took : %ds\n", round(time_pulseanalysis))
                save_time(path_file_txt_exe_times, 'Pulse Analysis', time_pulseanalysis)

                %exec time details
                fileID = fopen(path_file_txt_exe_times, 'a+');

                for i = 1:size(exec_times, 2)
                    fprintf(fileID, '\t%s : %.0fs \r\n', exec_times(1, i), exec_times(2, i));
                end

                fclose(fileID);
                clear exec_times

                % pulseVelocityTimer = tic;
                %
                % fprintf("\n----------------------------------\n")
                % fprintf("Pulse Velocity\n")
                % fprintf("----------------------------------\n")
                %
                % pulseVelocity(obj.M0_data_video, maskArtery, ToolBox, obj.directory)
                %
                % time_pulsevelocity = toc(pulseVelocityTimer);
                % fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
                % save_time(path_file_txt_exe_times, 'Pulse Velocity', time_pulsevelocity)

                if obj.flag_velocity_analysis
                    bloodFlowVelocityTimer = tic;

                    fprintf("\n----------------------------------\n")
                    fprintf("Blood Flow Velocity Calculation\n")
                    fprintf("----------------------------------\n")

                    bloodFlowVelocity(vRMS, vOneCycle, maskArtery, maskVein, obj.M0_disp_video, ToolBox, obj.directory)
                    % bloodFlowVelocityFullField(vRMS, vOneCycle, maskArtery, maskVein, obj.M0_data_video, ToolBox, obj.directory)

                    time_velo = toc(bloodFlowVelocityTimer);
                    fprintf("- Blood Flow Velocity calculation took : %ds\n", round(time_velo))
                    save_time(path_file_txt_exe_times, 'Blood Flow Velocity', time_velo)
                end

                if obj.flag_ARI_analysis
                    bloodARITimer = tic;

                    fprintf("\n----------------------------------\n")
                    fprintf("ARI and API Calculation\n")
                    fprintf("----------------------------------\n")

                    ArterialResistivityIndex(vOneCycle, obj.M0_disp_video, maskArtery, ToolBox, obj.directory);

                    time_ARIAPI = toc(bloodARITimer);
                    fprintf("- ARI & API calculation took : %ds\n", round(time_ARIAPI))
                    save_time(path_file_txt_exe_times, 'Arterial Resistivity Index', time_ARIAPI)
                end

                if obj.flag_bloodVolumeRate_analysis
                    bloodVolumeRateTimer = tic;

                    fprintf("\n----------------------------------\n")
                    fprintf("Blood Volume Rate Calculation\n")
                    fprintf("----------------------------------\n")

                    bloodVolumeRate(maskArtery, maskVein, vRMS, obj.M0_disp_video, ToolBox, obj.k, obj.directory, obj.flag_bloodVelocityProfile_analysis);
                    bloodVolumeRateForAllRadii(maskArtery, maskVein, vRMS, obj.M0_disp_video, ToolBox, obj.k, obj.directory, obj.flag_bloodVelocityProfile_analysis);

                    time_volumeRate = toc(bloodVolumeRateTimer);
                    fprintf("- Blood Volume rate calculation took : %ds\n", round(time_volumeRate))
                    save_time(path_file_txt_exe_times, 'Blood Volume rate', time_volumeRate)
                end

                tTotal = toc(totalTime);

                fprintf("\n----------------------------------\n")
                fprintf("Total Pulsewave timing : %ds\n" , round(tTotal))
                fileID = fopen(path_file_txt_exe_times, 'a+');
                fprintf(fileID, '\r\n=== Total : %.0fs \r\n\n----------\r\n', tTotal);
                fclose(fileID);
            end

            %% Spectrum Analysis
            % waitbar(0.9,progress_bar,"Spectrum analysis");

            if obj.flag_SH_analysis

                %% Import SH

                tmpname = strcat(ToolBox.main_foldername, '_SH');
                ext = '.raw';
                disp(['reading : ', fullfile(obj.directory, 'raw', [tmpname, ext])]);
                fileID = fopen(fullfile(obj.directory, 'raw', [tmpname, ext]));
                videoSH = fread(fileID, 'float32');
                fclose(fileID);
                SH_cube = reshape(videoSH, numX, numY, numFrames, []);

                tic
                spectrum_analysis(maskArtery, maskBackground, SH_cube, ToolBox, obj.M0_data_video);
                disp('Spectrum analysis :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrum_analysis', time)

                tic
                spectrogram(maskArtery, maskBackground, SH_cube, ToolBox);
                disp('Spectrogram timing :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrogram', time)
            end

            displaySuccessMsg(1);
            % close(progress_bar)
            close all

        end

    end

end
