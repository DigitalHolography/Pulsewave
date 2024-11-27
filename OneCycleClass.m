classdef OneCycleClass

    properties

        M0_data_video % M0 raw
        M1_data_video % M1 raw
        M2_data_video % M2 raw
        M0_disp_video % M0 AVI

        f_RMS_video % RMS M2/M0
        f_AVG_video % AVG M1/M0

        SH_data_hypervideo % SH raw
        
        maskArtery
        maskVein
        maskBackground
        maskSection
        maskNeighbors

        sysIdxList % list of frame indexes counting carciac cycles
        baryPos % x y position of the ONH
        vRMS % video estimate of velocity map in retinal vessels

        directory char % directory of input data (from HoloDoppler or HoloVibes)
        filenames char % name used for storing analysis data
        k double %interpolaton parameter % TODO move to pw params
        load_logs char
        
        flag_Segmentation
        flag_SH_analysis
        flag_PulseWave_analysis
        flag_velocity_analysis
        flag_ExtendedPulseWave_analysis
        flag_bloodVolumeRate_analysis
        flag_bloodVelocityProfile_analysis

        ToolBoxmaster ToolBoxClass

    end

    methods

        function obj = OneCycleClass(path)

            disp(path)

            if ~isfolder(path) % if th input path is a .holo file with moments inside 
                [filepath, name, ~] = fileparts(path);
                mkdir(fullfile(filepath, name));
                obj.directory = fullfile(filepath, name);
                obj.filenames = name;
            else
                obj.directory = path;
                tmp_idx = regexp(path, '\');
                obj.filenames = obj.directory(tmp_idx(end-1) + 1:end -1);
            end

            %% Parameters

            checkPulsewaveParamsFromJson(obj.directory); % checks compatibility between found PW params and Default PW params of this version of PW.
            PW_params = Parameters_json(obj.directory);

            obj.k = PW_params.k;
            obj.ToolBoxmaster = ToolBoxClass(obj.directory);

            obj.load_logs = '\n=== LOADING : \r\n';
            

            %% Video loading

            if ~isfolder(path) % if holo file with moments inside is the input
                disp(['reading moments in : ', strcat(obj.directory, '.holo')]);
                [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
                readMomentsFooter(obj.directory);
                obj.M0_disp_video = rescale(ff_correction(videoM0, 30))*255;
                obj.M0_data_video = videoM0;
                obj.M1_data_video = videoM1;
                obj.M2_data_video = videoM2;
            else
                obj = readRaw(obj);
            end

        end

        function obj = preprocessData(obj)
            % register
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Video Registering\n")
            fprintf("----------------------------------\n")
            obj = VideoRegistering(obj);
            fprintf("- Video Registering took : %ds\n", round(toc))

            % crop videos
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Video Cropping\n")
            fprintf("----------------------------------\n")
            obj = VideoCropping(obj);
            fprintf("- Video Cropping took : %ds\n", round(toc))

            % moment normalize
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Moment Normalizing\n")
            fprintf("----------------------------------\n")
            obj = VideoNormalizingLocally(obj);
            fprintf("- Moment Normalizing took : %ds\n", round(toc))

            % Video resize (preprocess interpolation interpolate)
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Video Resizing\n")
            fprintf("----------------------------------\n")
            obj = VideoResizing(obj);
            fprintf("- Video Resizing : %ds\n", round(toc))

            % interpolate
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Video Interpolation\n")
            fprintf("----------------------------------\n")
            obj = VideoInterpolating(obj);
            fprintf("- Video Interpolation : %ds\n", round(toc))

            % remove outliers
            tic
            fprintf("\n----------------------------------\n")
            fprintf("Video Outlier Cleaning\n")
            fprintf("----------------------------------\n")
            obj = VideoRemoveOutliers(obj);
            fprintf("- Video Outlier Cleaning took : %ds\n", round(toc))
        
        end

        function obj = onePulse(obj)
            %  ------- This is the app main routine. --------

            checkPulsewaveParamsFromJson(obj.directory); % checks compatibility between found PW params and Default PW params of this version of PW.
            PW_params = Parameters_json(obj.directory);
            totalTime = tic;

            obj.k = PW_params.k;
            obj.ToolBoxmaster = ToolBoxClass(obj.directory);
            ToolBox = obj.ToolBoxmaster;

            mkdir(ToolBox.PW_path_dir);
            mkdir(ToolBox.PW_path_png);
            mkdir(ToolBox.PW_path_eps);
            mkdir(ToolBox.PW_path_gif);
            mkdir(ToolBox.PW_path_txt);
            mkdir(ToolBox.PW_path_avi);
            mkdir(ToolBox.PW_path_mp4);
            mkdir(ToolBox.PW_path_json);
            mkdir(ToolBox.PW_path_log);

            fprintf("==================================\n")
            fprintf("File: %s\n", ToolBox.PW_folder_name)
            fprintf("==================================\n")
            fprintf("Loading Input Parameters\n")

            try
                copyfile(fullfile(obj.directory, 'log', 'RenderingParameters.json'), ToolBox.PW_path_log)
            catch
                warning('no rendering parameters were found')
            end

            % copyin the input parameters to the result folder
            path_dir_json = fullfile(obj.directory, 'pulsewave', 'json');
            path_file_json_params = fullfile(path_dir_json, 'InputPulsewaveParams.json');
            copyfile(path_file_json_params, ToolBox.PW_path_json);

            % saving times
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
                
                vers = readlines('version.txt');
                MessBranch = ['PulseWave GitHub version ', char(vers)];
            end

            gitHashCommand = 'git rev-parse HEAD';
            [statusHash, resultHash] = system(gitHashCommand);

            if statusHash == 0 %hash command was successful
                resultHash = strtrim(resultHash);
                MessHash = 'Latest Commit Hash : %s \r';
            else
                MessHash = '';
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

            if obj.flag_Segmentation

                fileID = fopen(path_file_txt_exe_times, 'a+');
                fprintf(fileID, '\r\n');
                fclose(fileID);
    
                createMasksTiming = tic;
    
                fprintf("\n----------------------------------\n")
                fprintf("Mask Creation\n")
                fprintf("----------------------------------\n")
    
                [obj.maskArtery, obj.maskVein, ~, obj.maskBackground, ~, ~, obj.maskSection, obj.maskNeighbors] = createMasks(obj.M0_disp_video, obj.f_AVG_video, obj.directory, ToolBox);
    
                time_create_masks = toc(createMasksTiming);
                fprintf("- Mask Creation took : %ds\n", round(time_create_masks))
                save_time(path_file_txt_exe_times, 'CreateMasks', time_create_masks)
    
                fileID = fopen(path_file_txt_exe_times, 'a+');
                fprintf(fileID, '\r\n----------\r\n');
                fclose(fileID);

            end

            %% PulseWave Analysis

            if obj.flag_PulseWave_analysis
                close all

                fileID = fopen(path_file_txt_exe_times, 'a+');
                fprintf(fileID, 'PULSEWAVE ANALYSIS : \r\n\n');
                fclose(fileID);

                findSystoleTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Find Systole\n")
                fprintf("----------------------------------\n")

                [obj.sysIdxList, ~] = find_systole_index(obj.M0_disp_video, obj.directory, obj.maskArtery);

                time_sys_idx = toc(findSystoleTimer);
                fprintf("- FindSystoleIndex took : %ds\n", round(time_sys_idx))
                save_time(path_file_txt_exe_times, 'Find Systole Index', time_sys_idx)

                pulseAnalysisTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Pulse Analysis\n")
                fprintf("----------------------------------\n")
                
                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));

                [obj.vRMS, exec_times] = pulseAnalysis(obj.f_RMS_video, f_AVG_mean, obj.M2_data_video, obj.M0_data_video, obj.M0_disp_video, obj.sysIdxList, obj.maskArtery, obj.maskVein, obj.maskBackground, obj.flag_ExtendedPulseWave_analysis, ToolBox, obj.directory);

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
            end

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

                bloodFlowVelocity(obj.vRMS, obj.maskArtery, obj.maskVein, obj.maskSection, obj.M0_disp_video, ToolBox, obj.directory)
                % bloodFlowVelocityFullField(vRMS, vOneCycle, maskArtery, maskVein, obj.M0_data_video, ToolBox, obj.directory)

                time_velo = toc(bloodFlowVelocityTimer);
                fprintf("- Blood Flow Velocity calculation took : %ds\n", round(time_velo))
                save_time(path_file_txt_exe_times, 'Blood Flow Velocity', time_velo)
            end

            if obj.flag_bloodVolumeRate_analysis
                bloodVolumeRateTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Blood Volume Rate Calculation\n")
                fprintf("----------------------------------\n")

                bloodVolumeRate(obj.maskArtery, obj.maskVein, obj.vRMS, obj.M0_disp_video, ToolBox, obj.k, obj.directory, obj.flag_bloodVelocityProfile_analysis);
                bloodVolumeRateForAllRadii(obj.maskArtery, obj.maskVein, obj.vRMS, obj.M0_disp_video, ToolBox, obj.k, obj.directory, obj.flag_bloodVelocityProfile_analysis,obj.sysIdxList);

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
            

            %% Spectrum Analysis

            if obj.flag_SH_analysis && isfile(fullfile(obj.directory, 'raw', [strcat(ToolBox.main_foldername, '_SH'),  '.raw']))

                %% Import SH


                tmpname = strcat(ToolBox.main_foldername, '_SH');
                ext = '.raw';
                disp(['reading : ', fullfile(obj.directory, 'raw', [tmpname, ext])]);
                fileID = fopen(fullfile(obj.directory, 'raw', [tmpname, ext]));
                videoSH = fread(fileID, 'float32');
                fclose(fileID);
                [numX, numY, numFrames] = size(obj.f_RMS_video);
                bin_x = 4;
                bin_y = 4;
                bin_w = 16;
                bin_t = 1;
                SH_cube = reshape(videoSH, ceil(numX/(2^obj.k*bin_x)), ceil(numY/(2^obj.k*bin_y)),[], ceil(numFrames/bin_t));

                tic
                spectrum_analysis(obj.maskArtery, obj.maskBackground, SH_cube, ToolBox, obj.M0_data_video);
                disp('Spectrum analysis :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrum_analysis', time)

                tic
                spectrogram(obj.maskArtery, obj.maskNeighbors,obj.maskSection, SH_cube, ToolBox);
                disp('Spectrogram timing :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrogram', time)
            end

            displaySuccessMsg(1);
            close all

        end

    end

end
