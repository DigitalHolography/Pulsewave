classdef OneCycleClass < handle
    properties
        M0_data_video % M0 raw
        M1_data_video % M1 raw
        M2_data_video % M2 raw
        M0_ff_video % M0 AVI

        is_preprocessed % tells if the data has been preprocessed

        f_RMS_video % RMS M2/M0
        f_AVG_video % AVG M1/M0

        SH_data_hypervideo % SH raw

        maskArtery
        maskVein
        maskSection
        maskNeighbors

        sysIdxList % list of frame indexes counting cardiac cycles
        xy_barycenter % x y position of the ONH
        vRMS % video estimate of velocity map in retinal vessels

        directory char % directory of input data (from HoloDoppler or HoloVibes)
        PW_params_names cell % filenames of all the current input parameters ('InputPulsewaveParams.json' for example by default)
        PW_param_name char % current filename
        filenames char % name id used for storing the measured rendered data

        flag_Segmentation
        flag_SH_analysis
        flag_PulseWave_analysis
        flag_velocity_analysis
        flag_bloodVolumeRate_analysis

        OverWrite logical
    end

     methods
        function obj = OneCycleClass(path)
            % Constructor for OneCycleClass.
            % Input: path - directory or .holo file path.

            if ~isfolder(path) % If the input path is a .holo file
                [filepath, name, ~] = fileparts(path);

                if ~isfolder(fullfile(filepath, name)) % Create a result folder
                    mkdir(fullfile(filepath, name));
                end

                obj.directory = fullfile(filepath, name);
                obj.filenames = name;
            else % If input is a directory
                obj.directory = path;
                tmp_idx = regexp(path, '\');
                obj.filenames = obj.directory(tmp_idx(end - 1) + 1:end - 1);
            end

            % Load parameters
            obj.PW_params_names = checkPulsewaveParamsFromJson(obj.directory);
            obj.PW_param_name = obj.PW_params_names{1}; % Default behavior

            % Load video data
            if ~isfolder(path) % If .holo file
                disp(['Reading moments in: ', strcat(obj.directory, '.holo')]);
                [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
                readMomentsFooter(obj.directory);
                obj.M0_ff_video = pagetranspose(improve_video(ff_correction(videoM0, 35), 0.0005, 2, 0));
                obj.M0_data_video = pagetranspose(videoM0);
                obj.M1_data_video = pagetranspose(videoM1 / 1e3); % Rescale M1
                obj.M2_data_video = pagetranspose(videoM2 / 1e6); % Rescale M2
            else
                obj = readRaw(obj);
            end

            obj.is_preprocessed = false;
        end


        function obj = preprocessData(obj)
            % Preprocess video data.

            fprintf("\n----------------------------------\nVideo PreProcessing\n----------------------------------\n");

            % Register video
            tic;
            fprintf("\n----------------------------------\nVideo Registering\n----------------------------------\n");
            obj = VideoRegistering(obj);
            fprintf("- Video Registering took: %ds\n", round(toc));

            % Crop video
            tic;
            fprintf("\n----------------------------------\nVideo Cropping\n----------------------------------\n");
            obj = VideoCropping(obj);
            fprintf("- Video Cropping took: %ds\n", round(toc));

            % Normalize moments
            tic;
            fprintf("\n----------------------------------\nMoment Normalizing\n----------------------------------\n");
            obj = VideoNormalizingLocally(obj);
            fprintf("- Moment Normalizing took: %ds\n", round(toc));

            % Resize video
            tic;
            fprintf("\n----------------------------------\nVideo Resizing\n----------------------------------\n");
            obj = VideoResizing(obj);
            fprintf("- Video Resizing took: %ds\n", round(toc));

            % Interpolate video
            tic;
            fprintf("\n----------------------------------\nVideo Interpolation\n----------------------------------\n");
            obj = VideoInterpolating(obj);
            fprintf("- Video Interpolation took: %ds\n", round(toc));

            % Remove outliers
            tic;
            fprintf("\n----------------------------------\nVideo Outlier Cleaning\n----------------------------------\n");
            obj = VideoRemoveOutliers(obj);
            fprintf("- Video Outlier Cleaning took: %ds\n", round(toc));

            obj.is_preprocessed = true;
        end


        function obj = onePulse(obj)
            % Main routine for PulseWave analysis.

            % Initialize ToolBox and parameters
            ToolBox = ToolBoxClass(obj.directory, obj.PW_param_name, obj.OverWrite);
            PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
            totalTime = tic;
            saveGit;

            %% Mask Creation
            if obj.flag_Segmentation
                fprintf("\n----------------------------------\nMask Creation\n----------------------------------\n");
                createMasksTimer = tic;

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.maskArtery, obj.maskVein, obj.maskSection, obj.maskNeighbors, obj.xy_barycenter] = ...
                    createMasks(obj.M0_ff_video, f_AVG_mean);

                fprintf("- Mask Creation took: %ds\n", round(toc(createMasksTimer)));
            end

            %% PulseWave Analysis
            if obj.flag_PulseWave_analysis
                fprintf("\n----------------------------------\nFind Systole\n----------------------------------\n");
                findSystoleTimer = tic;

                [obj.sysIdxList, ~, sysMaxList, sysMinList] = find_systole_index(obj.M0_ff_video, obj.maskArtery);

                % Log systole results
                fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_PW_main_outputs.txt')), 'a');
                fprintf(fileID, 'Heart beat: %f (bpm) \r\n', 60 / mean(diff(obj.sysIdxList) * ToolBox.stride / ToolBox.fs / 1000));
                fprintf(fileID, 'Systole Indices: %s \r\n', strcat('[', sprintf("%d,", obj.sysIdxList), ']'));
                fprintf(fileID, 'Number of Cycles: %d \r\n', numel(obj.sysIdxList) - 1);
                fprintf(fileID, 'Max Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMaxList), ']'));
                fprintf(fileID, 'Min Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMinList), ']'));
                fprintf(fileID, 'Time diastolic min to systolic max derivative (ms): %f \r\n', ...
                    1000 * mean((obj.sysIdxList(2:end) - sysMinList) * ToolBox.stride / ToolBox.fs / 1000));
                fprintf(fileID, 'Time diastolic min to systolic max (ms): %f \r\n', ...
                    1000 * mean((sysMaxList(2:end) - sysMinList(1:end - 1)) * ToolBox.stride / ToolBox.fs / 1000));
                fclose(fileID);

                fprintf("- FindSystoleIndex took: %ds\n", round(toc(findSystoleTimer)));

                fprintf("\n----------------------------------\nPulse Analysis\n----------------------------------\n");
                pulseAnalysisTimer = tic;

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                obj.vRMS = pulseAnalysis(obj.f_RMS_video, obj.maskArtery, obj.maskVein, obj.maskSection, obj.maskNeighbors);

                if PW_params.params.PulseAnalysis.ExtendedFlag
                    extendedPulseAnalysis(obj.M0_ff_video, obj.f_RMS_video, f_AVG_mean, obj.vRMS, obj.maskArtery, obj.maskVein, obj.maskSection, obj.sysIdxList);
                end

                fprintf("- Pulse Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
            end

            %% Pulse Velocity Analysis
            %  if obj.flag_pulseVelocity_analysis
            %     fprintf("\n----------------------------------\nPulse Velocity Calculation\n----------------------------------\n");
            %     pulseVelocityTimer = tic;

            %     pulseVelocity(obj.M0_data_video, maskArtery)

            %     time_pulsevelocity = toc(pulseVelocityTimer);
            %     fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
            % end

            %% Blood Flow Velocity Analysis
            if obj.flag_velocity_analysis
                fprintf("\n----------------------------------\nBlood Flow Velocity Calculation\n----------------------------------\n");
                bloodFlowVelocityTimer = tic;

                bloodFlowVelocity(obj.vRMS, obj.maskArtery, obj.maskVein, obj.maskSection, obj.M0_ff_video, obj.xy_barycenter);

                fprintf("- Blood Flow Velocity calculation took: %ds\n", round(toc(bloodFlowVelocityTimer)));
            end

            %% Blood Volume Rate Analysis
            if obj.flag_bloodVolumeRate_analysis
                fprintf("\n----------------------------------\nBlood Volume Rate Calculation\n----------------------------------\n");
                bloodVolumeRateTimer = tic;

                bloodVolumeRate(obj.maskArtery, obj.maskVein, obj.vRMS, obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList);

                fprintf("- Blood Volume Rate calculation took: %ds\n", round(toc(bloodVolumeRateTimer)));
            end

            %% Spectral Analysis
            if obj.flag_SH_analysis && isfile(fullfile(ToolBox.PW_path, 'raw', [strcat(ToolBox.main_foldername, '_SH'), '.raw']))
                fprintf("\n----------------------------------\nSpectral Analysis\n----------------------------------\n");
                timeSpectralAnalysis = tic;

                % Import SH data
                tmpname = strcat(ToolBox.main_foldername, '_SH');
                fileID = fopen(fullfile(ToolBox.PW_path, 'raw', [tmpname, '.raw']));
                videoSH = fread(fileID, 'float32');
                fclose(fileID);

                [numX, numY, numFrames] = size(obj.M0_data_video);
                bin_x = 4; bin_y = 4; bin_t = 1;
                SH_cube = reshape(videoSH, ceil(numX / (bin_x)), ceil(numY / (bin_y)), [], ceil(numFrames / bin_t));

                % Spectrum Analysis
                fprintf("\n----------------------------------\nSpectrum Analysis\n----------------------------------\n");
                spectrumAnalysisTimer = tic;
                
                fprintf("\n----------------------------------\nSpectrum analysis\n----------------------------------\n")
                spectrum_analysis(SH_cube, obj.M0_data_video);
                
                time_spectrumAnalysis = toc(spectrumAnalysisTimer);
                fprintf("- Spectrum Analysis took : %ds\n", round(time_spectrumAnalysis))
                
                
                % Spectrogram
                fprintf("\n----------------------------------\nSpectrogram\n----------------------------------\n");
                spectrogramTimer = tic;
                spectrum_video(SH_cube, obj.maskArtery, obj.maskNeighbors);
                fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));

                fprintf("\n----------------------------------\nSpectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
            end

            %% Final Output
            tTotal = toc(totalTime);
            fprintf("\n----------------------------------\nTotal Pulsewave timing: %ds\n", round(tTotal));
            fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

            clear ToolBox;
            diary off;
            displaySuccessMsg(1);
        end
    end
end