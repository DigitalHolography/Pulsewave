classdef ExecutionClass < handle
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
        params_names cell % filenames of all the current input parameters ('InputEyeFlowParams.json' for example by default)
        param_name char % current filename
        filenames char % name id used for storing the measured rendered data

        flag_Segmentation
        flag_SH_analysis
        flag_Pulse_analysis
        flag_velocity_analysis
        flag_bloodVolumeRate_analysis

        OverWrite logical
        ToolBox ToolBoxClass
    end

    methods
        function obj = ExecutionClass(path)
            % Constructor for ExecutionClass.
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
            obj.params_names = checkEyeFlowParamsFromJson(obj.directory);
            obj.param_name = obj.params_names{1}; % Default behavior

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


        function obj = analyzeData(obj, app)
            % Main routine for EyeFlow analysis.

            % Initialize ToolBox and parameters
            TB = obj.ToolBox;
            params = TB.getParams;
            totalTime = tic;
            saveGit;

            %% Mask Creation
            if obj.flag_Segmentation
                fprintf("\n----------------------------------\nMask Creation\n----------------------------------\n");
                createMasksTimer = tic;

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.maskArtery, obj.maskVein, obj.maskSection, obj.maskNeighbors, obj.xy_barycenter] = ...
                    createMasks(obj.M0_ff_video, f_AVG_mean);

                M0_ff_img = rescale(mean(obj.M0_ff_video, 3));
                cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
                cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
                cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

                M0_Artery = setcmap(M0_ff_img, obj.maskArtery, cmapArtery);
                M0_Vein = setcmap(M0_ff_img, obj.maskVein, cmapVein);
                M0_AV = setcmap(M0_ff_img, obj.maskArtery & obj.maskVein, cmapAV);

                M0_RGB = (M0_Artery + M0_Vein) .* ~(obj.maskArtery & obj.maskVein) + M0_AV + rescale(M0_ff_img) .* ~(obj.maskArtery | obj.maskVein);
                app.ImageDisplay.ImageSource = mat2gray(M0_RGB); % Rescale the image for display
                ax = ancestor(app.ImageDisplay, 'axes');
                axis(ax, 'equal');

                fprintf("- Mask Creation took: %ds\n", round(toc(createMasksTimer)));
            end

            %% Pulse Analysis
            if obj.flag_Pulse_analysis
                fprintf("\n----------------------------------\nFind Systole\n----------------------------------\n");
                findSystoleTimer = tic;

                [obj.sysIdxList, ~, sysMaxList, sysMinList] = find_systole_index(obj.M0_ff_video, obj.maskArtery);

                % Check if the output vectors are long enough
                if numel(obj.sysIdxList) < 2 || numel(sysMaxList) < 2 || numel(sysMinList) < 2
                    warning('There isnt enough systoles.');
                else
                    % Log systole results
                    fileID = fopen(fullfile(TB.path_txt, strcat(TB.main_foldername, '_EF_main_outputs.txt')), 'a');
                    fprintf(fileID, 'Heart beat: %f (bpm) \r\n', 60 / mean(diff(obj.sysIdxList) * TB.stride / TB.fs / 1000));
                    fprintf(fileID, 'Systole Indices: %s \r\n', strcat('[', sprintf("%d,", obj.sysIdxList), ']'));
                    fprintf(fileID, 'Number of Cycles: %d \r\n', numel(obj.sysIdxList) - 1);
                    fprintf(fileID, 'Max Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMaxList), ']'));
                    fprintf(fileID, 'Min Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMinList), ']'));
                    fprintf(fileID, 'Time diastolic min to systolic max derivative (ms): %f \r\n', ...
                        1000 * mean((obj.sysIdxList(2:end) - sysMinList) * TB.stride / TB.fs / 1000));
                    fprintf(fileID, 'Time diastolic min to systolic max (ms): %f \r\n', ...
                        1000 * mean((sysMaxList(2:end) - sysMinList(1:end - 1)) * TB.stride / TB.fs / 1000));
                    fclose(fileID);
                end

                fprintf("- FindSystoleIndex took: %ds\n", round(toc(findSystoleTimer)));

                fprintf("\n----------------------------------\nPulse Analysis\n----------------------------------\n");
                pulseAnalysisTimer = tic;

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                obj.vRMS = pulseAnalysis(obj.f_RMS_video, obj.maskArtery, obj.maskVein, obj.maskSection, obj.maskNeighbors);

                if params.json.PulseAnalysis.ExtendedFlag
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
            if obj.flag_SH_analysis && isfile(fullfile(TB.EF_path, 'raw', [strcat(TB.main_foldername, '_SH'), '.raw']))
                fprintf("\n----------------------------------\nSpectral Analysis\n----------------------------------\n");
                timeSpectralAnalysis = tic;

                % Import SH data
                tmpname = strcat(TB.main_foldername, '_SH');
                fileID = fopen(fullfile(TB.EF_path, 'raw', [tmpname, '.raw']));
                videoSH = fread(fileID, 'float32');
                fclose(fileID);

                [numX, numY, numFrames] = size(obj.M0_data_video);
                bin_x = 4; bin_y = 4; bin_t = 1;
                SH_cube = reshape(videoSH, ceil(numX / (bin_x)), ceil(numY / (bin_y)), [], ceil(numFrames / bin_t));

                %% Spectrum Analysis
                fprintf("\n----------------------------------\nSpectrum Analysis\n----------------------------------\n");
                spectrumAnalysisTimer = tic;

                spectrum_analysis(SH_cube, obj.M0_data_video);

                time_spectrumAnalysis = toc(spectrumAnalysisTimer);
                fprintf("- Spectrum Analysis took : %ds\n", round(time_spectrumAnalysis))

                %% Spectrogram
                fprintf("\n----------------------------------\nSpectrogram\n----------------------------------\n");
                spectrogramTimer = tic;

                spectrum_video(SH_cube, obj.maskArtery, obj.maskNeighbors);

                fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));

                fprintf("\n----------------------------------\nSpectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
            end

            %% Final Output
            tTotal = toc(totalTime);
            fprintf("\n----------------------------------\nTotal EyeFlow timing: %ds\n", round(tTotal));
            fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

            clear TB;
            diary off;
            displaySuccessMsg(1);
        end
    end
end