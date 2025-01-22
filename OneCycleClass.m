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

        sysIdxList % list of frame indexes counting carciac cycles
        xy_barycenter % x y position of the ONH
        vRMS % video estimate of velocity map in retinal vessels

        directory char % directory of input data (from HoloDoppler or HoloVibes)
        PW_params_names cell  % filenames of all the current input parameters ('InputPulsewaveParams.json' for example by default))
        PW_param_name char % current filename
        filenames char % name id used for storing the measured rendered data

        flag_Segmentation
        flag_SH_analysis
        flag_PulseWave_analysis
        flag_velocity_analysis
        flag_ExtendedPulseWave_analysis
        flag_bloodVolumeRate_analysis
        flag_bloodVelocityProfile_analysis

        ToolBoxmaster ToolBoxClass
        OverWrite logical

    end

    methods

        function obj = OneCycleClass(path)
            % creation of an instance of the main class
            % input is the rendered measurement data
            % path is either a directory (from HoloDoppler) or a .holo file path (from HoloVibes)
            % rendered data mainly consists of the short time Doppler spectrum Moments videos

            if ~isfolder(path) % if the input path is a .holo file with moments inside
                [filepath, name, ~] = fileparts(path);
                if ~ isfolder(fullfile(filepath, name)) % creates a result folder in the location as the file
                    mkdir(fullfile(filepath, name));
                end
                obj.directory = fullfile(filepath, name);
                obj.filenames = name;
            else % if input is a directory from HoloDoppler
                obj.directory = path;
                tmp_idx = regexp(path, '\');
                obj.filenames = obj.directory(tmp_idx(end - 1) + 1:end -1);
            end

            %% Parameters

            % looks for existing parameters and checks compatibility between found PW params and Default PW params of this version of PW.
            obj.PW_params_names = checkPulsewaveParamsFromJson(obj.directory);
            obj.PW_param_name = obj.PW_params_names{1}; % default behavior takes the first found parameter file ('InputPulseWaveParameters.json')

            %% Video loading

            if ~isfolder(path) % if holo file with moments inside is the input
                disp(['reading moments in : ', strcat(obj.directory, '.holo')]);
                [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
                readMomentsFooter(obj.directory);
                obj.M0_ff_video = pagetranspose(improve_video(ff_correction(videoM0, 35),0.0005, 2, 0));
                obj.M0_data_video = pagetranspose(videoM0);
                obj.M1_data_video = pagetranspose(videoM1/1e3); % rescale because M1 in Hz in Holovibes and in kHz in Pulsewave
                obj.M2_data_video = pagetranspose(videoM2/1e6); % rescale because M2 in Hz² in Holovibes and in kHz² in Pulsewave
                
            else
                obj = readRaw(obj);
                
            end

            obj.is_preprocessed = false;

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

            obj.is_preprocessed = true;
        end

        function obj = onePulse(obj)
            %  ------- This is the app main routine. --------

            obj.ToolBoxmaster = ToolBoxClass(obj.directory,obj.PW_param_name,obj.OverWrite);
            setGlobalToolBox(obj.ToolBoxmaster)
            ToolBox = getGlobalToolBox;
            PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);

            totalTime = tic;

            saveGit;
            
            %% Creating Masks

            if obj.flag_Segmentation

                createMasksTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Mask Creation\n")
                fprintf("----------------------------------\n")

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.maskArtery, obj.maskVein, obj.maskSection, obj.xy_barycenter] = createMasks(obj.M0_ff_video, f_AVG_mean);

                time_create_masks = toc(createMasksTimer);
                fprintf("- Mask Creation took : %ds\n", round(time_create_masks))

            end

            %% PulseWave Analysis

            if obj.flag_PulseWave_analysis

                findSystoleTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Find Systole\n")
                fprintf("----------------------------------\n")

                [obj.sysIdxList, ~] = find_systole_index(obj.M0_ff_video, obj.maskArtery);

                time_sys_idx = toc(findSystoleTimer);
                fprintf("- FindSystoleIndex took : %ds\n", round(time_sys_idx))

                pulseAnalysisTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Pulse Analysis\n")
                fprintf("----------------------------------\n")

                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.vRMS] = pulseAnalysis(obj.f_RMS_video, obj.maskArtery, obj.maskVein, obj.maskSection);

                if obj.flag_ExtendedPulseWave_analysis
                extendedPulseAnalysis(obj.M0_ff_video, obj.f_RMS_video, f_AVG_mean, obj.vRMS, obj.maskArtery, obj.maskVein, obj.maskSection, obj.sysIdxList);
                end

                time_pulseanalysis = toc(pulseAnalysisTimer);
                fprintf("- Pulse Analysis took : %ds\n", round(time_pulseanalysis))
            end

            % pulseVelocityTimer = tic;
            %
            % fprintf("\n----------------------------------\n")
            % fprintf("Pulse Velocity\n")
            % fprintf("----------------------------------\n")
            %
            % pulseVelocity(obj.M0_data_video, maskArtery)
            %
            % time_pulsevelocity = toc(pulseVelocityTimer);
            % fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))

            if obj.flag_velocity_analysis
                bloodFlowVelocityTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Blood Flow Velocity Calculation\n")
                fprintf("----------------------------------\n")

                bloodFlowVelocity(obj.vRMS, obj.maskArtery, obj.maskVein, obj.maskSection, obj.M0_ff_video, obj.xy_barycenter)

                time_velo = toc(bloodFlowVelocityTimer);
                fprintf("- Blood Flow Velocity calculation took : %ds\n", round(time_velo))
            end

            if obj.flag_bloodVolumeRate_analysis
                bloodVolumeRateTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Blood Volume Rate Calculation\n")
                fprintf("----------------------------------\n")

                bloodVolumeRate(obj.maskArtery, obj.maskVein, obj.vRMS, obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList, obj.flag_bloodVelocityProfile_analysis);

                time_volumeRate = toc(bloodVolumeRateTimer);
                fprintf("- Blood Volume rate calculation took : %ds\n", round(time_volumeRate))
            end


            %% Spectral Analysis

            if obj.flag_SH_analysis && isfile(fullfile(ToolBox.PW_path, 'raw', [strcat(ToolBox.main_foldername, '_SH'), '.raw']))

                % Import SH

                timeSpectralAnalysis = tic;

                tmpname = strcat(ToolBox.main_foldername, '_SH');
                ext = '.raw';
                disp(['reading : ', fullfile(ToolBox.PW_path, 'raw', [tmpname, ext])]);
                fileID = fopen(fullfile(obj.directory, 'raw', [tmpname, ext]));

                k = PW_params.k;

                videoSH = fread(fileID, 'float32');
                fclose(fileID);
                [numX, numY, numFrames] = size(obj.f_RMS_video);
                bin_x = 4;
                bin_y = 4;
                % bin_w = 16;
                bin_t = 1;

                SH_cube = reshape(videoSH, ceil(numX / (2 ^ k * bin_x)), ceil(numY / (2 ^ k * bin_y)), [], ceil(numFrames / bin_t));

                % Spectrum Analysis
                spectrumAnalysisTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Spectrum analysis\n")
                fprintf("----------------------------------\n")

                spectrum_analysis(SH_cube, obj.M0_data_video);

                time_spectrumAnalysis = toc(spectrumAnalysisTimer);
                fprintf("- Spectrum Analysis took : %ds\n", round(time_spectrumAnalysis))


                % Spectrogram
                spectrogramTimer = tic;

                fprintf("\n----------------------------------\n")
                fprintf("Spectrogram\n")
                fprintf("----------------------------------\n")

                spectrogram(obj.maskArtery, obj.maskSection, SH_cube);

                time_spectrogram = toc(spectrogramTimer);
                fprintf("- Sprectrogram took : %ds\n", round(time_spectrogram))


                tSpectralAnalysis = toc(timeSpectralAnalysis);

                fprintf("\n----------------------------------\n")
                fprintf("Spectral Analysis timing : %ds\n", round(tSpectralAnalysis))
            end

            tTotal = toc(totalTime);

            fprintf("\n----------------------------------\n")
            fprintf("Total Pulsewave timing : %ds\n", round(tTotal))
            fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'))

            clear ToolBox
            diary off
            displaySuccessMsg(1);

        end

    end

end
