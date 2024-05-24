classdef OneCycleClass
    properties
        reference (1,:) cell % M0 AVI
        reference_norm (1,:) cell % M0 norm AVI
        reference_interp (1,:) cell % M0 AVI
        reference_norm_interp (1,:) cell % M0 norm AVI
        dataM2M0 (1,:) cell % RMS M2/M0
        dataM1M0 (1,:) cell % AVG M1/M0
        dataM2M0_interp (1,:) cell % RMS M2/M0
        dataM1M0_interp (1,:) cell % AVG M1/M0
        dataM0 (1,:) cell % M0 raw
        dataM0_interp (1,:) cell % M0 raw
        dataM1 (1,:) cell % M1 raw
        dataM1_interp (1,:) cell % M1 raw
        dataM2 (1,:) cell % M2 raw
        dataM2_interp (1,:) cell % M2 raw
        dataSH (1,:) cell % SH raw
        dataSH_interp (1,:) cell % SH raw
        directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
        k double %interpolaton parameter
        load_logs char
        % For sectioning
        video_loaded
        pictureSection
        videoSection
        flag_FlatField
        flag_SH_analysis
        flag_PulseWave_analysis
        isdone_flatfield
        %% FIXME : relancer a chaque rendering
        ToolBoxmaster ToolBoxClass
        
        
        
        
        
    end
    
    methods
        function obj = OneCycleClass(path)
            
            arguments
                path
            end
            
            obj.directory = path;
            
            tmp_idx = regexp(path, '\');
            obj.filenames = {path(tmp_idx(end-1)+1:end-1)};
            
            obj.nbFiles = 1 ;
            obj.reference = cell(1,obj.nbFiles) ;
            obj.reference_norm = cell(1,obj.nbFiles);
            obj.dataM2M0 = cell(1,obj.nbFiles) ;
            obj.dataM1M0 = cell(1,obj.nbFiles) ;
            obj.dataM2M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM0 = cell(1,obj.nbFiles) ;
            obj.dataM0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1 = cell(1,obj.nbFiles) ;
            obj.dataM1_interp = cell(1,obj.nbFiles) ;
            obj.dataM2 = cell(1,obj.nbFiles) ;
            obj.dataM2_interp = cell(1,obj.nbFiles) ;
            obj.dataSH = cell(1,obj.nbFiles) ;
            obj.dataSH_interp = cell(1,obj.nbFiles) ;
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
            obj.isdone_flatfield = 0;
            obj.ToolBoxmaster =  ToolBoxClass(obj.directory);
            
            obj.load_logs = '\n=== LOADING : \r\n';
            logs = obj.load_logs;
            
            for ii = 1 : obj.nbFiles
                
                %% Ref loading
                dir_path_avi = fullfile(obj.directory, 'avi');
                NameRefAviFile = strcat(obj.filenames{ii}, '_M0');
                RefAviFilePath = fullfile(dir_path_avi, NameRefAviFile);
                ext = '.avi';
                
                disp(['reading : ',RefAviFilePath]);
                RefAviFilePath(strfind(RefAviFilePath,'\'))='/';
                str_tosave = sprintf('reading : %s', RefAviFilePath);
                logs = strcat(logs,'\r', str_tosave);
                
                V = VideoReader(fullfile(dir_path_avi,[NameRefAviFile,ext]));
                video = zeros(V.Height, V.Width, V.NumFrames);
                for n = 1 : V.NumFrames
                    video(:,:,n) = rgb2gray(read(V, n));
                end
                obj.reference{ii} = video;
                refvideosize = size(obj.reference{ii});
                
                %% File loading
                dir_path_raw = fullfile(obj.directory ,'raw');
                if isempty(dir(dir_path_raw))
                    fprintf("No Raw Files, please select a Holowave folder with raw files exported")
                end
                NameRawFile = strcat(obj.filenames{ii}, '_moment0');
                ext = '.raw';
                
                % Import Moment 0
                disp(['reading : ',fullfile(dir_path_raw,[NameRawFile,ext])]);
                FilePathUnix = fullfile(dir_path_raw,[NameRawFile,ext]);
                FilePathUnix(strfind(FilePathUnix,'\'))='/';
                str_tosave = sprintf('reading : %s', FilePathUnix);
                logs = strcat(logs,'\r', str_tosave);
                
                fileID = fopen(fullfile(dir_path_raw,[NameRawFile,ext]));
                videoM0 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM0{ii} = reshape(videoM0,refvideosize);
                % Import Moment 1
                NameRawFile = strcat(obj.filenames{ii}, '_moment1');
                disp(['reading : ',fullfile(dir_path_raw,[NameRawFile,ext])]);
                FilePathUnix = fullfile(dir_path_raw,[NameRawFile,ext]);
                FilePathUnix(strfind(FilePathUnix,'\'))='/';
                str_tosave = sprintf('reading : %s', FilePathUnix);
                logs = strcat(logs,'\r', str_tosave);
                
                fileID = fopen(fullfile(dir_path_raw,[NameRawFile,ext]));
                videoM1 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM1{ii} = reshape(videoM1,refvideosize);
                % Import Moment 2
                NameRawFile = strcat(obj.filenames{ii}, '_moment2');
                disp(['reading : ',fullfile(dir_path_raw,[NameRawFile,ext])]);
                FilePathUnix = fullfile(dir_path_raw,[NameRawFile,ext]);
                FilePathUnix(strfind(FilePathUnix,'\'))='/';
                str_tosave = sprintf('reading : %s', FilePathUnix);
                logs = strcat(logs,'\r', str_tosave);
                
                fileID = fopen(fullfile(dir_path_raw,[NameRawFile,ext]));
                videoM2 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM2{ii} = reshape(videoM2,refvideosize);
                
                %                 % Import SH
                %
                %                 NameRawFile = strcat(obj.filenames{ii}, '_SH');
                %                 disp(['reading : ',fullfile(dir_path_raw,[NameRawFile,ext])]);
                %                 FilePathUnix = fullfile(dir_path_raw,[NameRawFile,ext]);
                %                 FilePathUnix(strfind(FilePathUnix,'\'))='/';
                %                 str_tosave = sprintf('reading : %s', FilePathUnix);
                %                 logs = strcat(logs,'\r', str_tosave);
                %
                %                 fileID = fopen(fullfile(dir_path_raw,[NameRawFile,ext]));
                %                 videoSH = fread(fileID,'float32');
                %                 fclose(fileID);
                %                 obj.dataSH{ii} = reshape(videoSH,refvideosize(1),refvideosize(2),[]);
                %
                
                obj.load_logs = logs;
                
            end
        end
        
        function obj = MomentNormalize(obj)
            for ii = 1 : obj.nbFiles
                avgM0 = mean(obj.dataM0{ii},[1 2]);
                avgRef = mean(obj.reference{ii},[1 2]);
                
                tmp_M2M0 = obj.dataM0{ii};
                tmp_M2 = obj.dataM2{ii};
                for jj = 1 : size(obj.dataM0{ii}, 3)
                    tmp_M2M0(:,:,jj) = sqrt((tmp_M2(:,:,jj))./avgM0(:,:,jj));
                end
                obj.dataM2M0{ii} = tmp_M2M0;
                clear tmp_M2M0 tmp_M2
                
                tmp_M1M0 = obj.dataM0{ii};
                tmp_M1 = obj.dataM1{ii};
                for jj = 1 : size(obj.dataM0{ii}, 3)
                    tmp_M1M0(:,:,jj) = (tmp_M1(:,:,jj))./avgM0(:,:,jj);
                end
                obj.dataM1M0{ii} = tmp_M1M0;
                clear tmp_M1M0 tmp_M1
                
                tmp_Ref = obj.reference{ii};
                for jj = 1 : size(obj.dataM0{ii}, 3)
                    tmp_Ref(:,:,jj) = tmp_Ref(:,:,jj)./avgRef(:,:,jj);
                end
                obj.reference_norm{ii} = tmp_Ref;
                clear tmp_Ref
            end
        end
        
        function obj = MomentFlatField(obj)
            PW_params = Parameters_json(obj.directory);
            for ii = 1 : obj.nbFiles
                height = size(obj.dataM0_interp{1}, 1);
                width = size(obj.dataM0_interp{1}, 2);
                length = size(obj.dataM0_interp{1}, 3);
                num_frames = size(obj.reference{1}, 3);
                gwRatio = PW_params.flatField_gwRatio;
                flatField_border = PW_params.flatField_border;
                
                tmp_dataM0 = zeros(height, width, length);
                tmp_calc_data_M0 = obj.dataM0_interp{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM0(:,:,i) =  flat_field_correction_old(tmp_calc_data_M0(:,:,i), gwRatio*height, flatField_border);
                end
                obj.dataM0_interp{1} = tmp_dataM0;
                clear tmp_dataM0 tmp_calc_data_M0
                
                tmp_dataM1 = zeros(height, width, length);
                tmp_calc_data_M1 = obj.dataM1_interp{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM1(:,:,i) =  flat_field_correction_old(tmp_calc_data_M1(:,:,i), gwRatio*height, flatField_border);
                end
                obj.dataM1_interp{1} = tmp_dataM1;
                clear tmp_dataM1 tmp_calc_data_M1
                
                tmp_dataM2 = zeros(height, width, length);
                tmp_calc_data_M2 = obj.dataM2_interp{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM2(:,:,i) =  flat_field_correction_old(tmp_calc_data_M2(:,:,i), gwRatio*height, flatField_border);
                end
                obj.dataM2_interp{1} = tmp_dataM2;
                clear tmp_dataM2 tmp_calc_data_M2
                
                tmp_dataM2M0 = zeros(height, width, length);
                tmp_calc_data = (obj.dataM2M0_interp{1});
                parfor i = 1 : num_frames
                    tmp_dataM2M0(:,:,i) = flat_field_correction_old(tmp_calc_data(:,:,i), gwRatio*height, flatField_border);
                end
                obj.dataM2M0_interp{1} = single(tmp_dataM2M0);
                clear tmp_dataM2M0 tmp_calc_data
                
            end
        end
        
        function obj = cropAllVideo(obj)
            %Crop a video (matrix dim 3)
            PW_params = Parameters_json(obj.directory);
            firstFrame = PW_params.videoStartFrameIndex;
            lastFrame = PW_params.videoEndFrameIndex;
            
            logs = obj.load_logs;
            
            if firstFrame > 0 && firstFrame < size(obj.dataM0{1},3) && lastFrame > firstFrame && lastFrame <= size(obj.dataM0{1},3)
                
                
                obj.reference{1} = obj.reference{1}(:,:,firstFrame:lastFrame);
                obj.reference_norm{1} = obj.reference{1}(:,:,firstFrame:lastFrame);
                obj.dataM0{1} = obj.dataM0{1}(:,:,firstFrame:lastFrame);
                obj.dataM1{1} = obj.dataM1{1}(:,:,firstFrame:lastFrame);
                obj.dataM2{1} = obj.dataM2{1}(:,:,firstFrame:lastFrame);
                
                
                disp(['Data cube frame: ',num2str(firstFrame), '/', num2str(size(obj.dataM0{1},3)), ' to ',num2str(lastFrame), '/', num2str(size(obj.dataM0{1},3))])
                
                str_tosave = sprintf('Data cube frame: %s/%s to %s/%s',num2str(firstFrame), num2str(size(obj.dataM0{1},3)),num2str(lastFrame), num2str(size(obj.dataM0{1},3)));
                logs = strcat(logs,'\r', str_tosave);
            else
                disp('Wrong value for the first frame. Set as 1.')
                disp('Wrong value for the last frame. Set as the end.')
                disp(['Data cube frame: 1/', num2str(size(obj.dataM0{1},3)), ' to ',num2str(size(obj.dataM0{1},3)), '/', num2str(size(obj.dataM0{1},3))])
                
                str_tosave = sprintf('Wrong value for the first frame. Set as 1. \rWrong value for the last frame. Set as the end. \rData cube frame: 1/%s to %s/%s', num2str(size(obj.dataM0{1},3)),num2str(size(obj.dataM0{1},3)), num2str(size(obj.dataM0{1},3)));
                logs = strcat(logs,'\r\n\n', str_tosave,'\n');
            end
            
            obj.load_logs = logs;
            
        end
        
        function obj = Interpolate(obj) %ref = TRUE indicates the object is the reference
            height = size(obj.reference{1}, 1);
            width = size(obj.reference{1}, 2);
            num_frames = size(obj.reference{1}, 3);
            k_interp = obj.k;
            height = (height-1)*(2^k_interp-1)+height;
            width = (width-1)*(2^k_interp-1)+width;
            
            tmp_ref = zeros(height,width, size(obj.dataM0{1}, 3));
            tmp_calc_ref = obj.reference{1};
            parfor i = 1 : num_frames
                tmp_ref(:,:,i) = interp2(tmp_calc_ref(:,:,i), k_interp);
            end
            obj.reference_interp{1} = tmp_ref;
            clear tmp_ref tmp_calc_ref
            
            tmp_ref_norm = zeros(height,width, size(obj.dataM0{1}, 3));
            tmp_calc_ref_norm = obj.reference_norm{1};
            parfor i = 1 : num_frames
                tmp_ref_norm(:,:,i) = interp2(tmp_calc_ref_norm(:,:,i), k_interp);
            end
            obj.reference_norm_interp{1} = tmp_ref_norm;
            clear tmp_ref_norm tmp_calc_ref_norm
            
            tmp_dataM2M0 = zeros(height, width, size(obj.dataM0{1}, 3));
            tmp_calc_data = obj.dataM2M0{1};
            parfor i = 1 : num_frames
                tmp_dataM2M0(:,:,i) = interp2(tmp_calc_data(:,:,i), k_interp);
            end
            obj.dataM2M0_interp{1} = single(tmp_dataM2M0);
            clear tmp_dataM2M0 tmp_calc_data
            
            tmp_dataM0 = zeros(height, width, size(obj.dataM0{1}, 3));
            tmp_calc_data_M0 = obj.dataM0{1};
            parfor i = 1 : num_frames % loop over frames
                tmp_dataM0(:,:,i) = interp2(tmp_calc_data_M0(:,:,i), k_interp);
            end
            obj.dataM0_interp{1} = tmp_dataM0;
            clear tmp_dataM0 tmp_calc_data_M0
            
            tmp_dataM1 = zeros(height, width, size(obj.dataM1{1}, 3));
            tmp_calc_data_M1 = obj.dataM1{1};
            parfor i = 1 : num_frames % loop over frames
                tmp_dataM1(:,:,i) = interp2(tmp_calc_data_M1(:,:,i), k_interp);
            end
            obj.dataM1_interp{1} = tmp_dataM1;
            clear tmp_dataM1 tmp_calc_data_M1
            
            tmp_dataM2 = zeros(height, width, size(obj.dataM2{1}, 3));
            tmp_calc_data_M2 = obj.dataM2{1};
            parfor i = 1 : num_frames % loop over frames
                tmp_dataM2(:,:,i) = interp2(tmp_calc_data_M2(:,:,i), k_interp);
            end
            obj.dataM2_interp{1} = tmp_dataM2;
            clear tmp_dataM2 tmp_calc_data_M2
            
            tmp_data_M1M0 = zeros(height, width, size(obj.dataM0{1}, 3));
            tmp_calc_data_M1M0 = obj.dataM1M0{1};
            parfor i = 1 : num_frames
                tmp_data_M1M0(:,:,i) = interp2(tmp_calc_data_M1M0(:,:,i), k_interp);
            end
            obj.dataM1M0_interp{1} = tmp_data_M1M0;
            clear tmp_data_M1M0 tmp_calc_data_M1M0
            
            % obj.reference{1} = []; %need n frames ; maybe stock in obj ?
            obj.reference_norm{1} = [];
            obj.dataM2M0{1} = [];
            obj.dataM0{1} = [];
            obj.dataM1M0{1} = [];
            
            %             if obj.flag_SH == 1
            %                 tmp_data_SH = zeros(height, width, size(obj.dataSH{1}, 3));
            %                 tmp_calc_data_SH = obj.dataSH{1};
            %                 for i = 1 : size(obj.dataSH{1}, 3)
            %                     tmp_data_SH(:,:,i) = interp2(tmp_calc_data_SH(:,:,i), k_interp);
            %                 end
            %                 obj.dataSH_interp{1} = tmp_data_SH;
            %
            %             end
            
        end
        
        function onePulse(obj, Ninterp)
            % PW_params = Parameters_json(obj.directory);
            % ToolBox = obj.ToolBoxmaster;
            checkPulsewaveParamsFromJson(obj.directory);
            PW_params = Parameters_json(obj.directory);
            
            obj.k = PW_params.k;
            obj.ToolBoxmaster =  ToolBoxClass(obj.directory);
            ToolBox = obj.ToolBoxmaster;
            
            meanIm = squeeze(mean(obj.reference_interp{1}, 3)); % Because highest intensities in CRA usually
            meanIm_M1M0 = squeeze(mean(obj.dataM1M0_interp{1}, 3)); % Because velocities coming from the CRA are out-of-plane
            blurred_mask = imgaussfilt(double(meanIm.*meanIm_M1M0),PW_params.gauss_filt_size_for_barycentre*size(meanIm.*meanIm_M1M0,1),'Padding',0);
            [ToolBox.y_barycentre,ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask,[],'all'));
            
            clear blurred_mask
            
            mkdir(ToolBox.PW_path_dir);
            mkdir(ToolBox.PW_path_png);
            mkdir(ToolBox.PW_path_eps);
            mkdir(ToolBox.PW_path_gif);
            mkdir(ToolBox.PW_path_txt);
            mkdir(ToolBox.PW_path_avi);
            mkdir(ToolBox.PW_path_mp4);
            mkdir(ToolBox.PW_path_json);
            mkdir(ToolBox.PW_path_log);
            
            %             path_dir_txt = fullfile(obj.directory,'txt');
            %             path_file_txt_params = fullfile(path_dir_txt,'InputPulsewaveParams.txt');
            %             copyfile(path_file_txt_params,ToolBox.PW_path_txt );
            
            path_dir_json = fullfile(obj.directory,'pulsewave','json');
            path_file_json_params = fullfile(path_dir_json,'InputPulsewaveParams.json');
            copyfile(path_file_json_params,ToolBox.PW_path_json );
            
            %saving times
            path_file_txt_exe_times = fullfile(ToolBox.PW_path_log, sprintf('%s_execution_times.txt',ToolBox.PW_folder_name));
            fileID = fopen(path_file_txt_exe_times,'w');
            fprintf(fileID,'EXECUTION TIMES : \r\n==================\n\r\n');
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
            
            
            fileID = fopen(path_file_log,'w') ;
            
            fprintf(fileID,'==================\rGIT VERSION :\r');
            fprintf(fileID,MessBranch, resultBranch);
            fprintf(fileID,MessHash, resultHash);
            fprintf(fileID,'==================\r\n ');
            
            fprintf(fileID, obj.load_logs);
            
            fprintf(fileID, '\r\n=== EXECUTION :\r\n\n');
            
            fclose(fileID);
            
            %% FlatField
            
            if obj.flag_FlatField && (obj.isdone_flatfield == 0)
                
                fileID = fopen(path_file_txt_exe_times,'a+') ;
                fprintf(fileID,'FLATFIELD : \r\n\n');
                fclose(fileID);
                
                disp("Correcting data with flat field")
                fileID = fopen(path_file_log,'a+') ;
                fprintf(fileID, 'Correcting data with flat field \r\n');
                fclose(fileID);
                
                for ii = 1 : obj.nbFiles
                    height = size(obj.dataM0_interp{1}, 1);
                    width = size(obj.dataM0_interp{1}, 2);
                    length = size(obj.dataM0_interp{1}, 3);
                    num_frames = size(obj.reference{1}, 3);
                    gwRatio = PW_params.flatField_gwRatio;
                    flatField_border = PW_params.flatField_border;
                    
                    tic
                    tmp_dataM0 = zeros(height, width, length);
                    tmp_calc_data_M0 = obj.dataM0_interp{1};
                    parfor i = 1 : num_frames % loop over frames
                        tmp_dataM0(:,:,i) =  flat_field_correction_old(tmp_calc_data_M0(:,:,i), gwRatio*height, flatField_border);
                    end
                    obj.dataM0_interp{1} = tmp_dataM0;
                    clear tmp_dataM0 tmp_calc_data_M0
                    timeM0 = toc;
                    save_time(path_file_txt_exe_times, 'FlatField M0', timeM0)
                    
                    tic
                    tmp_dataM2M0 = zeros(height, width, length);
                    tmp_calc_data = obj.dataM2M0_interp{1};
                    parfor i = 1 : num_frames
                        tmp_dataM2M0(:,:,i) = flat_field_correction_old(tmp_calc_data(:,:,i), gwRatio*height, flatField_border);
                    end
                    obj.dataM2M0_interp{1} = single(tmp_dataM2M0);
                    clear tmp_dataM2M0 tmp_calc_data
                    timeM2M0 = toc;
                    save_time(path_file_txt_exe_times, 'FlatField M2M0', timeM2M0)
                    
                    fileID = fopen(path_file_txt_exe_times,'a+') ;
                    fprintf(fileID,'\r\n=== Total : %.0fs \r\n\n----------\r\n',timeM2M0+timeM0);
                    fclose(fileID);
                end
                obj.isdone_flatfield = 1;
            else
                disp("no flatfield applied")
                fileID = fopen(path_file_log,'a+') ;
                fprintf(fileID, 'No flatfield applied \r\n\n');
                fclose(fileID);
            end
            
            %% Creating Masks
            
            fileID = fopen(path_file_txt_exe_times,'a+') ;
            fprintf(fileID,'\r\n');
            fclose(fileID);
            
            tic
            [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV,maskSectionArtery] = createMasks(obj.reference_norm_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            disp('CreatMasks timing :')
            time = toc;
            disp(time)
            
            fileID = fopen(path_file_log,'a+') ;
            fprintf(fileID, 'CreatMasks timing : % \r\n\n');
            fclose(fileID);
            
            save_time(path_file_txt_exe_times, 'createMasks', time)
            fileID = fopen(path_file_txt_exe_times,'a+') ;
            fprintf(fileID,'\r\n----------\r\n');
            fclose(fileID);
            
            tic
            % [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV,maskSectionArtery]= createMasks(obj.reference_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            % [mask_artery, mask_vein, mask_vessel,mask_background,mask_CRA,maskCRV] = createMasks(obj.reference_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            
            % disp('CreatMasks New timing :')
            % time = toc;
            % disp(time)
            %save_time(path_file_txt_exe_times, 'createMasksNew', time)
            
            %% PulseWave Analysis
            
            if obj.flag_PulseWave_analysis
                close all
                
                sys_index_list_cell = cell(obj.nbFiles) ;
                fullPulseWave_cell = cell(obj.nbFiles) ;
                for n = 1:obj.nbFiles
                    fileID = fopen(path_file_txt_exe_times,'a+') ;
                    fprintf(fileID,'PULSEWAVE ANALYSIS : \r\n\n');
                    fclose(fileID);
                    
                    time_sys_idx=0;
                    time_pulseanalysis=0;
                    time_vel_map=0;
                    time_hist=0;
                    time_arterial_res=0;
                    time_flowrate =0;
                    
                    tic
                    % [sys_index_list_cell{n}, fullPulseWave_cell{n}] = find_systole_index(obj.reference_norm_interp{n}, obj.directory,maskArtery);
                    [sys_index_list_cell{n}, fullPulseWave_cell{n}] = find_systole_index(obj.reference_interp{n}, obj.directory,maskArtery);
                    disp('FindSystoleIndex timing :')
                    time_sys_idx = toc;
                    disp(time_sys_idx)
                    save_time(path_file_txt_exe_times, 'Find Systole Index', time_sys_idx)
                    
                    [v_RMS_one_cycle,v_RMS_all,onePulseVideoM0,exec_times, total_time] = pulseAnalysis(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},obj.dataM2_interp{n},obj.dataM0_interp{n},obj.reference_interp{1},sys_index_list_cell{n},meanIm, maskArtery,maskVein,maskBackground ,ToolBox,obj.directory);
                    disp('PulseAnalysis timing :')
                    time_pulseanalysis = total_time;
                    disp(time_pulseanalysis)
                    save_time(path_file_txt_exe_times, 'Pulse Analysis', time_pulseanalysis)
                    %exec time details
                    fileID = fopen(path_file_txt_exe_times,'a+') ;
                    for i=1:size(exec_times,2)
                        fprintf(fileID,'\t%s : %.0fs \r\n',exec_times(1,i),exec_times(2,i));
                    end
                    fclose(fileID);
                    clear exec_times
                    
                    tic
                    velocity_map(maskArtery, maskVein, v_RMS_one_cycle, onePulseVideoM0, ToolBox); %FIXME Histo en trop, jsute pour faire la FLOWMAP ?
                    disp('Velocity map timing :')
                    time_vel_map = toc;
                    disp(time_vel_map)
                    save_time(path_file_txt_exe_times, 'Velocity map', time_vel_map)
                    
                    tic
                    ImgM0 = rescale(mean(obj.reference_interp{1},3));
                    velocityHistogramm(v_RMS_all, maskArtery,maskVein ,ImgM0 , ToolBox, path)
                    disp('Velocity Histogramm timing :')
                    time_hist = toc;
                    disp(time_hist)
                    save_time(path_file_txt_exe_times, 'Velocity Histogramm', time_hist)
                    
                    
                    %                     tic
                    %                     BKGHistogramm(obj.dataM2M0_interp{1}, maskBackground ,ToolBox)
                    %                     disp('BKG Histogramm timing :')
                    %                     toc
                    
                    
                    
                    tic
                    ArterialResistivityIndex(v_RMS_one_cycle, obj.reference_interp{1}, maskArtery,  ToolBox);
                    disp('ArterialResistivityIndex timing :')
                    time_arterial_res = toc;
                    disp(time_arterial_res)
                    save_time(path_file_txt_exe_times, 'Arterial Resistivity Index', time_arterial_res)
                    
                    
                    %[v] = pulseAnalysisTest(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},obj.reference_interp{1},sys_index_list_cell{n},maskArtery,maskVessel,maskVein,maskBackground ,ToolBox,obj.directory);
                    %                     detectElasticWave(datacube, maskArtery, maskCRA);
                    tic
                    % try
                    disp('flowrate :')
                    flow_rate(maskArtery, maskVein, maskCRA, maskSectionArtery, v_RMS_one_cycle, obj.dataM0_interp{1}, obj.reference_interp{n}, ToolBox, obj.k,obj.directory);
                    disp('FlowRate timing :')
                    time_flowrate = toc;
                    disp(time_flowrate)
                    save_time(path_file_txt_exe_times, 'Flow rate', time_flowrate)
                    % catch
                    % disp('no flow rate')
                    
                    %                     try
                    %                         flow_rate_old(maskArtery, maskVein, maskCRA, v_RMS, ToolBox, obj.k,obj.directory);
                    %                     catch
                    %                     end
                    
                    
                    
                    %                         tic
                    %                         spectrogram(maskArtery,maskBackground,  obj.dataSH_interp{n}, ToolBox);
                    %                         disp('Spectrogram timing :')
                    %                         time = toc;
                    %                         disp(time)
                    %                         save_time(path_file_txt_exe_times, '\tSpectrogram', time)
                    %                     else
                    
                    fileID = fopen(path_file_txt_exe_times,'a+') ;
                    fprintf(fileID,'\r\n=== Total : %.0fs \r\n\n----------\r\n',time_sys_idx+time_pulseanalysis+time_vel_map+time_hist+time_arterial_res+time_flowrate);
                    fclose(fileID);
                    % end
                end
            end
            
            %% Spectrum Analysis
            
            if obj.flag_SH_analysis
                
                %% Import SH
                cubeSize = 256 ;
                cubeFreqLength = 32 ;
                tmpname  = strcat(ToolBox.main_foldername, '_SH');
                ext = '.raw';
                disp(['reading : ',fullfile(obj.directory,'raw',[tmpname,ext])]);
                fileID = fopen(fullfile(obj.directory,'raw',[tmpname,ext]));
                videoSH = fread(fileID,'float32');
                fclose(fileID);
                SH_cube=  reshape(videoSH,cubeSize,cubeSize,cubeFreqLength,[]);
                
                tic
                spectrum_analysis(maskArtery,maskBackground ,SH_cube ,ToolBox,obj.dataM0_interp{1});
                disp('Spectrum analysis :')
                time = toc;
                disp(time)
                save_time(path_file_txt_exe_times, 'spectrum_analysis', time)
                
                tic
                spectrogram(maskArtery,maskBackground ,SH_cube ,ToolBox);
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