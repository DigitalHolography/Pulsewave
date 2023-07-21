classdef OneCycleClass
    properties
        reference (1,:) cell % M0 AVI
        reference_interp (1,:) cell % M0 AVI
        dataM2M0 (1,:) cell % RMS M2/M0
        dataM1M0 (1,:) cell % AVG M1/M0
        dataM2M0_interp (1,:) cell % RMS M2/M0
        dataM1M0_interp (1,:) cell % AVG M1/M0
        dataM0 (1,:) cell % M0 raw
        dataM0_interp (1,:) cell % M0 raw
        dataM1 (1,:) cell % M1 raw
        dataM2 (1,:) cell % M2 raw
        dataSH (1,:) cell % SH raw
        dataSH_interp (1,:) cell % SH raw
        directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
        k double %interpolaton parameter
        % For sectioning
        video_loaded
        pictureSection
        videoSection


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
            obj.dataM2M0 = cell(1,obj.nbFiles) ;
            obj.dataM1M0 = cell(1,obj.nbFiles) ;
            obj.dataM2M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1M0_interp = cell(1,obj.nbFiles) ;
            obj.dataM0 = cell(1,obj.nbFiles) ;
            obj.dataM0_interp = cell(1,obj.nbFiles) ;
            obj.dataM1 = cell(1,obj.nbFiles) ;
            obj.dataM2 = cell(1,obj.nbFiles) ;
            obj.dataSH = cell(1,obj.nbFiles) ;
            obj.dataSH_interp = cell(1,obj.nbFiles) ;
            obj.k = 0; 

            checkPulsewaveParamsFromTxt(obj.directory);
            PW_params = Parameters(obj.directory);
            obj.k = PW_params.k;

            for ii = 1 : obj.nbFiles

                %% Ref loading
                dir_path_avi = fullfile(obj.directory, 'avi');
                NameRefAviFile = strcat(obj.filenames{ii}, '_M0');
                RefAviFilePath = fullfile(dir_path_avi, NameRefAviFile);
                ext = '.avi';

                disp(['reading : ',RefAviFilePath]);
                V = VideoReader(fullfile(dir_path_avi,[NameRefAviFile,ext]));
                video = zeros(V.Height, V.Width, V.NumFrames);
                for n = 1 : V.NumFrames
                    video(:,:,n) = rgb2gray(read(V, n));
                end
                obj.reference{ii} = video;
                refvideosize = size(obj.reference{ii});

                %% File loading
                dir_path_raw = fullfile(obj.directory ,'raw');
                RawFilePath = fullfile(dir_path_raw,strcat(obj.filenames{ii}, '_DopplerRMS.raw'));
                NameRawFile = strcat(obj.filenames{ii}, '_DopplerRMS');
                ext = '.raw';

                % sqrt M2/M0 : DopplerRMS
                disp(['reading : ',fullfile(dir_path_raw,NameRawFile)]);
                fileID = fopen(RawFilePath);
                video = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM2M0{ii} = reshape(video,refvideosize);
                % M1/M0 : DopplerAvg
                tmpname = NameRawFile;
                tmpname(end-2:end) = 'AVG';
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM1M0 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM1M0{ii} = reshape(videoM1M0,refvideosize);
                % Import Moment 0
                tmpname = strcat(NameRawFile(1:end-10), 'moment0');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM0 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM0{ii} = reshape(videoM0,refvideosize);
                % Import Moment 1
                tmpname = strcat(NameRawFile(1:end-10), 'moment1');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM1 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM1{ii} = reshape(videoM1,refvideosize);
                % Import Moment 2
                tmpname = strcat(NameRawFile(1:end-10), 'moment2');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoM2 = fread(fileID,'float32');
                fclose(fileID);
                obj.dataM2{ii} = reshape(videoM2,refvideosize);
                % Import SH
                tmpname = strcat(NameRawFile(1:end-10), 'SH');
                disp(['reading : ',fullfile(dir_path_raw,[tmpname,ext])]);
                fileID = fopen(fullfile(dir_path_raw,[tmpname,ext]));
                videoSH = fread(fileID,'float32');
                fclose(fileID);
                obj.dataSH{ii} = reshape(videoSH,refvideosize(1),refvideosize(2),[]);


            end
        end

        function obj = MomentNormalize(obj)
            for ii = 1 : obj.nbFiles
                avgM0 = mean(obj.dataM0{ii},[1 2]);
                avgRef = mean(obj.reference{ii},[1 2]);


                %donnée temporaire
                tmp_M2M0 = obj.dataM2M0{ii};
                tmp_M2 = obj.dataM2{ii};
                tmp_M1 = obj.dataM1{ii};
                tmp_M1M0 = obj.dataM1M0{ii};
                tmp_Ref = obj.reference{ii};

                for jj = 1 : size(obj.dataM2M0{ii}, 3)
                    tmp_M2M0(:,:,jj) = sqrt((tmp_M2(:,:,jj))./avgM0(:,:,jj));
                    tmp_M1M0(:,:,jj) = (tmp_M1(:,:,jj))./avgM0(:,:,jj);
                    tmp_Ref(:,:,jj) = tmp_Ref(:,:,jj)./avgRef(:,:,jj);

                end
                obj.dataM2M0{ii} = tmp_M2M0;
                obj.dataM1M0{ii} = tmp_M1M0;
                obj.reference{ii} = tmp_Ref;
            end
        end

        function obj = cropAllVideo(obj)
            %Crop a video (matrix dim 3)
            PW_params = Parameters(obj.directory);
            firstFrame = PW_params.videoStartFrameIndex;
            lastFrame = PW_params.videoEndFrameIndex;

            if firstFrame > 0 && firstFrame < size(obj.dataM2M0{1},3) && lastFrame > firstFrame && lastFrame <= size(obj.dataM2M0{1},3)
                
                
                obj.reference{1} = obj.reference{1}(:,:,firstFrame:lastFrame);
                obj.dataM2M0{1} = obj.dataM2M0{1}(:,:,firstFrame:lastFrame);
                obj.dataM1M0{1} = obj.dataM1M0{1}(:,:,firstFrame:lastFrame);
                obj.dataM0{1} = obj.dataM0{1}(:,:,firstFrame:lastFrame);
                obj.dataM1{1} = obj.dataM1{1}(:,:,firstFrame:lastFrame);
                obj.dataM2{1} = obj.dataM2{1}(:,:,firstFrame:lastFrame);
                

                disp(['Data cube frame: ',num2str(firstFrame), '/', num2str(size(obj.dataM2M0{1},3)), ' to ',num2str(lastFrame), '/', num2str(size(obj.dataM2M0{1},3))])
            else
                disp('Wrong value for the first frame. Set as 1')
                disp('Wrong value for the last frame. Set as the end.')
                disp(['Data cube frame: 1/', num2str(size(obj.dataM2M0{1},3)), ' to ',num2str(size(obj.dataM2M0{1},3)), '/', num2str(size(obj.dataM2M0{1},3))])
            end
        end

        function obj = Interpolate(obj) %ref = TRUE indicates the object is the reference
            height = size(obj.reference{1}, 1);
            width = size(obj.reference{1}, 2);
            num_frames = size(obj.reference{1}, 3);
            k_interp = obj.k;
            height = (height-1)*(2^k_interp-1)+height;
            width = (width-1)*(2^k_interp-1)+width;

      
                tmp_ref = zeros(height,width, size(obj.dataM2M0{1}, 3));
                tmp_calc_ref = obj.reference{1}; 
                tmp_dataM0 = zeros(height, width, size(obj.dataM2M0{1}, 3));
                tmp_dataM2M0 = zeros(height, width, size(obj.dataM2M0{1}, 3));
                tmp_data_M1M0 = zeros(height, width, size(obj.dataM2M0{1}, 3));
                tmp_data_SH = zeros(height, width, size(obj.dataSH{1}, 3));
                tmp_calc_data = obj.dataM2M0{1};
                tmp_calc_data_M0 = obj.dataM0{1};
                tmp_calc_data_M1M0 = obj.dataM1M0{1};
                tmp_calc_data_SH = obj.dataSH{1};
                parfor i = 1 : num_frames % loop over frames
                    tmp_dataM0(:,:,i) = interp2(tmp_calc_data_M0(:,:,i), k_interp);
                    tmp_dataM2M0(:,:,i) = interp2(tmp_calc_data(:,:,i), k_interp);
                    tmp_data_M1M0(:,:,i) = interp2(tmp_calc_data_M1M0(:,:,i), k_interp);
                    tmp_ref(:,:,i) = interp2( tmp_calc_ref(:,:,i), k_interp);
                end
                for i = 1 : size(obj.dataSH{1}, 3)
                    tmp_data_SH(:,:,i) = interp2(tmp_calc_data_SH(:,:,i), k_interp);
                end
                obj.reference_interp{1} = tmp_ref;
                obj.dataM2M0_interp{1} = tmp_dataM2M0;
                obj.dataM1M0_interp{1} = tmp_data_M1M0;
                obj.dataM0_interp{1} = tmp_dataM0;
                obj.dataSH_interp{1} = tmp_data_SH;


            
        end

%         function [sys_index_list_cell, mask_cell,maskVessel, maskArtery, fullPulseWave_cell] = getSystole(obj)
%             PW_params = Parameters(obj.global_directory);
%             sys_index_list_cell = cell(obj.nbFiles) ;
%             for i = 1:obj.nbFiles
%                 datacube = obj.dataM2M0_interp{i}; % choix du cube sur lequel travailler 
%                 [sys_index_list_cell{i}, mask_cell{i},maskVessel, maskArtery, fullPulseWave_cell{i}] = find_systole_index(datacube, obj.global_directory);
%             end
%         end

        function onePulse(obj, Ninterp, add_infos)
            PW_params = Parameters(obj.directory);
            ToolBox = ToolBoxClass(obj.directory);
          
            meanIm = squeeze(mean(obj.reference_interp{1}, 3));
            blurred_mask = imgaussfilt(double(meanIm),0.05*size(meanIm,1),'Padding',0);
            [ToolBox.y_barycentre,ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask,[],'all'));

            mkdir(ToolBox.PW_path_dir);
            mkdir(ToolBox.PW_path_png);
            mkdir(ToolBox.PW_path_eps);
            mkdir(ToolBox.PW_path_txt);
            mkdir(ToolBox.PW_path_avi);
            mkdir(ToolBox.PW_path_mp4);
 
            path_dir_txt = fullfile(obj.directory,'txt');
            path_file_txt_params = fullfile(path_dir_txt,'InputPulsewaveParams.txt');
            copyfile(path_file_txt_params,ToolBox.PW_path_txt );


            tic
            [maskArtery, maskVein, maskVessel,maskBackground,maskCRA,maskCRV] = createMasks(obj.reference_interp{1} ,obj.dataM1M0_interp{1}, obj.directory, ToolBox);
            disp('CreatMasks timing :')
            toc

            if add_infos
                close all

                sys_index_list_cell = cell(obj.nbFiles) ;
                fullPulseWave_cell = cell(obj.nbFiles) ;
                for n = 1:obj.nbFiles
                    tic
                    [sys_index_list_cell{n}, fullPulseWave_cell{n}] = find_systole_index(obj.reference_interp{n}, obj.directory,maskArtery);
                    disp('FindSystoleIndex timing :')
                    toc
                    %                 % FIXME flatfield correction does not work
                    % %                 for j = 1:size(one_cycle_video_to_save,3)
                    % %                     one_cycle_video_to_save(:,:,j) = ...
                    % %                         flat_field_correction(squeeze(one_cycle_video_to_save(:,:,j)),0.07,0.25);
                    % %                 end
                    tic
                    [v_RMS] = pulseAnalysis(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},sys_index_list_cell{n},maskArtery,maskVein,maskBackground ,ToolBox,obj.directory);
                    disp('PulseAnalysis timing :')
                    toc
                    
                    tic
                    ArterialResistivityIndex(v_RMS,obj.dataM2M0_interp{n}, maskArtery,  ToolBox);
                    disp('ArterialResistivityIndex timing :')
                    toc
                    %[v] = pulseAnalysisTest(Ninterp,obj.dataM2M0_interp{n},obj.dataM1M0_interp{n},sys_index_list_cell{n},maskArtery,maskVessel,maskVein,maskBackground ,ToolBox,obj.directory);
                    %                     detectElasticWave(datacube, maskArtery, maskCRA);
                    tic
                    flow_rate(maskArtery, maskVein, maskCRA, v_RMS, ToolBox, obj.k,obj.directory);
                    disp('FlowRate timing :')
                    toc

                    tic
                    spectrogram(maskArtery,maskBackground,  obj.dataSH_interp{n}, ToolBox);
                    disp('Spectrogram timing :')
                    toc
                    % add_infos
                end
            end
            displaySuccessMsg(1);
        end
    end
end