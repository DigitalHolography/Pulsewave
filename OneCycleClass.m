classdef OneCycleClass
    properties
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
        global_directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
        k double %interpolaton parameter
        % For sectioning
        video_loaded
        pictureSection
        videoSection


    end

    methods
        function obj = OneCycleClass(files,path,refvideosize)
            arguments
                files
                path
                refvideosize = []
            end
            obj.directory = path;
            obj.filenames = files;
            obj.nbFiles = length(files) ;
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

            dir_path_cell = split(obj.directory,'\');
            dir_path_cell = dir_path_cell(1:(size(dir_path_cell)-2));
            dir_path_cell = join(dir_path_cell,'\');
            obj.global_directory = dir_path_cell{1};

            for ii = 1 : obj.nbFiles
                currentFilePath = fullfile(obj.directory,obj.filenames{ii});
                [filepath,name,ext] = fileparts(currentFilePath);

                if (ext == '.avi')
                    %FIXME : refvideosize & double appel de obj.data{i}
                    % dans la boucle avi et raw
                    
                    % sqrt M2/M0 : DopplerRMS
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    V = VideoReader(fullfile(path,files{ii}));
                    video = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        video(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.dataM2M0{ii} = video;
                    %obj.k = 1;
                    checkPulsewaveParamsFromTxt(obj.global_directory);
                    PW_params = Parameters(obj.global_directory);
                    obj.k = PW_params.k;
                   % [~,obj.k] =  getPulsewaveParamFromTxt(path,'Value of the interpolation parameter :');
                    
                    
                elseif (ext == '.raw')
                    % sqrt M2/M0 : DopplerRMS
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    fileID = fopen(currentFilePath);
                    video = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataM2M0{ii} = reshape(video,refvideosize);
                    % M1/M0 : DopplerAvg
                    tmpname = name;
                    tmpname(end-2:end) = 'AVG';
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    fileID = fopen(fullfile(filepath,[tmpname,ext]));
                    videoM1M0 = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataM1M0{ii} = reshape(videoM1M0,refvideosize);
                    % Import Moment 0
                    tmpname = strcat(name(1:end-10), 'moment0');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    fileID = fopen(fullfile(filepath,[tmpname,ext]));
                    videoM0 = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataM0{ii} = reshape(videoM0,refvideosize);
                    % Import Moment 1
                    tmpname = strcat(name(1:end-10), 'moment1');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    fileID = fopen(fullfile(filepath,[tmpname,ext]));
                    videoM1 = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataM1{ii} = reshape(videoM1,refvideosize);
                    % Import Moment 2
                    tmpname = strcat(name(1:end-10), 'moment2');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    fileID = fopen(fullfile(filepath,[tmpname,ext]));
                    videoM2 = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataM2{ii} = reshape(videoM2,refvideosize);
                     % Import SH
                    tmpname = strcat(name(1:end-10), 'SH');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    fileID = fopen(fullfile(filepath,[tmpname,ext]));
                    videoSH = fread(fileID,'float32');
                    fclose(fileID);
                    obj.dataSH{ii} = reshape(videoSH,refvideosize(1),refvideosize(2),[]);

                    checkPulsewaveParamsFromTxt(obj.global_directory);
                    PW_params = Parameters(obj.global_directory);
                    obj.k = PW_params.k;
   %                 [~,obj.k] =  getPulsewaveParamFromTxt(path,'Value of the interpolation parameter :');
                  
   
                else
                    disp([filepath,name,ext,' : non recognized video format']);
                end
            end
        end

        function obj = MomentNormalize(obj)
            for ii = 1 : obj.nbFiles
                avgM0 = mean(obj.dataM0{ii},[1 2]);


                %donnée temporaire
                tmp_data = obj.dataM2M0{ii};
                tmp_M2 = obj.dataM2{ii};
                tmp_M1 = obj.dataM1{ii};
                tmp_M1M0 = obj.dataM1M0{ii};

                for jj = 1 : size(obj.dataM2M0{ii}, 3)
                    tmp_data(:,:,jj) = sqrt((tmp_M2(:,:,jj))./avgM0(:,:,jj));
                    tmp_M1M0(:,:,jj) = (tmp_M1(:,:,jj))./avgM0(:,:,jj);
%                       obj.data{ii}(:,:,jj) = obj.data{ii}(:,:,jj);
                end
                obj.dataM2M0{ii} = tmp_data;
                obj.dataM1M0{ii} = tmp_M1M0;
            end
        end

        function obj = Interpolate(obj,height,width,ref) %ref = TRUE indicates the object is the reference
            
            num_frames = size(obj.dataM2M0{1}, 3);
            k_interp = obj.k;
            height = (height-1)*(2^k_interp-1)+height;
            width = (width-1)*(2^k_interp-1)+width;

            if ref 
                tmp_ref = zeros(height,width, size(obj.dataM2M0{1}, 3));
                tmp_calc_ref = obj.dataM2M0{1}; 
                
                for i = 1 : num_frames % loop over frames
                    tmp_ref(:,:,i) = interp2( tmp_calc_ref(:,:,i), k_interp);
                end
                obj.dataM2M0_interp{1} = tmp_ref;


            else
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
                end
                for i = 1 : size(obj.dataSH{1}, 3)
                    tmp_data_SH(:,:,i) = interp2(tmp_calc_data_SH(:,:,i), k_interp);
                end
                obj.dataM2M0_interp{1} = tmp_dataM2M0;
                obj.dataM1M0_interp{1} = tmp_data_M1M0;
                obj.dataM0_interp{1} = tmp_dataM0;
                obj.dataSH_interp{1} = tmp_data_SH;

            end
            
        end

        function [sys_index_list_cell, mask_cell,maskVessel, maskArtery, fullPulseWave_cell] = getSystole(obj)
            sys_index_list_cell = cell(obj.nbFiles) ;
            for i = 1:obj.nbFiles
                datacube = obj.dataM2M0_interp{i}; % choix du cube sur lequel travailler 
                [sys_index_list_cell{i}, mask_cell{i},maskVessel, maskArtery, fullPulseWave_cell{i}] = find_systole_index(datacube, obj.global_directory);
            end
        end

        function onePulse(obj, fullPulseWave_cell, mask_cell, maskVessel, maskArtery, sys_index_list_cell, Ninterp, add_infos, k)
            % k corresponds to interpolation 2^k-1
            idx = 0 ;
            [~,file_name,~] = fileparts(obj.global_directory);

            path_dir_txt = fullfile(obj.global_directory,'txt');
            path_file_txt_params = fullfile(path_dir_txt,'InputPulsewaveParams.txt');

            path_dir_pulswave = fullfile(obj.global_directory,'pulsewave');

            if ~exist(path_dir_pulswave,'dir')
                mkdir(path_dir_pulswave);
            end
            folder_name = strcat( file_name, '_pulsewave');
            while (exist(fullfile(path_dir_pulswave, sprintf('%s_%d', folder_name,idx)), 'dir'))
                idx = idx + 1 ;
            end

            %FIXME Remplacer One_cycle_dir par une cell qui contient les
            %chemin ves tous les dossiers
            one_cycle_dir = fullfile(path_dir_pulswave, sprintf('%s_%d', folder_name,idx)) ;
            mkdir(one_cycle_dir);
            mkdir(fullfile(one_cycle_dir, 'png'));
            mkdir(fullfile(one_cycle_dir, 'eps'));
            mkdir(fullfile(one_cycle_dir, 'avi'));
            mkdir(fullfile(one_cycle_dir, 'txt'));
            one_cycle_dir_png = fullfile(one_cycle_dir, 'png');
            one_cycle_dir_eps = fullfile(one_cycle_dir, 'eps');
            one_cycle_dir_txt = fullfile(one_cycle_dir, 'txt');
            one_cycle_dir_avi = fullfile(one_cycle_dir, 'avi');

            copyfile(path_file_txt_params,one_cycle_dir_txt );

            PW_params = Parameters(obj.global_directory);

%             maskBackground = not(maskVessel);
%             maskVein = double(maskVessel) - double(maskArtery); 
%             maskVein = maskVein > 0; 
% 
% 
%             imwrite(mat2gray(single(maskArtery)),fullfile(one_cycle_dir_png,'maskArtery.png'),'png') ;
%             imwrite(mat2gray(single(maskVein)),fullfile(one_cycle_dir_png,'maskVein.png'),'png') ;
%             imwrite(mat2gray(single(maskVessel)),fullfile(one_cycle_dir_png,'maskVessel.png'),'png') ;
%             imwrite(mat2gray(single(maskBackground)),fullfile(one_cycle_dir_png,'maskBackground.png'),'png') ;


            for n = 1:obj.nbFiles
                % FIXME maybe mask_cell{n} unnecessary. compute mask only once?
                datacube = obj.dataM2M0_interp{n}; % choix du cube sur lequel travailler 
%                 avgM0 = mean(obj.dataM0{n},[1 2]);
%                 datacube = sqrt(obj.dataM2{n}/avgM0);
%                 one_cycle_video = create_one_cycle(datacube, mask_cell{n}, sys_index_list_cell{n}, Ninterp) ;
%                 % FIXME : si l'image de depart est raw, sauver un .raw
%                 one_cycle_video_to_save = mat2gray(one_cycle_video); %1st normalization
                [~,name,~] = fileparts(obj.filenames{n}) ;
%                 w = VideoWriter(fullfile(one_cycle_dir_avi,strcat(name,'_one_cycle'))) ;
%                 open(w)
%                 % FIXME flatfield correction does not work
% %                 for j = 1:size(one_cycle_video_to_save,3)
% %                     one_cycle_video_to_save(:,:,j) = ...
% %                         flat_field_correction(squeeze(one_cycle_video_to_save(:,:,j)),0.07,0.25);
% %                 end
%                 one_cycle_video_to_save = mat2gray(one_cycle_video_to_save);% 2nd normalization
%                 for j = 1:size(one_cycle_video_to_save,3)
%                     writeVideo(w,one_cycle_video_to_save(:,:,j)) ;
%                 end
%                 close(w);

                if add_infos
                    close all
                    datacube = obj.dataM2M0_interp{n}; % choix du cube sur lequel travailler
                    datacube_freq = obj.dataSH{n};
%                     avgM0 = mean(obj.dataM0{n},[1 2]);
%                     datacube = sqrt(obj.dataM2{n}/avgM0);
%                      [maskVein, maskArteryInPlane, maskCRA, v_RMS] = pulseAnalysis(Ninterp,datacube,obj.dataM1M0_interp{n},one_cycle_video,one_cycle_dir,name,sys_index_list_cell{n}, mask_cell{n},maskArtery);
                    [maskVein, maskArteryInPlane, maskCRA, v_RMS] = pulseAnalysis(Ninterp,datacube,obj.dataM0_interp{n},obj.dataM1M0_interp{n},[],one_cycle_dir,name,sys_index_list_cell{n}, mask_cell{n},maskArtery,maskVessel,obj.global_directory);
%                     detectElasticWave(datacube, maskArtery, maskCRA);
                    tic
                   [flowVideoRGB] = flow_rate(maskArtery, maskVein, maskCRA, v_RMS, one_cycle_dir, name, obj.k,obj.global_directory);
                    toc
                  [SpectrogramVideo] = spectrogram(maskArtery, maskVein, maskCRA, obj.dataSH_interp{1}, one_cycle_dir, name, obj.k);

                end % add_infos
            end
            displaySuccessMsg(1);
        end
    end
end