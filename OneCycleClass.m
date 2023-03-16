classdef OneCycleClass
    properties
        data (1,:) cell % RMS M2/M0
        dataM1M0 (1,:) cell % AVG M1/M0
        dataM0 (1,:) cell % M0 raw
        dataM1 (1,:) cell % M1 raw
        dataM2 (1,:) cell % M2 raw
        directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
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
            obj.data = cell(1,obj.nbFiles) ;
            obj.dataM1M0 = cell(1,obj.nbFiles) ;
            obj.dataM0 = cell(1,obj.nbFiles) ;
            obj.dataM1 = cell(1,obj.nbFiles) ;
            obj.dataM2 = cell(1,obj.nbFiles) ;
            for ii = 1 : obj.nbFiles
                currentFilePath = fullfile(obj.directory,obj.filenames{ii});
                [filepath,name,ext] = fileparts(currentFilePath);
                if (ext == '.avi')
                    % sqrt M2/M0 : DopplerRMS
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    V = VideoReader(fullfile(path,files{ii}));
                    video = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        video(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.data{ii} = video;
                    % Import Moment 0
                    tmpname = strcat(name(1:end-2), 'moment0');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    V = VideoReader(fullfile(filepath,[tmpname,ext]));
                    videoM0 = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        videoM0(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.dataM0{ii} = videoM0;
                    % Import Moment 1
                    tmpname = strcat(name(1:end-2), 'moment1');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    V = VideoReader(fullfile(filepath,[tmpname,ext]));
                    videoM1 = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        videoM1(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.dataM1{ii} = videoM1;
                    % Import Moment 2
                    tmpname = strcat(name(1:end-2), 'moment2');
                    disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
                    V = VideoReader(fullfile(filepath,[tmpname,ext]));
                    videoM2 = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        videoM2(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.dataM2{ii} = videoM2;

                    % M1/M0 : DopplerAvg
                    directorySection = path;
                    directorySection(end-3:end) = 'avi\';
                    filenamesSection = files{ii};
                    filenamesSection = filenamesSection(1:end-6);
                    filenamesSection = strcat(filenamesSection, 'DopplerRMS.avi');
                    disp(['reading : ',strcat(directorySection, filenamesSection)]);
                    
                    V_section = VideoReader(strcat(directorySection, filenamesSection));
                    if isempty(V_section)
                        obj.video_loaded = false;
                    else
                        obj.video_loaded = true;
                    end
                    obj.videoSection = zeros(V_section.Height, V_section.Width, V_section.NumFrames);
                    for n = 1 : V_section.NumFrames
                        obj.videoSection(:,:,n) = rgb2gray(read(V_section, n));
                    end
                    obj.pictureSection = mean(mat2gray(obj.videoSection), 3);
                    


                elseif (ext == '.raw')
                    % sqrt M2/M0 : DopplerRMS
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    fileID = fopen(currentFilePath);
                    video = fread(fileID,'float32');
                    fclose(fileID);
                    obj.data{ii} = reshape(video,refvideosize);
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
                else
                    disp([filepath,name,ext,' : non recognized video format']);
                end
            end
        end

        function obj = MomentNormalize(obj)
            for ii = 1 : obj.nbFiles
                avgM0 = mean(obj.dataM0{ii},[1 2]);
                for jj = 1 : size(obj.data{ii}, 3)
                    obj.data{ii}(:,:,jj) = sqrt((obj.dataM2{ii}(:,:,jj))./avgM0(:,:,jj)); 
%                       obj.data{ii}(:,:,jj) = obj.data{ii}(:,:,jj);
                end
            end
        end

        function [sys_index_list_cell, mask_cell, maskArtery, fullPulseWave_cell] = getSystole(obj)
            sys_index_list_cell = cell(obj.nbFiles) ;
            for i = 1:obj.nbFiles
                datacube = obj.data{i}; % choix du cube sur lequel travailler 
                [sys_index_list_cell{i}, mask_cell{i}, maskArtery, fullPulseWave_cell{i}] = find_systole_index(datacube);
            end
        end

        function onePulse(obj, fullPulseWave_cell, mask_cell, maskArtery, sys_index_list_cell, Ninterp, add_infos)
            idx = 0 ;
            while (exist(fullfile(obj.directory, sprintf("one_cycle_%d",idx)), 'dir'))
                idx = idx + 1 ;
            end
            one_cycle_dir = fullfile(obj.directory, sprintf("one_cycle_%d",idx)) ;
            mkdir(one_cycle_dir);
            for n = 1:obj.nbFiles
                % FIXME maybe mask_cell{n} unnecessary. compute mask only once?
                datacube = obj.data{n}; % choix du cube sur lequel travailler 
%                 avgM0 = mean(obj.dataM0{n},[1 2]);
%                 datacube = sqrt(obj.dataM2{n}/avgM0);
                one_cycle_video = create_one_cycle(datacube, mask_cell{n}, sys_index_list_cell{n}, Ninterp) ;
                % FIXME : si l'image de depart est raw, sauver un .raw
                one_cycle_video_to_save = mat2gray(one_cycle_video); %1st normalization
                [~,name,ext] = fileparts(obj.filenames{n}) ;
                w = VideoWriter(fullfile(one_cycle_dir,strcat(name,'_one_cycle',ext))) ;
                open(w)
                % FIXME flatfield correction does not work
%                 for j = 1:size(one_cycle_video_to_save,3)
%                     one_cycle_video_to_save(:,:,j) = ...
%                         flat_field_correction(squeeze(one_cycle_video_to_save(:,:,j)),0.07,0.25);
%                 end
                one_cycle_video_to_save = mat2gray(one_cycle_video_to_save);% 2nd normalization
                for j = 1:size(one_cycle_video_to_save,3)
                    writeVideo(w,one_cycle_video_to_save(:,:,j)) ;
                end
                close(w);

                if add_infos
                    datacube = obj.data{n}; % choix du cube sur lequel travailler
%                     avgM0 = mean(obj.dataM0{n},[1 2]);
%                     datacube = sqrt(obj.dataM2{n}/avgM0);
                    [maskVein, maskArteryInPlane, maskCRA, v_RMS] = pulseAnalysis(Ninterp,datacube,obj.dataM1M0{n},one_cycle_video,one_cycle_dir,name,sys_index_list_cell{n}, mask_cell{n},maskArtery);
                    [flowVideoRGB] = flow_rate(maskArtery, maskVein, maskCRA, v_RMS, one_cycle_dir, name);

                    w = VideoWriter(fullfile(one_cycle_dir,strcat(name,'_flowVideo',ext))) ;
                    open(w)
                    flow_video_to_save = mat2gray(flowVideoRGB);% 2nd normalization
                    for jj = 1:size(flow_video_to_save,3)
                        writeVideo(w,squeeze(flow_video_to_save(:,:,jj,:))) ;
                    end
                    close(w);

                end % add_infos
            end
            displaySuccessMsg(1);
        end
    end
end