classdef OneCycleClass
    properties
        data (1,:) cell
        dataM1M0 (1,:) cell
        directory char
        nbFiles {mustBeNumeric , mustBePositive}
        filenames (1,:) cell
        % For sectioning
<<<<<<< HEAD
        video_loaded
        pictureSection
        videoSection
=======
        directorySection char
        filenamesSection char
        pictureSection 
        videoSection 
        video_loaded
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0


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
            for ii = 1 : obj.nbFiles
                currentFilePath = fullfile(obj.directory,obj.filenames{ii});
                [filepath,name,ext] = fileparts(currentFilePath);
                if (ext == '.avi')
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    V = VideoReader(fullfile(path,files{ii}));
                    video = zeros(V.Height, V.Width, V.NumFrames);
                    for n = 1 : V.NumFrames
                        video(:,:,n) = rgb2gray(read(V, n));
                    end
                    obj.data{ii} = video;
                    %% importation fichier M0 pour pulse section
<<<<<<< HEAD
                    directorySection = path;
                    directorySection(end-3:end) = 'avi\';
                    filenamesSection = files{ii};

                    filenamesSection = filenamesSection(1:end-14);
                    filenamesSection = strcat(filenamesSection, 'DopplerRMS.avi');
                    disp(['reading : ',strcat(directorySection, filenamesSection)]);
                    
                    V_section = VideoReader(strcat(directorySection, filenamesSection));
=======
                    obj.directorySection = path;
                    obj.directorySection(end-3:end) = 'avi\';
                    obj.filenamesSection = files;

                    obj.filenamesSection = obj.filenamesSection(1:end-14);
                    obj.filenamesSection = strcat(obj.filenamesSection, 'DopplerRMS.avi');
                    disp(['reading : ',strcat(obj.directorySection, obj.filenamesSection)]);
                    
                    V_section = VideoReader(strcat(obj.directorySection, obj.filenamesSection));
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0
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
<<<<<<< HEAD
                    
=======
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0


                elseif (ext == '.raw')
                    % sqrt M2/M0 : DopplerRMS
                    disp(['reading : ',fullfile(filepath,[name,ext])]);
                    fileID = fopen(currentFilePath);
                    video = fread(fileID,'float32');
                    fclose(fileID);
                    obj.data{ii} = reshape(video,refvideosize);
%                     % M1/M0 : DopplerAvg
%                     tmpname = name;
%                     tmpname(end-2:end) = 'AVG';
%                     disp(['reading : ',fullfile(filepath,[tmpname,ext])]);
%                     fileID = fopen(fullfile(filepath,[tmpname,ext]));
%                     videoM1M0 = fread(fileID,'float32');
%                     fclose(fileID);
%                     obj.dataM1M0{ii} = reshape(videoM1M0,refvideosize);
                else
                    disp([filepath,name,ext,' : non recognized video format']);
                end
            end
        end

        function [sys_index_list_cell, mask_cell, fullPulseWave_cell] = getSystole(obj)
            sys_index_list_cell = cell(obj.nbFiles) ;
            for i = 1:obj.nbFiles
                [sys_index_list_cell{i}, mask_cell{i}, fullPulseWave_cell{i}] = find_systole_index(obj.data{i});
            end
        end

        function onePulse(obj, fullPulseWave_cell, mask_cell, sys_index_list_cell, Ninterp, add_infos)
            idx = 0 ;
            while (exist(fullfile(obj.directory, sprintf("one_cycle_%d",idx)), 'dir'))
                idx = idx + 1 ;
            end
            one_cycle_dir = fullfile(obj.directory, sprintf("one_cycle_%d",idx)) ;
            mkdir(one_cycle_dir);
            for n = 1:obj.nbFiles
                % FIXME maybe mask_cell{n} unnecessary. compute mask only once?
                one_cycle_video = create_one_cycle(obj.data{n}, mask_cell{n}, sys_index_list_cell{n}, Ninterp) ;
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
                    pulseAnalysis(Ninterp,obj.data{n},obj.dataM1M0{n},one_cycle_video,one_cycle_dir,name,sys_index_list_cell{n}, mask_cell{n});
                end
            end
            displaySuccessMsg(1);
        end
    end
end