classdef OneCycleClass
    properties
        data (1,:) cell
        directory char
        filenames (1,:) cell
        nbFiles {mustBeNumeric , mustBePositive}  
    end

    methods
        function obj = OneCycleClass(files,path)
            obj.directory = path;
            obj.filenames = files ;
            obj.nbFiles = length(files) ; 
            obj.data = cell(1,obj.nbFiles) ;
            for i = 1 : obj.nbFiles
                V = VideoReader(fullfile(path,files{i})) ;
                video = zeros(V.Height, V.Width, V.NumFrames);
                for n = 1 : V.NumFrames
                    video(:,:,n) = rgb2gray(read(V, n));
                end
                obj.data{i} = video;
            end
        end

        function sys_index_list_cell = getSystole(obj)
            sys_index_list_cell = cell(obj.nbFiles) ;
            for i = 1:obj.nbFiles
                sys_index_list_cell{i} = find_systole_index(obj.data{i});
            end
        end

        function writeCut(obj,sys_index_list_cell, Ninterp, add_infos)
            idx = 0 ;
            while (exist(fullfile(obj.directory, sprintf("one_cycle_%d",idx)), 'dir'))
                idx = idx + 1 ;
            end
            one_cycle_dir = fullfile(obj.directory, sprintf("one_cycle_%d",idx)) ; 
            mkdir(one_cycle_dir);
            for n = 1 : obj.nbFiles
                one_cycle_video = create_one_cycle(obj.data{n}, sys_index_list_cell{n}, Ninterp) ;
                [~,name,ext]=fileparts(obj.filenames{n}) ;
                w = VideoWriter(fullfile(one_cycle_dir,strcat(name,'_one_cycle',ext))) ;
                open(w)
                for j = 1:size(one_cycle_video,3)
                    writeVideo(w,one_cycle_video(:,:,j)) ;
                end
                close(w);
                if add_infos
                    SegmentationAV(one_cycle_video,one_cycle_dir, name) ; 
                end 
            end
        end
    end
end