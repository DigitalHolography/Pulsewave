classdef OneCycleClass
    properties
        data (1,:) cell
        directory char
        filenames (1,:) cell

    end

    methods
        function obj = OneCycleClass(files,path)
            obj.directory = path;
            obj.filenames = files ;
            obj.data = cell(1,length(files)) ;
            for i = 1 : length(files)
                V = VideoReader(fullfile(path,files{i})) ;
                tmp = im2double(read(V, [1 Inf]));
                video = zeros(size(tmp, 1), size(tmp, 2), size(tmp, 4));
                for frame = 1 : size(video, 3)
                    video(:,:,frame) = squeeze(rgb2gray(tmp(:,:,:, frame)));
                end
                obj.data{i} = video;
            end
        end

        function sys_index_list_cell = getSystole(obj)
            sys_index_list_cell = cell(length(obj.data)) ;
            for i = 1:length(obj.data)
                sys_index_list_cell{i} = find_systole_index(obj.data{i});
            end
        end

        function writeCut(obj,sys_index_list_cell,Ninterp)
            for n = 1 : length(obj.data)
                one_cycle_video = one_cycle_2(obj.data{n},sys_index_list_cell{n},Ninterp) ;
                if (~exist(fullfile(obj.directory, 'one_cycle'), 'dir'))
                    mkdir(fullfile(obj.directory, 'one_cycle'));
                end
                [~,name,ext]=fileparts(obj.filenames{n}) ;
                w = VideoWriter(fullfile(obj.directory,'one_cycle',strcat(name,'_one_cycle',ext))) ;
                open(w)
                for j = 1:size(one_cycle_video,3)
                    writeVideo(w,one_cycle_video(:,:,j)) ;
                end
                close(w)

            end
        end
    end
end