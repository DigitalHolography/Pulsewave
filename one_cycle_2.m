function one_cycle_video = one_cycle_2(video, sys_index_list, Ninterp)
    %   one_cycle() : identifies pulse cycles and average them to one video
    %   sys_index_list : list of systole indexes in video
    arguments
        video 
        sys_index_list {mustBeNumeric, mustBeNonnegative} = []
        Ninterp {mustBeNumeric, mustBePositive} = 64
    end

    %Video retrieval frame by frame
    
    M = length(sys_index_list)-1;
    A1_video = zeros(size(video,1), size(video,2), Ninterp, M);
    for ii = 1:M
        for id_x = 1 : size(video,1)
            for id_y = 1 : size(video,2)
                interp_range = linspace(sys_index_list(ii),sys_index_list(ii+1)-1,Ninterp);
                A1_video(id_x,id_y,:,ii) = interp1((sys_index_list(ii):sys_index_list(ii+1)-1),squeeze(video(id_x, id_y,sys_index_list(ii):sys_index_list(ii+1)-1)),interp_range);
            end
        end
    end
    shift = floor(size(A1_video, 3) / 16);
    A1_video = circshift(A1_video, shift, 3);
    figure(43)
    plot(squeeze(mean(A1_video(200:300,200:300,:,:), [1 2])));

    one_cycle_video = mat2gray(mean(A1_video, 4));

%     [path, name, ext] = fileparts(video_pathname);
%     if (~exist(fullfile(path, 'one_cycle'), 'dir'))
%         mkdir(fullfile(path, 'one_cycle'));
%     end
%     V = VideoWriter(fullfile(path, 'one_cycle', [name, '_oneCycle', ext]));
%     open(V);
%     for i = 1:size(one_cycle_video,3)
%         writeVideo(V, one_cycle_video(:,:,i));
%     end
%     close(V);
end