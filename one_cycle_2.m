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

    one_cycle_video = mat2gray(mean(A1_video, 4));
end