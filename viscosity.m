function video = viscosity(mask_cell , velocity_cell, angle_list, one_cycle_dir, filename)

num_vessels = size(mask_cell, 2);
n_interp = 100;
%interpolation parameter
k = 2;

velocity_profiles = zeros(n_interp, size(velocity_cell{1}, 3), num_vessels);

for ii = 1 : num_vessels
    subImg = mask_cell{ii};
    subVideo = velocity_cell{ii};
    tilt_angle = (angle_list(ii));
    
    interp_size = 4*size(subImg, 1)-3;
    % subImg_interp = zeros(interp_size, interp_size);
    subVideo_interp = zeros(interp_size, interp_size, size(subVideo, 3));
    
    %% interpolate
    subImg_interp = interp2(subImg, k);
    for frame = 1 : size(subVideo, 3)
        subVideo_interp(:,:, frame) = interp2(subVideo(:,:, frame), k);
    end

    subImg_interp = imrotate(subImg_interp, tilt_angle ,'bilinear','crop');
    for frame = 1 : size(subVideo, 3)
        subVideo_interp(:,:, frame) = imrotate(subVideo_interp(:,:, frame), tilt_angle ,'bilinear','crop');
    end
    

    % find vessels width
    projImg = squeeze(sum(subImg_interp, 1));
    projImg = projImg/max(projImg, [], "all");
    projImg = (projImg > 0.1);
    projVideo = squeeze(sum(subVideo_interp, 1)/size(subVideo_interp,1));

    width_masked = nnz(projImg);
    list = find(projImg);
    
    % projVideo_cropped = projVideo(list(1) : list(end), :);
    
    x = 1:size(projVideo,1);
    x_wall2wall = list(1) : list(end);
    xinterp_wall2wall = linspace(list(1),list(end),n_interp);

    velocityProfileInterp = zeros(n_interp, size(projVideo,2));
    for tt = 1:size(projVideo,2)
        velocityProfileInterp(:,tt) = interp1(x, projVideo(:,tt), xinterp_wall2wall);
    end

    % figure(1234)
    % plot(projImg)

    velocity_profiles(:, :, ii) = velocityProfileInterp;
    1;

end% ii (artery #)

average_velocity_profile = squeeze(mean(velocity_profiles, 3));

v = VideoWriter(fullfile(one_cycle_dir, strcat(filename,'_velocity_profile.avi')));
open(v);
mimin = min(average_velocity_profile(:));
mamax = max(average_velocity_profile(:));
for tt = 1 : size(velocity_cell{1}, 3)
    fifig = figure(899);
    plot(squeeze(average_velocity_profile(:,tt)),'-k', 'LineWidth',2) ;
    title('average wall-to-wall arterial velocity profile');
    % legend('arterial pulse','background', 'venous signal') ;
    fontsize(gca,12,"points") ;
    xticks(0:10:100);
    xticklabels({'0','wall start', '20','30', '40','50','60','70', '80', '90','wall end'});
    xlabel('section','FontSize',14) ;
    pbaspect([1.618 1 1]) ;
    set(gca, 'LineWidth', 2);
    axis tight;
    ylim([mimin mamax]);
    ylabel('quantitative velocity mm/s','FontSize',14) ;
    writeVideo(v, getframe(fifig));
end
close(v)
video = subVideo;



end