function video = viscosity(mask_cell , velocity_cell, angle_list, one_cycle_dir, filename)

num_vessels = size(mask_cell, 2);
num_frames = size(velocity_cell{1}, 3);
n_interp = 100;
%interpolation parameter
k = 2;

velocity_profiles = zeros(n_interp, num_frames, num_vessels);

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
    
    
    projVideo = sum(subVideo_interp, 3);
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
% x for normalize wall length for fiting
x = linspace(-1,1,length(average_velocity_profile(:,1)));

Vmax_list = zeros(num_frames,1);
alpha_list = zeros(num_frames,1);
beta_list = zeros(num_frames,1);
eta_list = zeros(num_frames,1);
viscosity_list = zeros(num_frames,1);

for tt = 1 : num_frames
    tmp_velocity_profile = squeeze(average_velocity_profile(:,tt));

    % Use the defined function as an input to fit the function of viscosity
    tmp_fittype = fittype('Vmax .* (1-(1-alpha).* (abs(x).^beta))',...
    'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    'coefficients',{'Vmax','alpha','beta'});
    % tmp_fittype = fittype('Vmax .* (1-(1-0.13).* (abs(x).^beta))',...
    % 'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    % 'coefficients',{'Vmax','beta'});
    [tmp_fit, R2_tmp_fit] = fit(x', tmp_velocity_profile, tmp_fittype, 'StartPoint', [15 0.7 2],'Lower', [5 0.3 1], 'Upper', [30 0.95 3]);
    R2_tmp_fit = R2_tmp_fit.rsquare;
    
    fifig = figure(899);
    plot(tmp_velocity_profile,'-k', 'LineWidth',2) ;

    hold on
    plot(tmp_fit(x), '-r', 'LineWidth',2);
    title('average wall-to-wall arterial velocity profile');
    legend(strcat('R² = ',string(R2_tmp_fit),' Vmax = ', string(tmp_fit.Vmax),' alpha = ', string(tmp_fit.alpha),' beta = ', string(tmp_fit.beta)));
    fontsize(gca,12,"points") ;
    xticks(x);
    xticklabels({'-1','wall start', '-0.6','-0.4', '-0.2','0','0.2','0.4', '0.6', '0.8','wall end'});
    xlabel('section','FontSize',14) ;
    pbaspect([1.618 1 1]) ;
    set(gca, 'LineWidth', 2);
    axis tight;
    ylim([mimin mamax]);
    ylabel('quantitative velocity mm/s','FontSize',14) ;
    hold off
    writeVideo(v, getframe(fifig));


    Vmax_list(tt) = tmp_fit.Vmax;
    alpha_list(tt) = tmp_fit.alpha;
    % alpha_list(tt) = 0.13;
    beta_list(tt) = tmp_fit.beta;
    eta_list(tt) = (tmp_fit.beta + 1)/(tmp_fit.beta + tmp_fit.alpha);
    viscosity_list(tt) = -(eta_list(tt)-1.459)/0.017;

end
close(v)
video = subVideo;
figure(666)
plot(viscosity_list)
pbaspect([1.618 1 1]);
xlabel('Frame','FontSize',14);
ylabel('Viscosity (cP)','FontSize',14); 
set(gca, 'LineWidth', 2);
axis tight;
end