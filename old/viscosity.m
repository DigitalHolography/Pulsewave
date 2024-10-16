function video = viscosity(mask_cell , velocity_cell, angle_list, one_cycle_dir, filename)

one_cycle_dir_png = fullfile(one_cycle_dir, 'png');
one_cycle_dir_eps = fullfile(one_cycle_dir, 'eps');
one_cycle_dir_txt = fullfile(one_cycle_dir, 'txt');
one_cycle_dir_avi = fullfile(one_cycle_dir, 'avi');
one_cycle_dir_mp4 = fullfile(one_cycle_dir, 'mp4');

numVessels = size(mask_cell, 2);
numFrames = size(velocity_cell{1}, 3);
numInterpFrames = 100;
%interpolation parameter
k = 2;

velocity_profiles = zeros(numInterpFrames, numFrames, numVessels);

for ii = 1 : numVessels
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
    xinterp_wall2wall = linspace(list(1),list(end),numInterpFrames);

    velocityProfileInterp = zeros(numInterpFrames, size(projVideo,2));
    for tt = 1:size(projVideo,2)
        velocityProfileInterp(:,tt) = interp1(x, projVideo(:,tt), xinterp_wall2wall);
    end

    % figure(1234)
    % plot(projImg)

    velocity_profiles(:, :, ii) = velocityProfileInterp;
    1;

end% ii (artery #)

average_velocity_profile = squeeze(mean(velocity_profiles, 3));

v = VideoWriter(fullfile(one_cycle_dir_avi, strcat(filename,'_velocity_profile.avi')));% avi
vMP4 = VideoWriter(fullfile(one_cycle_dir_mp4, strcat(filename,'_velocity_profile.mp4')),'MPEG-4');% mp4
open(v);
open(vMP4);
mimin = min(average_velocity_profile(:));
[mamax, idx_mamax] = max(average_velocity_profile(:));
[~,idx_syst] = ind2sub(size(average_velocity_profile),idx_mamax);
% x for normalize wall length for fiting
x = linspace(-1,1,length(average_velocity_profile(:,1)));

Vmax_list = zeros(numFrames,1);
alpha_list = zeros(numFrames,1);
beta_list = zeros(numFrames,1);
eta_list = zeros(numFrames,1);
viscosity_list = zeros(numFrames,1);

average_velocity_profile_systole = average_velocity_profile(:,idx_syst);
average_velocity_profile_diastole = average_velocity_profile(:,end);

for tt = 1 : numFrames
    tmp_velocity_profile = squeeze(average_velocity_profile(:,tt));

    % Use the defined function as an input to fit the function of viscosity
    tmp_fittype = fittype('Vmax .* (1-(1-alpha).* (abs(0.7*x).^beta))',...
    'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    'coefficients',{'Vmax','alpha','beta'});
    % tmp_fittype = fittype('Vmax .* (1-(1-0.13).* (abs(x).^beta))',...
    % 'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
    % 'coefficients',{'Vmax','beta'});
    [tmp_fit, R2_tmp_fit] = fit(x', tmp_velocity_profile, tmp_fittype, 'StartPoint', [40 0.7 2],'Lower', [10 -5 0], 'Upper', [80 3 6]);
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
    writeVideo(vMP4, getframe(fifig));


    Vmax_list(tt) = tmp_fit.Vmax;
    alpha_list(tt) = tmp_fit.alpha;
    % alpha_list(tt) = 0.13;
    beta_list(tt) = tmp_fit.beta;
    eta_list(tt) = (tmp_fit.beta + 1)/(tmp_fit.beta + tmp_fit.alpha);
    viscosity_list(tt) = -(eta_list(tt)-1.459)/0.017;

end
close(v)
close(vMP4)
video = subVideo;

% Systole/Diastole velocity profile

x_section = linspace(-0.7,0.7,length(squeeze(average_velocity_profile_systole)));
fit_velocity_profile_systole = Vmax_list(idx_syst)*(1-(1-alpha_list(idx_syst)).*abs(x_section).^beta_list(idx_syst));
fit_velocity_profile_diastole = Vmax_list(end)*(1-(1-alpha_list(end)).*abs(x_section).^beta_list(end));
figure(668)
plot(x_section,average_velocity_profile_systole,'-k',...
x_section,average_velocity_profile_diastole,'-k',...
    x_section,fit_velocity_profile_systole,'-r',...
    x_section,fit_velocity_profile_diastole,'-r', 'LineWidth',2)
title('Systole and diastole arterial velocity profile');
fontsize(gca,12,"points") ;
% xticks(x);
xticklabels({'section'});
xlabel('section','FontSize',14) ;
pbaspect([1.618 1 1]) ;
set(gca, 'LineWidth', 2);
axis tight;
ylim([0.9*mimin 1.1*mamax]);
ylabel('velocity (mm/s)','FontSize',14) ;


figure(666)
plot(viscosity_list)
pbaspect([1.618 1 1]);
xlabel('Frame','FontSize',14);
ylabel('Viscosity (cP)','FontSize',14); 
set(gca, 'LineWidth', 2);
axis tight;

% png
print('-f668','-dpng',fullfile(one_cycle_dir_png,strcat(filename,'_velocity_cross_section.png'))) ;
print('-f666','-dpng',fullfile(one_cycle_dir_png,strcat(filename,'_viscosity_in_time.png'))) ;
% eps
print('-f668','-depsc',fullfile(one_cycle_dir_eps,strcat(filename,'_velocity_cross_section.eps'))) ;
print('-f666','-depsc',fullfile(one_cycle_dir_eps,strcat(filename,'_viscosity_in_time.eps'))) ;

list_fig_close = [666, 899];
for ii=1:length(list_fig_close)
    close(list_fig_close(ii));
end

end