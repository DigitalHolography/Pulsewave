function [] = bloodSectionProfile(SubImage_cell, SubVideo_cell, type_of_vessel, circle_num)

ToolBox = getGlobalToolBox;
nb_section = size(SubImage_cell, 2);
numFrames = size(SubVideo_cell{1}, 3);
n_interp = 100;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);
%interpolation parameter
k = 2;

velocity_profiles = zeros(n_interp, numFrames, nb_section);
velocity_profiles_std = zeros(n_interp, numFrames, nb_section);

for ii = 1:nb_section
    subImg = SubImage_cell{ii};
    subVideo = SubVideo_cell{ii};
    
    for ii = 1:nb_section
        
        if nargin > 4
            subImg = SubImage_cell{circle_num, ii};
            subVideo = SubVideo_cell{circle_num, ii};
            nameFig = sprintf("%s_%s", type_of_vessel, circle_num);
        else
            subImg = SubImage_cell{ii};
            subVideo = SubVideo_cell{ii};
        end
        
        %% interpolate
        interp_size = 4 * size(subImg, 1) - 3;
        
        subVideo_interp = zeros(interp_size, interp_size, size(subVideo, 3));
        subImg_interp = interp2(subImg, k);
        
        for frame = 1:size(subVideo, 3)
            subVideo_interp(:, :, frame) = interp2(subVideo(:, :, frame), k);
        end
        
        avg_profile = squeeze(sum(subImg_interp, 1) / (size(subImg_interp, 1)));
        projVideo = squeeze(sum(subVideo_interp, 1) / size(subVideo_interp, 1));
        projVideo_std = squeeze(std(subVideo_interp, 0, 1));
        list = find(avg_profile > (0.1 * max(avg_profile, [], 'all')));
        
        if isempty(list)
            continue
        end
        
        x = 1:size(projVideo, 1);
        xinterp_wall2wall = linspace(list(1), list(end), n_interp);
        
        velocityProfileInterp = zeros(n_interp, size(projVideo, 2));
        
        for frameIdx = 1:size(projVideo, 2)
            velocityProfileInterp(:, frameIdx) = interp1(x, projVideo(:, frameIdx), xinterp_wall2wall);
        end
        
        velocity_profiles(:, :, ii) = velocityProfileInterp;
        
        velocityProfileInterp_std = zeros(n_interp, size(projVideo, 2));
        
        for frameIdx = 1:size(projVideo, 2)
            velocityProfileInterp_std(:, frameIdx) = interp1(x, projVideo_std(:, frameIdx), xinterp_wall2wall);
        end
        
        velocity_profiles_std(:, :, ii) = velocityProfileInterp_std;
        
    end % ii (artery #)
    
    average_velocity_profile = squeeze(mean(velocity_profiles, 3));
    average_velocity_profile_std = squeeze(mean(velocity_profiles_std, 3));
    mimin = min(average_velocity_profile(:));
    [mamax, idx_mamax] = max(average_velocity_profile(:));
    video_name = strcat('_velocity_profile_', nameFig);
    v = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, strcat(video_name, '.avi')))); % avi
    vMP4 = VideoWriter(fullfile(ToolBox.PW_path_mp4, strcat(ToolBox.main_foldername, strcat(video_name, '.mp4'))), 'MPEG-4'); % mp4
    open(v);
    open(vMP4);
    
    [~, idx_syst] = ind2sub(size(average_velocity_profile), idx_mamax);
    % x for normalize wall length for fiting
    x = linspace(-1, 1, length(average_velocity_profile(:, 1)));
    
    Vmax_list = zeros(numFrames, 1);
    alpha_list = zeros(numFrames, 1);
    beta_list = zeros(numFrames, 1);
    eta_list = zeros(numFrames, 1);
    viscosity_list = zeros(numFrames, 1);
    
    average_velocity_profile_systole = average_velocity_profile(:, idx_syst);
    average_velocity_profile_diastole = average_velocity_profile(:, end);
    
    %
    % curve1 = total_avg_blood_volume_rate_vein+0.5*total_std_blood_volume_rate_vein;
    % curve2 = total_avg_blood_volume_rate_vein-0.5*total_std_blood_volume_rate_vein;
    % fullTime2 = [fullTime, fliplr(fullTime)];
    % inBetween = [curve1', fliplr(curve2')];
    % fill(fullTime2, inBetween, Color_std);
    % hold on;
    % plot(fullTime,curve1,"Color",Color_std, 'LineWidth', 2);
    % plot(fullTime, curve2, "Color",Color_std, 'LineWidth', 2);
    % plot(fullTime,total_avg_blood_volume_rate_vein,'-k','LineWidth',1);
    % axis tight ;
    % hold off
    
    Color_std = [0.8 0.8 0.8];
    fullTime = 1:n_interp;
    
    if PW_params.exportVideos
        
        gifWriter = GifWriter(sprintf("velocityProfile_%s", nameFig), numFrames);
        
        for frameIdx = 1:numFrames
            tmp_velocity_profile = squeeze(average_velocity_profile(:, frameIdx));
            tmp_velocity_profile_plus_std = tmp_velocity_profile + 0.5 * average_velocity_profile_std(:, frameIdx);
            tmp_velocity_profile_minus_std = tmp_velocity_profile - 0.5 * average_velocity_profile_std(:, frameIdx);
            inBetween = [tmp_velocity_profile_plus_std', fliplr(tmp_velocity_profile_minus_std')];
            fullTime2 = [fullTime, fliplr(fullTime)];
            
            % Use the defined function as an input to fit the function of viscosity
            tmp_fittype = fittype('Vmax .* (1-(1-alpha).* (abs(0.7*x).^beta))', ...
                'dependent', {'tmp_velocity_profile'}, 'independent', {'x'}, ...
                'coefficients', {'Vmax', 'alpha', 'beta'});
            % tmp_fittype = fittype('Vmax .* (1-(1-0.13).* (abs(x).^beta))',...
            % 'dependent',{'tmp_velocity_profile'},'independent',{'x'},...
            % 'coefficients',{'Vmax','beta'});
            [tmp_fit, R2_tmp_fit] = fit(x', tmp_velocity_profile, tmp_fittype, 'StartPoint', [40 0.7 2], 'Lower', [0 -5 1.5], 'Upper', [80 3 6]);
            R2_tmp_fit = R2_tmp_fit.rsquare;
            
            fifig = figure("Visible","off");
            
            fill(fullTime2, inBetween, Color_std);
            hold on
            plot(tmp_velocity_profile, '-k', 'LineWidth', 2);
            plot(tmp_velocity_profile_plus_std, "Color", Color_std, 'LineWidth', 2)
            plot(tmp_velocity_profile_minus_std, "Color", Color_std, 'LineWidth', 2)
            plot(tmp_fit(x), '-r', 'LineWidth', 2);
            
            if strcmp(type_of_vessel, 'artery')
                title('average wall-to-wall arterial velocity profile');
            else
                title('average wall-to-wall venous velocity profile');
            end
            
            legend(strcat('R² = ', string(round(R2_tmp_fit, 2)), ' Vmax = ', string(round(tmp_fit.Vmax, 1)), ' alpha = ', string(round(tmp_fit.alpha, 1)), ' beta = ', string(round(tmp_fit.beta, 1))));
            fontsize(gca, 12, "points");
            xticks(x);
            xticklabels({'-1', 'wall start', '-0.6', '-0.4', '-0.2', '0', '0.2', '0.4', '0.6', '0.8', 'wall end'});
            xlabel('section', 'FontSize', 14);
            pbaspect([1.618 1 1]);
            set(gca, 'LineWidth', 2);
            axis tight;
            ylim([0.9 * mimin 1.1 * mamax]);
            ylabel('quantitative velocity mm/s', 'FontSize', 14);
            hold off
            drawnow
            
            frame = getframe(fifig);
            gifWriter.write(frame, frameIdx);
            
            writeVideo(v, getframe(fifig));
            writeVideo(vMP4, getframe(fifig));
            
            Vmax_list(frameIdx) = tmp_fit.Vmax;
            alpha_list(frameIdx) = tmp_fit.alpha;
            % alpha_list(tt) = 0.13;
            beta_list(frameIdx) = tmp_fit.beta;
            eta_list(frameIdx) = (tmp_fit.beta + 1) / (tmp_fit.beta + tmp_fit.alpha);
            viscosity_list(frameIdx) =- (eta_list(frameIdx) - 1.459) / 0.017;
            
        end
        
        gifWriter.generate();
        gifWriter.delete();
        
        close(v)
        close(vMP4)
        
        close 2302
        
    end
    
    % video = subVideo;
    
    % Systole/Diastole velocity profile
    average_velocity_profile_systole_plus_std = average_velocity_profile_systole + 0.5 * average_velocity_profile_std(:, idx_syst);
    average_velocity_profile_systole_minus_std = average_velocity_profile_systole - 0.5 * average_velocity_profile_std(:, idx_syst);
    average_velocity_profile_diastole_plus_std = average_velocity_profile_diastole + 0.5 * average_velocity_profile_std(:, end);
    average_velocity_profile_diastole_minus_std = average_velocity_profile_diastole - 0.5 * average_velocity_profile_std(:, end);
    inBetween_syst = [average_velocity_profile_systole_plus_std', fliplr(average_velocity_profile_systole_minus_std')];
    inBetween_diast = [average_velocity_profile_diastole_plus_std', fliplr(average_velocity_profile_diastole_minus_std')];
    fullTime2 = [fullTime, fliplr(fullTime)];
    
    x_section = linspace(-0.7, 0.7, length(squeeze(average_velocity_profile_systole)));
    fit_velocity_profile_systole = Vmax_list(idx_syst) * (1 - (1 - alpha_list(idx_syst)) .* abs(x_section) .^ beta_list(idx_syst));
    fit_velocity_profile_diastole = Vmax_list(end) * (1 - (1 - alpha_list(end)) .* abs(x_section) .^ beta_list(end));
    
    figure("Visible","off")
    % plot(x_section,average_velocity_profile_systole,'-k',...
    % x_section,average_velocity_profile_diastole,'-k',...
    %     x_section,fit_velocity_profile_systole,'-r',...
    %     x_section,fit_velocity_profile_diastole,'-r', 'LineWidth',2)
    fill(fullTime2, inBetween_syst, Color_std, 'FaceAlpha', 0.5);
    hold on
    fill(fullTime2, inBetween_diast, Color_std, 'FaceAlpha', 0.5);
    
    plot(average_velocity_profile_systole, '--k', 'LineWidth', 2);
    plot(average_velocity_profile_systole_plus_std, "Color", Color_std, 'LineWidth', 2)
    plot(average_velocity_profile_systole_minus_std, "Color", Color_std, 'LineWidth', 2)
    
    plot(average_velocity_profile_diastole, ':k', 'LineWidth', 2);
    plot(average_velocity_profile_diastole_plus_std, "Color", Color_std, 'LineWidth', 2)
    plot(average_velocity_profile_diastole_minus_std, "Color", Color_std, 'LineWidth', 2)
    
    plot(fit_velocity_profile_diastole, '-r', 'LineWidth', 2);
    plot(fit_velocity_profile_systole, '-r', 'LineWidth', 2);
    
    if strcmp(type_of_vessel, 'artery')
        title('Systole and diastole arterial velocity profile');
    else
        title('Systole and diastole venous velocity profile');
    end
    
    fontsize(gca, 12, "points");
    % xticks(x);
    xticklabels({'section'});
    xlabel('section', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;
    ylim([0.9 * mimin 1.1 * mamax]);
    ylabel('velocity (mm/s)', 'FontSize', 14);
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('viscosity_%s.png', nameFig))))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('viscosity_%s.eps', nameFig))))
    
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

    figure("Visible","off")
    plot(fullTime, viscosity_list, 'k-', 'LineWidth', 2)
    pbaspect([1.618 1 1]);
    
    if strcmp(type_of_vessel, 'artery')
        title("Viscosity in arteries")
    else
        title("Viscosity in veins")
    end
    
    ylim([-150 150])
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Viscosity (cP)', 'FontSize', 14);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('velocityCrossSection_%s.png', nameFig))))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, sprintf('velocityCrossSection_%s.png', nameFig))))

    
end
