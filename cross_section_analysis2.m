function [avg_blood_volume_rate, std_blood_volume_rate, cross_section_area, avg_blood_velocity, cross_section_mask, total_avg_blood_volume_rate, total_std_blood_volume_rate,velocity_profiles,std_velocity_profiles,subImg_cell] = cross_section_analysis2(locs, width, mask, v_RMS, slice_half_thickness, k, ToolBox, path, type_of_vessel, flagBloodVelocityProfile,circle,force_width)
    % validate_cross_section
    %   Detailed explanation goes here FIXME

    if strcmp(type_of_vessel, 'artery')
        fig_idx_start = 70;
        name_section = 'A';
    else
        fig_idx_start = 90;
        name_section = 'V';
    end

    insert = '';
    if ~isempty(circle)
        insert = sprintf('_circle_%d',circle);
    end

    nb_section = size(locs, 1);

    PW_params = Parameters_json(path);
    subImg_cell = cell([1 nb_section]);
    subVideo_cell = cell([1 nb_section]);

    [M, N, T_max] = size(v_RMS);
    width_cross_section = zeros(nb_section, 1);
    cross_section_area = zeros(nb_section, 1);
    avg_blood_velocity = zeros(nb_section, T_max);
    avg_blood_volume_rate = zeros(nb_section, T_max);
    std_blood_velocity = zeros(nb_section, T_max);
    std_blood_volume_rate = zeros(nb_section, T_max);
    cross_section_mask = zeros(size(mask));
    mask_sections = zeros(M, N, nb_section);
    total_avg_blood_volume_rate = zeros(T_max, 1);
    total_std_blood_volume_rate = zeros(T_max, 1);

    % %% VARIABLES FOR VELOCITY PROFILE VIDEO

    mkdir(ToolBox.PW_path_png, 'crossSection')
    mkdir(ToolBox.PW_path_png, 'projection')

    tilt_angle_list = zeros(1, length(locs));

    img_v_artery = squeeze(mean(v_RMS, 3)) .* mask;
    v_RMS_masked = v_RMS .* mask;
    velocity_profiles = cell([1 nb_section]);
    for section_idx = 1:nb_section % section_idx: vessel_number

        if width(section_idx) > 2
            subImgHW = round(width(section_idx) * PW_params.cropSection_scaleFactorWidth); % default 1% of interp image size and modified by tunable factor

            % select a sub image around the center of the section
            xRange = max(round(-subImgHW / 2) + locs(section_idx, 2),1):min(round(subImgHW / 2) + locs(section_idx, 2),N);
            yRange = max(round(-subImgHW / 2) + locs(section_idx, 1),1):min(round(subImgHW / 2) + locs(section_idx, 1),M);
            subImg = img_v_artery(yRange, xRange);

            subImg = cropCircle(subImg);

            angles = linspace(0, 180, 181);
            projx = zeros(size(subImg, 1), length(angles));
            projy = zeros(size(subImg, 2), length(angles));
            Video_subIm_rotate = zeros(size(subImg, 1), size(subImg, 2), length(angles));
            
            % calculate the sum of columns and rows
            for theta = 1:length(angles)
                tmpImg = imrotate(subImg, angles(theta), 'bilinear', 'crop');
                Video_subIm_rotate(:, :, theta) = tmpImg;
                projx(:, theta) = squeeze(sum(tmpImg, 1));
                projy(:, theta) = squeeze(sum(tmpImg, 2));
            end

            % save a video to enjoy
            w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername,insert, ['_' name_section num2str(section_idx) '.avi'])));
            tmp_video = mat2gray(Video_subIm_rotate);
            open(w)
            for theta = 1:length(angles)
                writeVideo(w, tmp_video(:, :, theta));
            end
            close(w);

            figure(3001)
            imagesc(projx)

            figure(3002)
            imagesc(projy)

            % select angle of max of sum of columns 
            projx_bin = (projx == 0);
            list_x = squeeze(sum(projx_bin, 1));
            [~, idc] = max(list_x);
            tilt_angle_list(section_idx) = idc(1);
            % rotate the sub image
            subImg = imrotate(subImg, tilt_angle_list(section_idx), 'bilinear', 'crop');
            subImg_cell{section_idx} = subImg;
            subVideo = v_RMS_masked(yRange, xRange, :);
            
            % rotate every sub frames
            for tt = 1:T_max
                subVideo(:, :, tt) = imrotate(cropCircle(subVideo(:, :, tt)), tilt_angle_list(section_idx), 'bilinear', 'crop');
            end

            subVideo_cell{section_idx} = subVideo;
            
            profile = mean(subImg,1); % mean velocity profile along the length of the section

            section_cut = projx(:, tilt_angle_list(section_idx));
            
            
            [~,centt] = max(profile);
            central_range = max(1,centt-round(subImgHW/6)):min(length(profile),centt+round(subImgHW/6));
            r_range = (central_range - centt) * PW_params.cropSection_pixelSize / 2 ^ k;
            f = fit(r_range',profile(central_range)','poly2');
            figure(section_idx)
            plot(f,r_range,profile(central_range));
            r = roots([f.p1,f.p2,f.p3]);
            width_cross_section(section_idx) = r(1)-r(2);

            set(gca, 'PlotBoxAspectRatio', [1, 1.618, 1]);
            f = getframe(gca); %# Capture the current window

            imwrite(f.cdata, fullfile(ToolBox.PW_path_png, 'projection', strcat(ToolBox.main_foldername,insert, ['_proj_' name_section num2str(section_idx) '.png'])));


            tmp_section = (section_cut ./ max(section_cut)) * size(section_cut, 1);

            figure(1000 + section_idx)
            xAx = linspace(0, size(projx, 1), size(projx, 1));
            xAy = linspace(0, size(projx, 2), size(projx, 2));
            imagesc(xAy, xAx, projx)
            colormap("gray")
            axis image
            hold on
            x = [tilt_angle_list(section_idx) tilt_angle_list(section_idx)];
            y = [0 size(projx, 1)];
            line(x, y, 'Color', 'red', 'LineStyle', ':', 'LineWidth', 3)
            axis off;
            hold off;
            set(gca, 'PlotBoxAspectRatio', [1, 1.618, 1]);
            f = getframe(gca); %# Capture the current window

            imwrite(f.cdata, fullfile(ToolBox.PW_path_png, 'projection', strcat(ToolBox.main_foldername,insert, ['_proj_' name_section num2str(section_idx) '.png'])));

            

            figure(fig_idx_start + section_idx)
            xAx = linspace(0, size(section_cut, 1), size(subImg, 1));
            imagesc(xAx, xAx, subImg)
            colormap("gray")
            axis image
            hold on
            p = plot(xAx, tmp_section);
            p.LineWidth = 2;
            p.Color = 'red';
            p.LineStyle = ':';
            set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);
            x = [round(size(subImg, 1) / 2) - round(width_cross_section(section_idx) / 2) round(size(subImg, 1) / 2) + round(width_cross_section(section_idx) / 2)];
            y = [round(size(subImg, 1) / 2) round(size(subImg, 1) / 2)];
            line(x, y, 'Color', 'red', 'LineWidth', 3)
            axis off;
            f = getframe(gca); %# Capture the current
            imwrite(f.cdata, fullfile(ToolBox.PW_path_png, 'crossSection', strcat(ToolBox.main_foldername,insert, ['_' name_section num2str(section_idx) '.png'])));

            mask_slice_subImg = false(size(subImg, 1), size(subImg, 2));
            slice_center = round(size(subImg, 1) / 2);

            slice_half_thickness_tmp = min(slice_half_thickness, floor(subImgHW / 2));
            mask_slice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp), :) = true;
            mask_slice_subImg = imrotate(double(mask_slice_subImg), -tilt_angle_list(section_idx), 'bilinear', 'crop');
            mask_slice_subImg = mask_slice_subImg > PW_params.cropSection_maskThreshold;

            %% Average the blood flow calculation over a circle before dividing by the section

            mask_current_slice = false(size(mask));
            mask_current_slice(1:size(mask_slice_subImg, 1), 1:size(mask_slice_subImg, 2)) = mask_slice_subImg;

            shift_x = locs(section_idx, 2) - round(size(mask_slice_subImg, 1) / 2);
            shift_y = locs(section_idx, 1) - round(size(mask_slice_subImg, 2) / 2);

            mask_current_slice = circshift(mask_current_slice, [shift_y shift_x]);
            mask_current_slice = mask_current_slice .* mask;
            % mask_current_slice_inverse = ~mask_current_slice;
            mask_current_slice = double(mask_current_slice);
            % mask_current_slice_inverse = double(mask_current_slice_inverse);
            cross_section_mask = cross_section_mask + mask_current_slice;
            mask_sections(:, :, section_idx) = mask_current_slice;

        end
        
        if ~isempty(force_width)
            width_cross_section(section_idx) = force_width * PW_params.cropSection_pixelSize / 2 ^ k;
        end
        cross_section_area(section_idx) = pi * ((width_cross_section(section_idx) / 2)) ^ 2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    end


    %% Blood Volume Rate computation

    for section_idx = 1:nb_section
        std_profils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)),T_max],'single');
        profils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)),T_max],'single');
        for tt = 1:T_max

            current_frame = v_RMS(:, :, tt);

            %FIXME mean(v_RMS,3)
            tmp = current_frame .* mask_sections(:, :, section_idx);

            %tmp_velocity = zeros(1,size(nnz(tmp(:))));
            xRange = round(-subImgHW / 2) + locs(section_idx, 2):round(subImgHW / 2) + locs(section_idx, 2);
            yRange = round(-subImgHW / 2) + locs(section_idx, 1):round(subImgHW / 2) + locs(section_idx, 1);
            subFrame = tmp(yRange, xRange);
            subFrame = cropCircle(subFrame);
            subFrame = imrotate(subFrame, tilt_angle_list(section_idx), 'bilinear', 'crop');
            avg_profil = mean(subFrame, 1);
            profils(:,tt) = avg_profil;
            

            for ll = 1:size(subFrame, 1)
                subFrame(ll, :) = subFrame(ll, :) - avg_profil;

            end

            %FIXME calcul std avg avec des v = 0
            %avg_blood_velocity(section_idx,tt) = sum(tmp(:))/nnz(tmp(:));
            avg_blood_velocity(section_idx, tt) = mean(tmp(tmp ~= 0));

            if isnan(avg_blood_velocity(section_idx, tt))
                avg_blood_velocity(section_idx, tt) = 0;
            end

            %std_blood_velocity(section_idx,tt) = std(tmp(tmp~=0));
            std_profils(:,tt) = std(subFrame,[],1);
            std_blood_velocity(section_idx, tt) = mean(std(subFrame,[],1)); % mean of std along first dimension (columns)

            if isnan(std_blood_velocity(section_idx, tt))
                std_blood_velocity(section_idx, tt) = 0;
            end

            avg_blood_volume_rate(section_idx, tt) = avg_blood_velocity(section_idx, tt) * cross_section_area(section_idx) * 60; % microL/min
            std_blood_volume_rate(section_idx, tt) = std_blood_velocity(section_idx, tt) * cross_section_area(section_idx) * 60; % microL/min

            %     figure(101)
            %     plot(plot_values);
            %     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
            %     % text(locs,pks,num2str((1:numel(pks))'))
            %     text(locs,pks,string(round(avg_blood_rate,3)))
            %     title("Peaks of luminosity")
            %     pbaspect([1.618 1 1]);

            %print(['-f' num2str(70+section_idx)],'-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(section_idx) '.png']))) ;
            %print(['-f' num2str(1000+section_idx)],'-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,['_Proj_Artery_Section_' num2str(section_idx) '.png']))) ;

        end

        velocity_profiles{section_idx} = profils;
        std_velocity_profiles{section_idx} = std_profils;

        avg_blood_volume_rate(section_idx, :) = filloutliers(avg_blood_volume_rate(section_idx, :), 'linear');
        std_blood_volume_rate(section_idx, :) = filloutliers(std_blood_volume_rate(section_idx, :), 'linear');

    end % section_idx

    total_avg_blood_volume_rate = sum(avg_blood_volume_rate, 1);
    total_std_blood_volume_rate = mean(std_blood_volume_rate, 1); % propagation des incertitudes -> devrait être la somme également

    if isempty(circle) && flagBloodVelocityProfile % only for the main circle (not all circles)

        bloodSectionProfile(subImg_cell, subVideo_cell, type_of_vessel, ToolBox);
        % viscosity_video = viscosity(subImg_cell, subVideo_cell, tilt_angle_list, ToolBox.PW_path_dir, ToolBox.main_foldername);

    end

    
    
    close all
end
