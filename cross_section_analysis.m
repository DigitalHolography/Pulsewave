function [avg_blood_volume_rate, std_blood_volume_rate, cross_section_area, avg_blood_velocity, cross_section_mask, total_avg_blood_volume_rate, total_std_blood_volume_rate,velocity_profiles,subImg_cell] = cross_section_analysis(locs, width, mask, v_RMS, slice_half_thickness, k, ToolBox, path, type_of_vessel, flagBloodVelocityProfile,circle,force_width)
    % validate_cross_section
    %   Detailed explanation goes here FIXME

    if strcmp(type_of_vessel, 'artery')
        fig_idx_start = 70;
        name_section = 'A';
    else
        fig_idx_start = 80;
        name_section = 'V';
    end

    numSections = size(locs, 1);

    PW_params = Parameters_json(path);
    subImg_cell = cell([1 numSections]);
    subVideo_cell = cell([1 numSections]);

    [numX, numY, numFrames] = size(v_RMS);
    width_cross_section = zeros(numSections, 1);
    cross_section_area = zeros(numSections, 1);
    avg_blood_velocity = zeros(numSections, numFrames);
    avg_blood_volume_rate = zeros(numSections, numFrames);
    std_blood_velocity = zeros(numSections, numFrames);
    std_blood_volume_rate = zeros(numSections, numFrames);
    cross_section_mask = zeros(size(mask));
    mask_sections = zeros(numX, numY, numSections);
    total_avg_blood_volume_rate = zeros(numFrames, 1);
    total_std_blood_volume_rate = zeros(numFrames, 1);

    % %% VARIABLES FOR VELOCITY PROFILE VIDEO

    mkdir(ToolBox.PW_path_png, 'crossSection')
    mkdir(ToolBox.PW_path_png, 'projection')

    tilt_angle_list = zeros(1, length(locs));

    img_v_artery = squeeze(mean(v_RMS, 3)) .* mask;
    v_RMS_masked = v_RMS .* mask;
    velocity_profiles = cell([1 numSections]);
    for section_idx = 1:numSections % section_idx: vessel_number

        if width(section_idx) > 2
            subImgHW = round(width(section_idx) * PW_params.cropSection_scaleFactorWidth);
            %FIXME bords d IMG,

            xRange = max(round(-subImgHW / 2) + locs(section_idx, 2),1):min(round(subImgHW / 2) + locs(section_idx, 2),numX);
            yRange = max(round(-subImgHW / 2) + locs(section_idx, 1),1):min(round(subImgHW / 2) + locs(section_idx, 1),numY);
            subImg = img_v_artery(yRange, xRange);

            %make disk mask
            %FIXME img anamorphique
            subImg = cropCircle(subImg);

            angles = linspace(0, 180, 181);
            projx = zeros(size(subImg, 1), length(angles));
            projy = zeros(size(subImg, 2), length(angles));
            Video_subIm_rotate = zeros(size(subImg, 1), size(subImg, 2), length(angles));

            for theta = 1:length(angles)
                tmpImg = imrotate(subImg, angles(theta), 'bilinear', 'crop');
                Video_subIm_rotate(:, :, theta) = tmpImg;
                projx(:, theta) = squeeze(sum(tmpImg, 1));
                projy(:, theta) = squeeze(sum(tmpImg, 2));
            end

            figure(3001)
            imagesc(projx)

            figure(3002)
            imagesc(projy)
            % avi

            % [max_projx,tilt_idx] = max(projx(:),[],'all','linear');
            % [row,col] = ind2sub(size(squeeze(projx)),tilt_idx);
            % tilt_angle{section_idx} = col;
            %     [~,tilt_angle] = find(projx == max_projx); %x_max angle de rotation pour une coupe normale au vaisseau

            projx_bin = (projx == 0);
            list_x = squeeze(sum(projx_bin, 1));
            [~, idc] = max(list_x);
            tilt_angle_list(section_idx) = idc(1);
            subImg = imrotate(subImg, tilt_angle_list(section_idx), 'bilinear', 'crop');
            subImg_cell{section_idx} = subImg;
            subVideo = v_RMS_masked(yRange, xRange, :);

            for tt = 1:numFrames
                subVideo(:, :, tt) = imrotate(cropCircle(subVideo(:, :, tt)), tilt_angle_list(section_idx), 'bilinear', 'crop');
            end

            subVideo_cell{section_idx} = subVideo;
            section_cut = projx(:, tilt_angle_list(section_idx));

            for zz = 1:length(section_cut)

                if section_cut(zz) < 0
                    section_cut(zz) = 0;
                end

            end

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
            
            insert = '';
            if ~isempty(circle)
                insert = sprintf('_circle_%d',circle);
            end
            imwrite(f.cdata, fullfile(ToolBox.PW_path_png, 'projection', strcat(ToolBox.main_foldername,insert, ['_proj_' name_section num2str(section_idx) '.png'])));

            % Video_subIm_rotate = circshift(Video_subIm_rotate,[0 0 -tilt_angle_list(section_idx)]);
            w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, ['_' name_section num2str(section_idx) '.avi'])));
            tmp_video = mat2gray(Video_subIm_rotate);
            open(w)

            for theta = 1:length(angles)
                writeVideo(w, tmp_video(:, :, theta));
            end

            close(w);

            % [ ~, ~, tmp_0, ~] = findpeaks(section_cut,1:size(subImg,1), 'MinPeakWidth', round(PW_params.cropSection_scaleFactorSize*size(mask,1)));
            tmp = nnz(section_cut);
            % if tmp contains more than 1 element, we select the first one
            % if tmp is empty, we select 0 as width because then it's a 1-2 pixel noise peak
            if isempty(tmp)
                width_cross_section(section_idx) = 0;
                % display_width = 0;
            else
                width_cross_section(section_idx) = tmp(1);
                % display_width = tmp_0(1);
            end

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

            %bords blancs
            imwrite(f.cdata, fullfile(ToolBox.PW_path_png, 'crossSection', strcat(ToolBox.main_foldername, ['_' name_section num2str(section_idx) '.png'])));

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

        cross_section_area(section_idx) = pi * ((width_cross_section(section_idx) / 2) * (PW_params.cropSection_pixelSize / 2 ^ k)) ^ 2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    end


    %% Blood Volume Rate computation

    for section_idx = 1:numSections
        profils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)), numFrames],'single');
        for tt = 1:numFrames

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
            std_blood_velocity(section_idx, tt) = std(subFrame(subFrame ~= 0));

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
