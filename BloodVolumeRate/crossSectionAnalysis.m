function [avgVolumeRate, stdVolumeRate, crossSectionArea, avgVelocity, stdVelocity, crossSectionMask, velocityProfiles, stdVelocityProfiles, subImg_cell, crossSectionWidth, stdCrossSectionWidth] = crossSectionAnalysis(locs, width, mask, v_RMS, slice_half_thickness, type_of_vessel, flagBloodVelocityProfile, circle, force_width, flag_show_fig)
% validate_cross_section
%   Detailed explanation goes here FIXME

TB = getGlobalToolBox;
if strcmp(type_of_vessel, 'artery')
    name_section = 'A';
else
    name_section = 'V';
end

insert = '';

if ~isempty(circle)
    insert = sprintf('_circle_%d', circle);
end

numSections = size(locs, 1);

params = getParams;
k = params.k;
subImg_cell = cell([1 numSections]);
subVideo_cell = cell([1 numSections]);

[numX, numY, numFrames] = size(v_RMS);
crossSectionWidth = zeros(numSections, 1);
stdCrossSectionWidth = zeros(numSections, 1);
crossSectionArea = zeros(numSections, 1);
stdCrossSectionArea = zeros(numSections, 1);
avgVelocity = zeros(numSections, numFrames);
avgVolumeRate = zeros(numSections, numFrames);
stdVelocity = zeros(numSections, numFrames);
stdVolumeRate = zeros(numSections, numFrames);
crossSectionMask = zeros(numX, numY);
mask_sections = zeros(numX, numY, numSections);

% %% VARIABLES FOR VELOCITY PROFILE VIDEO

tilt_angle_list = zeros(1, length(locs));

v_RMS_mean_masked = squeeze(mean(v_RMS, 3)) .* mask;
v_RMS_masked = v_RMS .* mask;
velocityProfiles = cell([1 numSections]);
stdVelocityProfiles = cell([1 numSections]);

for sectionIdx = 1:numSections % sectionIdx: vessel_number
    
    if width(sectionIdx) > 2
        subImgHW = round(width(sectionIdx) * params.cropSection_scaleFactorWidth);
        %FIXME bords d IMG,
        
        xRange = max(round(-subImgHW / 2) + locs(sectionIdx, 2), 1):min(round(subImgHW / 2) + locs(sectionIdx, 2), numX);
        yRange = max(round(-subImgHW / 2) + locs(sectionIdx, 1), 1):min(round(subImgHW / 2) + locs(sectionIdx, 1), numY);
        subImg = v_RMS_mean_masked(yRange, xRange);
        
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
        
        if flag_show_fig
            figure(2200)
            imagesc(projx)
            
            figure(2201)
            imagesc(projy)
        end
        
        % avi
        
        % [max_projx,tilt_idx] = max(projx(:),[],'all','linear');
        % [row,col] = ind2sub(size(squeeze(projx)),tilt_idx);
        % tilt_angle{sectionIdx} = col;
        %     [~,tilt_angle] = find(projx == max_projx); %x_max angle de rotation pour une coupe normale au vaisseau
        
        projx_bin = (projx == 0);
        list_x = squeeze(sum(projx_bin, 1));
        [~, idc] = max(list_x);
        tilt_angle_list(sectionIdx) = idc(1);
        subImg = imrotate(subImg, tilt_angle_list(sectionIdx), 'bilinear', 'crop');
        subImg_cell{sectionIdx} = subImg;
        subVideo = v_RMS_masked(yRange, xRange, :);
        
        for tt = 1:numFrames
            subVideo(:, :, tt) = imrotate(cropCircle(subVideo(:, :, tt)), tilt_angle_list(sectionIdx), 'bilinear', 'crop');
        end
        
        subVideo_cell{sectionIdx} = subVideo;
        section_cut = projx(:, tilt_angle_list(sectionIdx));
        
        section_cut(section_cut < 0) = 0;
        
        tmp_section = (section_cut ./ max(section_cut)) * size(section_cut, 1);
        
        if flag_show_fig
            figure('Visible', 'off')
            xAx = linspace(0, size(projx, 1), size(projx, 1));
            xAy = linspace(0, size(projx, 2), size(projx, 2));
            imagesc(xAy, xAx, projx)
            colormap("gray")
            axis image
            hold on
            x = [tilt_angle_list(sectionIdx) tilt_angle_list(sectionIdx)];
            y = [0 size(projx, 1)];
            line(x, y, 'Color', 'red', 'LineStyle', ':', 'LineWidth', 3)
            axis off;
            hold off;
            set(gca, 'PlotBoxAspectRatio', [1, 1.618, 1]);
            f = getframe(gca); %# Capture the current window
            
            imwrite(f.cdata, fullfile(TB.path_png, 'volumeRate', 'projection', strcat(TB.main_foldername, insert, ['_proj_' name_section num2str(sectionIdx) '.png'])));
        end
        
        % Video_subIm_rotate = circshift(Video_subIm_rotate,[0 0 -tilt_angle_list(sectionIdx)]);
        w = VideoWriter(fullfile(TB.path_avi, strcat(TB.main_foldername, ['_' name_section num2str(sectionIdx) '.avi'])));
        tmp_video = mat2gray(Video_subIm_rotate);
        open(w)
        
        for theta = 1:length(angles)
            writeVideo(w, tmp_video(:, :, theta));
        end
        
        close(w);
        
        % [ ~, ~, tmp_0, ~] = findpeaks(section_cut,1:size(subImg,1), 'MinPeakWidth', round(params.cropSection_scaleFactorSize*size(mask,1)));
        tmp = nnz(section_cut);
        crossSectionWidth(sectionIdx) = mean(sum(subImg ~= 0, 2));
        stdCrossSectionWidth(sectionIdx) = std(sum(subImg ~= 0, 2));
        
        if flag_show_fig
            figure('Visible', 'off')
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
            x = [round(size(subImg, 1) / 2) - round(crossSectionWidth(sectionIdx) / 2) round(size(subImg, 1) / 2) + round(crossSectionWidth(sectionIdx) / 2)];
            y = [round(size(subImg, 1) / 2) round(size(subImg, 1) / 2)];
            line(x, y, 'Color', 'red', 'LineWidth', 3)
            axis off;
            f = getframe(gca); %# Capture the current
            
            %bords blancs
            imwrite(f.cdata, fullfile(TB.path_png, 'volumeRate', 'crossSection', strcat(TB.main_foldername, insert, ['_' name_section num2str(sectionIdx) '.png'])));
        end
        
        maskSlice_subImg = false(size(subImg, 1), size(subImg, 2));
        slice_center = round(size(subImg, 1) / 2);
        
        slice_half_thickness_tmp = min(slice_half_thickness, floor(subImgHW / 2));
        maskSlice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp), :) = true;
        maskSlice_subImg = imrotate(double(maskSlice_subImg), -tilt_angle_list(sectionIdx), 'bilinear', 'crop');
        maskSlice_subImg = maskSlice_subImg > params.cropSection_maskThreshold;
        
        %% Average the blood flow calculation over a circle before dividing by the section
        
        maskCurrentSlice = false(size(mask));
        maskCurrentSlice(1:size(maskSlice_subImg, 1), 1:size(maskSlice_subImg, 2)) = maskSlice_subImg;
        
        shift_x = locs(sectionIdx, 2) - round(size(maskSlice_subImg, 1) / 2);
        shift_y = locs(sectionIdx, 1) - round(size(maskSlice_subImg, 2) / 2);
        
        maskCurrentSlice = circshift(maskCurrentSlice, [shift_y shift_x]);
        maskCurrentSlice = maskCurrentSlice .* mask;
        % maskCurrentSlice_inverse = ~maskCurrentSlice;
        maskCurrentSlice = double(maskCurrentSlice);
        % maskCurrentSlice_inverse = double(maskCurrentSlice_inverse);
        crossSectionMask = crossSectionMask + maskCurrentSlice;
        mask_sections(:, :, sectionIdx) = maskCurrentSlice;
        
    end
    
    if ~isempty(force_width)
        crossSectionWidth(sectionIdx) = force_width;
    end
    
    crossSectionArea(sectionIdx) = pi * ((crossSectionWidth(sectionIdx) / 2) * (params.cropSection_pixelSize / 2 ^ k)) ^ 2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    stdCrossSectionArea(sectionIdx) = pi * (1/2 * (params.cropSection_pixelSize / 2 ^ k)) ^ 2 * sqrt(stdCrossSectionWidth(sectionIdx) ^ 4 + 2 * stdCrossSectionWidth(sectionIdx) ^ 2 * crossSectionWidth(sectionIdx) ^ 2);
end

%% Blood Volume Rate computation

for sectionIdx = 1:numSections
    stdProfils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)), numFrames], 'single');
    profils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)), numFrames], 'single');
    
    for tt = 1:numFrames
        
        current_frame = v_RMS(:, :, tt);
        
        %FIXME mean(v_RMS,3)
        tmp = current_frame .* mask_sections(:, :, sectionIdx);
        
        %tmp_velocity = zeros(1,size(nnz(tmp(:))));
        xRange = round(-subImgHW / 2) + locs(sectionIdx, 2):round(subImgHW / 2) + locs(sectionIdx, 2);
        yRange = round(-subImgHW / 2) + locs(sectionIdx, 1):round(subImgHW / 2) + locs(sectionIdx, 1);
        subFrame = tmp(yRange, xRange);
        subFrame = cropCircle(subFrame);
        subFrame = imrotate(subFrame, tilt_angle_list(sectionIdx), 'bilinear', 'crop');
        avg_profil = mean(subFrame, 1);
        profils(:, tt) = avg_profil;
        
        for ll = 1:size(subFrame, 1)
            subFrame(ll, :) = subFrame(ll, :) - avg_profil;
            
        end
        
        %FIXME calcul std avg avec des v = 0
        %avgVelocity(sectionIdx,tt) = sum(tmp(:))/nnz(tmp(:));
        avgVelocity(sectionIdx, tt) = mean(tmp(tmp ~= 0));
        
        if isnan(avgVelocity(sectionIdx, tt))
            avgVelocity(sectionIdx, tt) = 0;
        end
        
        %stdVelocity(sectionIdx,tt) = std(tmp(tmp~=0));
        stdProfils(:, tt) = std(subFrame, [], 1);
        stdVelocity(sectionIdx, tt) = mean(std(subFrame, [], 1)); % mean of std along first dimension (columns)
        
        if isnan(stdVelocity(sectionIdx, tt))
            stdVelocity(sectionIdx, tt) = 0;
        end
        
        avgVolumeRate(sectionIdx, tt) = avgVelocity(sectionIdx, tt) * crossSectionArea(sectionIdx) * 60; % microL/min
        stdVolumeRate(sectionIdx, tt) = sqrt(stdVelocity(sectionIdx, tt) ^ 2 * stdCrossSectionArea(sectionIdx) ^ 2 + stdVelocity(sectionIdx, tt) ^ 2 * crossSectionArea(sectionIdx) ^ 2 + stdCrossSectionArea(sectionIdx) ^ 2 * avgVelocity(sectionIdx, tt) ^ 2) * 60; % microL/min
        
        %     figure(101)
        %     plot(plot_values);
        %     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
        %     % text(locs,pks,num2str((1:numel(pks))'))
        %     text(locs,pks,string(round(avg_blood_rate,3)))
        %     title("Peaks of luminosity")
        %     pbaspect([1.618 1 1]);
        
        %print(['-f' num2str(70+sectionIdx)],'-dpng',fullfile(ToolBox.path_png,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(sectionIdx) '.png']))) ;
        %print(['-f' num2str(1000+sectionIdx)],'-dpng',fullfile(ToolBox.path_png,strcat(ToolBox.main_foldername,['_Proj_Artery_Section_' num2str(sectionIdx) '.png']))) ;
        
    end
    
    velocityProfiles{sectionIdx} = profils;
    stdVelocityProfiles{sectionIdx} = stdProfils;
    
    avgVolumeRate(sectionIdx, :) = filloutliers(avgVolumeRate(sectionIdx, :), 'linear');
    stdVolumeRate(sectionIdx, :) = filloutliers(stdVolumeRate(sectionIdx, :), 'linear');
    
end % sectionIdx

% propagation des incertitudes -> devrait être la somme également

if isempty(circle) && flagBloodVelocityProfile && flag_show_fig % only for the main circle (not all circles)
    
    bloodSectionProfile(subImg_cell, subVideo_cell, type_of_vessel);
    % viscosity_video = viscosity(subImg_cell, subVideo_cell, tilt_angle_list, ToolBox.path_dir, ToolBox.main_foldername);
    
end

close all
end
