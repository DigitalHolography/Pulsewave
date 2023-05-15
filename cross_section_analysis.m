function [avg_blood_rate, cross_section_area, avg_blood_velocity, cross_section_mask] = cross_section_analysis(locs, width, mask, cx, cy, v_RMS, slice_half_thickness, k, one_cycle_dir, filename, vessel_type)
% validate_cross_section
%   Detailed explanation goes here FIXME

img_v_artery = squeeze(mean(v_RMS,3)).* mask;
width_cross_section = zeros(length(locs),1);
avg_blood_rate = zeros(length(locs),1);
cross_section_area = zeros(length(locs),1);
avg_blood_velocity = zeros(length(locs),1);
cross_section_mask = zeros(size(mask));

% %% VARIABLES FOR VELOCITY PROFILE VIDEO
subImg_cell = cell(1,length(locs));
subVideo_cell = cell(1,length(locs));
tilt_angle_list = zeros(1,length(locs));

for ii = 1:size(locs)
    subImgHW = round(width(ii)*3);
    %FIXME bords d IMG, 

    xRange = round(-subImgHW/2) + cx(locs(ii)) : round(subImgHW/2) + cx(locs(ii));
    yRange = round(-subImgHW/2) + cy(locs(ii)) : round(subImgHW/2) + cy(locs(ii));
    subImg = img_v_artery(yRange ,xRange);

    %make disk mask
    %FIXME img anamorphique
    subImg = cropCircle(subImg);

    subImg_cell{ii} = subImg;
    subVideo_cell{ii} = cropCircle(v_RMS(yRange, xRange, :));

    theta = linspace(0,180,181);
    projx = zeros(size(subImg,1),length(theta));
    projy = zeros(size(subImg,2),length(theta));
    for tt = 1:length(theta)
        tmpImg = imrotate(subImg,theta(tt),'bilinear','crop');
        projx(:,tt) = squeeze(sum(tmpImg,1));
        projy(:,tt) = squeeze(sum(tmpImg,2));
    end
    figure(3001)
    imagesc(projx)

    figure(3002)
    imagesc(projy)


    % [max_projx,tilt_idx] = max(projx(:),[],'all','linear');
    % [row,col] = ind2sub(size(squeeze(projx)),tilt_idx);
    % tilt_angle{ii} = col;
%     [~,tilt_angle] = find(projx == max_projx); %x_max angle de rotation pour une coupe normale au vaisseau
    
    projx_bin = (projx == 0);
    list_x = squeeze(sum(projx_bin, 1));
    [~, idc] = max(list_x);
    tilt_angle_list(ii) = idc(1);
    subImg = imrotate(subImg,tilt_angle_list(ii),'bilinear','crop');
    section_cut = projx(:,tilt_angle_list(ii));

    figure(1013)
    plot(section_cut);

%     [ ~, ~, tmp, ~] = findpeaks(section_cut,1:size(subImg,1),'MinPeakProminence',std(section_cut));
    [ ~, ~, tmp, ~] = findpeaks(section_cut,1:size(subImg,1), 'MinPeakWidth', round(0.004*size(mask,1)));
    % if tmp contains more than 1 element, we select the first one
    % if tmp is empty, we select 0 as width because then it's a 1-2 pixel noise peak
    if isempty(tmp) 
        width_cross_section(ii) = 0;
    else
        width_cross_section(ii) = tmp(1);
    end

    
    mask_slice_subImg = false(size(subImg,1),size(subImg,2));
    slice_center = round(size(subImg,1)/2);

    slice_half_thickness_tmp = min(slice_half_thickness,floor(subImgHW/2));
    mask_slice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp),:) = true;
    mask_slice_subImg = imrotate(double(mask_slice_subImg),-tilt_angle_list(ii),'bilinear','crop');
    mask_slice_subImg = mask_slice_subImg>0.01;

    %% Average the blood flow calculation over a circle before dividing by the section

    mask_current_slice = false(size(mask));
    mask_current_slice(1:size(mask_slice_subImg,1), 1:size(mask_slice_subImg,2)) = mask_slice_subImg;

    shift_x = cx(locs(ii)) - round(size(mask_slice_subImg,1)/2);
    shift_y = cy(locs(ii)) - round(size(mask_slice_subImg,2)/2);

    mask_current_slice = circshift(mask_current_slice, [shift_y shift_x]);
    mask_current_slice = mask_current_slice.*mask;
    mask_current_slice_inverse = ~mask_current_slice;
    mask_current_slice = double(mask_current_slice);
    mask_current_slice_inverse = double(mask_current_slice_inverse);
    cross_section_mask = cross_section_mask + mask_current_slice;

    %FIXME mean(v_RMS,3)
    tmp = squeeze(mean(v_RMS,3));
    tmp = tmp.*mask_current_slice;

    avg_blood_velocity(ii) = sum(tmp(:))/nnz(tmp(:));
    cross_section_area(ii) = pi*(width_cross_section(ii)/2*0.0102/2^k)^2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    avg_blood_rate(ii) = avg_blood_velocity(ii)*cross_section_area(ii); % mm^3/s

%     figure(101)
%     plot(plot_values);
%     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
%     % text(locs,pks,num2str((1:numel(pks))'))
%     text(locs,pks,string(round(avg_blood_rate,3)))
%     title("Peaks of luminosity")
%     pbaspect([1.618 1 1]);


end % ii

if strcmp(vessel_type,'artery')
    viscosity_video = viscosity(subImg_cell, subVideo_cell, tilt_angle_list, one_cycle_dir, filename);
end

list_fig_close = [3001,3002,1013];
for ii=1:length(list_fig_close)
    close(list_fig_close(ii));
end

end