function [avg_blood_rate, cross_section_area, avg_blood_velocity, cross_section_mask] = cross_section_analysis_new(locs, width, mask, v_RMS, slice_half_thickness, k,ToolBox, path)
% validate_cross_section
%   Detailed explanation goes here FIXME
PW_params = Parameters(path);

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
    if width(ii)>2
    subImgHW = round(width(ii)*PW_params.cropSection_scaleFactorWidth);
    %FIXME bords d IMG, 

    xRange = round(-subImgHW/2) + locs(ii,2) : round(subImgHW/2) + locs(ii,2);
    yRange = round(-subImgHW/2) + locs(ii,1) : round(subImgHW/2) + locs(ii,1);
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
    tmp_section = (section_cut./max(section_cut))*size(section_cut,1);


%     figure(1013)
%     plot(section_cut);






%     [ ~, ~, tmp, ~] = findpeaks(section_cut,1:size(subImg,1),'MinPeakProminence',std(section_cut));
    [ ~, ~, tmp, ~] = findpeaks(section_cut,1:size(subImg,1), 'MinPeakWidth', round(PW_params.cropSection_scaleFactorSize*size(mask,1)));
    % if tmp contains more than 1 element, we select the first one
    % if tmp is empty, we select 0 as width because then it's a 1-2 pixel noise peak
    if isempty(tmp) 
        width_cross_section(ii) = 0;
    else
        width_cross_section(ii) = tmp(1);
    end


    figure(70+ii)
xAx = linspace(0,size(section_cut,1),size(subImg,1));
    imagesc(subImg)
    colormap("gray")
    axis image 
    hold on 
    p = plot(xAx, tmp_section);
    p.LineWidth = 2;
    p.Color = 'red';
    p.LineStyle = ':';
    set(gca,'PlotBoxAspectRatio',  [1,1,1]);
    x = [round(size(subImg,1)/2)-round(width_cross_section(ii)/2) round(size(subImg,1)/2)+round(width_cross_section(ii)/2)];
    y = [round(size(subImg,1)/2) round(size(subImg,1)/2)];
    line(x,y,'Color','red','LineWidth',3)




    
    mask_slice_subImg = false(size(subImg,1),size(subImg,2));
    slice_center = round(size(subImg,1)/2);

    slice_half_thickness_tmp = min(slice_half_thickness,floor(subImgHW/2));
    mask_slice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp),:) = true;
    mask_slice_subImg = imrotate(double(mask_slice_subImg),-tilt_angle_list(ii),'bilinear','crop');
    mask_slice_subImg = mask_slice_subImg>PW_params.cropSection_maskThreshold;
    

    %% Average the blood flow calculation over a circle before dividing by the section

    mask_current_slice = false(size(mask));
    mask_current_slice(1:size(mask_slice_subImg,1), 1:size(mask_slice_subImg,2)) = mask_slice_subImg;

    shift_x = locs(ii,2) - round(size(mask_slice_subImg,1)/2);
    shift_y = locs(ii,1) - round(size(mask_slice_subImg,2)/2);

    mask_current_slice = circshift(mask_current_slice, [shift_y shift_x]);
    mask_current_slice = mask_current_slice.*mask;
    mask_current_slice_inverse = ~mask_current_slice;
    mask_current_slice = double(mask_current_slice);
    mask_current_slice_inverse = double(mask_current_slice_inverse);
    cross_section_mask = cross_section_mask + mask_current_slice;

    %FIXME mean(v_RMS,3)
    tmp = squeeze(mean(v_RMS,3));
    tmp = tmp.*mask_current_slice;

    avg_blood_velocity(ii) = sum(tmp(:))/nnz(tmp(:)); %% Atention, la moyenn est peu être fausse 
    cross_section_area(ii) = pi*((width_cross_section(ii)/2)*(PW_params.cropSection_pixelSize/2^k))^2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    avg_blood_rate(ii) = avg_blood_velocity(ii)*cross_section_area(ii)*60; % microL/min

%     figure(101)
%     plot(plot_values);
%     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
%     % text(locs,pks,num2str((1:numel(pks))'))
%     text(locs,pks,string(round(avg_blood_rate,3)))
%     title("Peaks of luminosity")
%     pbaspect([1.618 1 1]);

    end
end % ii

%viscosity_video = viscosity(subImg_cell, subVideo_cell, tilt_angle_list, ToolBox.PW_path_dir, ToolBox.main_foldername);



end