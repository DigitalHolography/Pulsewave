function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, ~, v_RMS, dataM0, videoM0_from_holowaves, ToolBox, k, path, flagBloodVelocityProfile)

    PW_params = Parameters_json(path);

    veins_analysis = PW_params.veins_analysis;

    mkdir(ToolBox.PW_path_png, 'bloodVolumeRate')
    mkdir(ToolBox.PW_path_eps, 'bloodVolumeRate')

    [N, M] = size(maskArtery);
    N_frame = size(v_RMS, 3);
    [x, y] = meshgrid(1:M, 1:N);

    %maskArtery = imdilate(maskArtery, strel('disk', 5));
    %FIXME function velocity map

    v_RMS_AVG = mean(v_RMS, 3);
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);

    %% change mask section ?

    % radius_ratio = 0.18;
    % radius_gap = 0.05;
    % radius1 = (radius_ratio-radius_gap)* (M+N)/2;
    % radius2 = (radius_ratio+radius_gap)* (M+N)/2;

%     
    dataM0 = mat2gray(dataM0);
    maskOnes = ones(size(maskArtery, 1), size(maskArtery, 2), size(dataM0, 3));
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);
    mean_M0 = rescale(mean(videoM0_from_holowaves, 3));

    ratio_etiquette = 1.2;

    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;
    
    %% All circles testing

    % for the all circles output
    nbCircles = PW_params.nbCircles;
    maskSectionCircles = cell(1,nbCircles);
    deltr = (PW_params.velocity_bigRadiusRatio - PW_params.velocity_smallRadiusRatio) * (M + N)/2 /nbCircles; %PW_params.radius_gap 
    for i = 1:nbCircles
        rad1 = (PW_params.velocity_smallRadiusRatio) * (M + N) / 2 + (i-1) * deltr ; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * deltr ;
        rad2 = rad1 + deltr ;
        c1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad1;
        c2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= rad2;
        maskSectionCircles(i) = {xor(c1,c2)};
        
        % save mask image
        createMaskSection(mean_M0 , maskArtery,rad1,rad2,sprintf('_mask_artery_section_circle_%d.png',i), ToolBox, path);
    end
    close(156);

    % for all circles output 

    SubImg_locs_artery_Circles = cell(nbCircles);
    SubImg_width_artery_Circles = cell(nbCircles);
    nb_sections_artery = zeros(1,nbCircles);
    for i = 1:nbCircles
        maskSectionArtery = maskSectionCircles{i} .* maskArtery;

        [maskSectionArtery,n_] = bwlabel(maskSectionArtery);
    
        nb_sections_artery(i) = n_;
        masksSectionsArtery = zeros(N, M, nb_sections_artery(i));
    
        parfor section_idx = 1:nb_sections_artery(i)
            masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);
        end
        SubImg_locs_artery = zeros(nb_sections_artery(i),2);
        SubImg_width_artery = zeros(nb_sections_artery(i));
        for section_idx = 1:nb_sections_artery(i)
            [row, col] = find(masksSectionsArtery(:, :, section_idx));
            SubImg_locs_artery(section_idx, 1) = round(mean(row));
            SubImg_locs_artery(section_idx, 2) = round(mean(col));
            SubImg_width_artery(section_idx) = 0.01 * size(maskArtery, 1);
        end
        SubImg_width_artery_Circles{i} = SubImg_width_artery;
        SubImg_locs_artery_Circles{i} = SubImg_locs_artery;
    end

    % for all circles output 

    avg_bloodVolumeRate_artery_r = zeros(nbCircles,max(nb_sections_artery), N_frame,'single');
    std_bloodVolumeRate_artery_r = zeros(nbCircles,max(nb_sections_artery), N_frame,'single');
    cross_section_area_artery_r = zeros(nbCircles,max(nb_sections_artery),'single');
    cross_section_mask_artery_r = zeros(nbCircles, M,N,'single');
    velocity_profiles_r = cell([nbCircles  max(nb_sections_artery)]);
    sub_images_r = cell([nbCircles  max(nb_sections_artery)]);

    for i = 1:nbCircles
        [avg_bloodVolumeRate_artery, std_bloodVolumeRate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery, total_avg_bloodVolumeRate_artery, total_std_bloodVolumeRate_artery,velocity_profiles,subImg_cell] = cross_section_analysis(SubImg_locs_artery_Circles{i}, SubImg_width_artery_Circles{i}, maskArtery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile,i);
        if length(avg_bloodVolumeRate_artery)<1
            continue
        end
        avg_blood_velocity_artery_r(i,1:nb_sections_artery(i),:) = reshape(avg_blood_velocity_artery,1,nb_sections_artery(i),N_frame);
        avg_bloodVolumeRate_artery_r(i,1:nb_sections_artery(i),:) = reshape(avg_bloodVolumeRate_artery,1,nb_sections_artery(i),N_frame);
        std_bloodVolumeRate_artery_r(i,1:nb_sections_artery(i),:) = reshape(std_bloodVolumeRate_artery,1,nb_sections_artery(i),N_frame);
        cross_section_area_artery_r(i,1:nb_sections_artery(i)) = reshape(cross_section_area_artery,1,nb_sections_artery(i));
        cross_section_mask_artery_r(i,:,:) = reshape(cross_section_mask_artery,1,N,M);
        for j=1:nb_sections_artery(i)
            velocity_profiles_r{i,j} = velocity_profiles{j};
            sub_images_r{i,j} = rescale(subImg_cell{j});
        end

    end
    colors = lines(nbCircles);

    imgRGB = repmat(mean_M0,1,1,3);
    for i =1:nbCircles
        indxs = find(cross_section_mask_artery_r(i,:,:)>0);
        imgRGB(indxs) = colors(i,1);
        imgRGB(N*M+indxs) = colors(i,2);
        imgRGB(2*N*M+indxs) = colors(i,3);

        if i>1 % intersections should be drawn in white
            indxs = find(cross_section_mask_artery_r(i,:,:)>0&cross_section_mask_artery_r(i-1,:,:)>0);
            imgRGB(indxs) = 1;
            imgRGB(N*M+indxs) = 1;
            imgRGB(2*N*M+indxs) = 1;
        end
    end
    figure(16774)
    imshow(imgRGB)
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'ateries_sections.png')))
    
    figure(11174)
    montage(sub_images_r(1:nbCircles,1:max(nb_sections_artery)),"Size",[nbCircles,max(nb_sections_artery)])
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'all_sections_with_increasing_radius.png')))


    figure(16796)
    cross_section_hist = histogram(cross_section_area_artery_r(cross_section_area_artery_r~=0),100);
    title('histogram of cross sections');
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'histogram_of_cross_sections.png')))
    
%     discarded_cross_sections = ~(2e-3<cross_section_area_artery_r<8e-3);
%     cross_section_mask_artery_r_new = cross_section_mask_artery_r;
%     for i =1:nbCircles
%         [labels,nbsection] = bwlabel(squeeze(cross_section_mask_artery_r_new(i,:,:)));
%         for j =1:nbsection
%             if discarded_cross_sections(i,j)
%                 mask = cross_section_mask_artery_r_new(i,:,:);
%                 mask (labels==j) = 0;
%                 cross_section_mask_artery_r_new(i,:,:)=mask;
%             end
%         end
%     end
% 
%     imgRGB = repmat(mean_M0,1,1,3);
%     for i =1:nbCircles
%         indxs = find(cross_section_mask_artery_r_new(i,:,:)>0);
%         imgRGB(indxs) = colors(i,1);
%         imgRGB(N*M+indxs) = colors(i,2);
%         imgRGB(2*N*M+indxs) = colors(i,3);
% 
%         if i>1 % intersections should be drawn in white
%             indxs = find(cross_section_mask_artery_r_new(i,:,:)>0&cross_section_mask_artery_r_new(i-1,:,:)>0);
%             imgRGB(indxs) = 1;
%             imgRGB(N*M+indxs) = 1;
%             imgRGB(2*N*M+indxs) = 1;
%         end
%     end
%     figure(164)
%     imshow(imgRGB)
%     exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'ateries_sections_after_discard.png')))
% 

    plot_bvr_full_field = figure(1676);

    Color_std = [0.7 0.7 0.7];
    rad = ((PW_params.velocity_smallRadiusRatio * (M + N) / 2 )+deltr/2:deltr : (PW_params.velocity_bigRadiusRatio * (M + N) / 2)-deltr/2)'';
    bvr_r = sum(avg_bloodVolumeRate_artery_r,2);
    std_bvr_r= sum(std_bloodVolumeRate_artery_r,2);
    mean_bvr_r = squeeze(mean(bvr_r,3))';
    mean_std_bvr_r = squeeze(mean(std_bvr_r,3))';
    curve1 = mean_bvr_r + 0.5 * mean_std_bvr_r;
    curve2 = mean_bvr_r - 0.5 * mean_std_bvr_r;
    rad2 = [rad, fliplr(rad)];
    inBetween = [curve1, fliplr(curve2)]';
    
    fill(rad2, inBetween, Color_std);
    hold on;
    plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(rad, mean_bvr_r, '-k', 'LineWidth', 2);
    
    axis tight;
    aa = axis;
    aa(3)=0;
    aa(4)=1.3*aa(4);
    axis(aa);
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('radius in pixels')
    title("Mean over time Total Blood Volume Rate with the radius to center CRA")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'meanbloodVolumeRatexradius.png')))

    plot_bvr_r_variance = figure(1677);

    hold on;
    for i = 1:nbCircles
        plot(fullTime, squeeze(bvr_r(i,:,:)), 'LineWidth', 2);
    end
    axis tight;
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Total Blood Volume Rate over time in each artery sections")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'bloodVolumeRatevariancextime.png')))

    plot_bvr_t = figure(1579);

    mean_bvr_t = squeeze(mean(bvr_r,1))';
    mean_std_bvr_t = squeeze(mean(std_bvr_r,1))';

    
    hold off
    curve1 = mean_bvr_t + 0.5 * mean_std_bvr_t;
    curve2 = mean_bvr_t - 0.5 * mean_std_bvr_t;
    ft2 = [fullTime, fliplr(fullTime)];
    inBetween = [curve1, fliplr(curve2)]';
    
    fill(ft2, inBetween, Color_std);
    hold on;
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
    plot(fullTime, mean_bvr_t, '-k', 'LineWidth', 2);
    
    axis tight;
    aa = axis;
    aa(3)=0;
    aa(4)=1.3*aa(4);
    axis(aa);
    hold off
    
    ylabel('Blood Volume Rate (µL/min)')
    xlabel('time (s)')
    title("Average over all radii Total Blood Volume Rate over time")
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'bloodVolumeRateallradxtime.png')))

end
