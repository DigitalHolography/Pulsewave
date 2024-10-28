function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, v_RMS, M0_disp_video, ToolBox, k, path, flagBloodVelocityProfile,systolesIndexes)

PW_params = Parameters_json(path);

mkdir(ToolBox.PW_path_png, 'bloodVolumeRate')
mkdir(ToolBox.PW_path_eps, 'bloodVolumeRate')

[numX, numY, numFrames] = size(v_RMS);
[X, Y] = meshgrid(1:numY, 1:numX);

fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
M0_disp_image = rescale(mean(M0_disp_video, 3));

%% All circles testing

% for the all circles output
numCircles = PW_params.nbCircles;
maskSectionCircles = cell(1, numCircles);
delta_rad = (PW_params.velocityBigRadiusRatio - PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 / numCircles; %PW_params.radius_gap

for i = 1:numCircles
    rad_in = (PW_params.velocitySmallRadiusRatio) * (numY + numX) / 2 + (i - 1) * delta_rad; %PW_params.radius_gap) * (M + N) / 2 + (i-1) * delta_rad ;
    rad_out = rad_in + delta_rad;
    c1 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad_in;
    c2 = sqrt((X - ToolBox.x_barycentre) .^ 2 + (Y - ToolBox.y_barycentre) .^ 2) <= rad_out;
    maskSectionCircles(i) = {xor(c1, c2)};

    % save mask image
    createMaskSection(M0_disp_image, maskArtery, rad_in, rad_out, sprintf('_mask_artery_section_circle_%d.png', i), ToolBox, path);
end

close(156);

% for all circles output

SubImg_locs_artery_Circles = cell(numCircles);
SubImg_width_artery_Circles = cell(numCircles);
nb_sections_artery = zeros(1, numCircles);

for i = 1:numCircles
    maskSectionArtery = maskSectionCircles{i} .* maskArtery;

    [maskSectionArtery, n_] = bwlabel(maskSectionArtery);

    nb_sections_artery(i) = n_;
    masksSectionsArtery = zeros(numX, numY, nb_sections_artery(i));

    parfor section_idx = 1:nb_sections_artery(i)
        masksSectionsArtery(:, :, section_idx) = (maskSectionArtery == section_idx);
    end

    SubImg_locs_artery = zeros(nb_sections_artery(i), 2);
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

avg_bloodVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
std_bloodVolumeRateArteryR = zeros(numCircles, max(nb_sections_artery), numFrames, 'single');
cross_section_area_artery_r = zeros(numCircles, max(nb_sections_artery), 'single');
cross_section_mask_artery_r = zeros(numCircles, numY, numX, 'single');
stdCrossSectionWidthR = zeros(numCircles, max(nb_sections_artery), 'single');
velocity_profiles_r = cell([numCircles max(nb_sections_artery)]);
std_velocity_profiles_r = cell([numCircles max(nb_sections_artery)]);
sub_images_r = cell([numCircles max(nb_sections_artery)]);
force_width = [];
if ~isempty(PW_params.forcewidth)
    force_width = PW_params.forcewidth;
end
for i = 1:numCircles
    [avg_bloodVolumeRate_artery, std_bloodVolumeRate_artery, cross_section_area_artery, ~, ~, cross_section_mask_artery, velocity_profiles,std_velocity_profiles, subImg_cell,~,stdCrossSectionWidth] = crossSectionAnalysis(SubImg_locs_artery_Circles{i}, SubImg_width_artery_Circles{i}, maskArtery, v_RMS, PW_params.flowRate_sliceHalfThickness, k, ToolBox, path, 'artery', flagBloodVelocityProfile, i,force_width);

    if length(avg_bloodVolumeRate_artery) < 1
        continue
    end

    avg_bloodVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(avg_bloodVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
    std_bloodVolumeRateArteryR(i, 1:nb_sections_artery(i), :) = reshape(std_bloodVolumeRate_artery, 1, nb_sections_artery(i), numFrames);
    cross_section_area_artery_r(i, 1:nb_sections_artery(i)) = reshape(cross_section_area_artery, 1, nb_sections_artery(i));
    stdCrossSectionWidthR(i, 1:nb_sections_artery(i)) = reshape(stdCrossSectionWidth, 1, nb_sections_artery(i));
    cross_section_mask_artery_r(i, :, :) = reshape(cross_section_mask_artery, 1, numX, numY);

    for j = 1:nb_sections_artery(i)
        velocity_profiles_r{i, j} = velocity_profiles{j};
        std_velocity_profiles_r{i, j} = std_velocity_profiles{j};
        sub_images_r{i, j} = rescale(subImg_cell{j});
    end

end
if isempty(PW_params.forcewidth)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end
colors = lines(numCircles);

imgRGB = repmat(M0_disp_image,1,1,3);
for i =1:numCircles
    indxs = find(cross_section_mask_artery_r(i,:,:)>0);
    imgRGB(indxs) = colors(i,1);
    imgRGB(numY*numX+indxs) = colors(i,2);
    imgRGB(2*numY*numX+indxs) = colors(i,3);

    if i>1 % intersections should be drawn in white
        indxs = find(cross_section_mask_artery_r(i,:,:)>0&cross_section_mask_artery_r(i-1,:,:)>0);
        imgRGB(indxs) = 1;
        imgRGB(numY*numX+indxs) = 1;
        imgRGB(2*numY*numX+indxs) = 1;
    end
end
figure(16774)
imshow(imgRGB)
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'ateries_sections.png')))

figure(11174)
% fill with zero images the zeros parts
subimage_size = size(sub_images_r{1,1},1);
for i=1:numCircles
    for j=1:max(nb_sections_artery)
        if isempty(sub_images_r{i,j})
            sub_images_r{i,j} = zeros(subimage_size,'single');
        end
    end
end
montage(sub_images_r(1:numCircles,1:max(nb_sections_artery)),"Size",[max(nb_sections_artery),numCircles])
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'all_sections_with_increasing_radius.png')))


figure(16796)
cross_section_hist = histogram(2*sqrt(cross_section_area_artery_r(cross_section_area_artery_r~=0)/pi )*1000,50,FaceColor='k');
aa = axis;
aa(4) = aa(4)*1.14;
axis(aa);
title('Histogram of sections width (µm)');
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'histogram_of_section_width.png')))
writematrix(2*sqrt(cross_section_area_artery_r/pi )*1000,fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername,'section_widths.txt')));
writematrix(stdCrossSectionWidthR*PW_params.cropSection_pixelSize/(2^PW_params.k)*1000,fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername,'standard_deviation_section_width.txt')));

plot_bvr_full_field = figure(1676);

Color_std = [0.7 0.7 0.7];
rad = ((PW_params.velocitySmallRadiusRatio * (numX + numY) / 2 )+delta_rad/2:delta_rad : (PW_params.velocityBigRadiusRatio * (numX + numY) / 2)-delta_rad/2)'';
bvr_r = sum(avg_bloodVolumeRateArteryR,2);
std_bvr_r= sqrt(sum(std_bloodVolumeRateArteryR.^2,2)); % sqrt of the sum of variances
mean_bvr_r = squeeze(mean(bvr_r(:,:,index_start:index_end),3))';
mean_std_bvr_r = squeeze(rms(std_bvr_r(:,:,index_start:index_end),3))'; % quadratic mean
curve1 = mean_bvr_r + 0.5 * mean_std_bvr_r;
curve2 = mean_bvr_r - 0.5 * mean_std_bvr_r;
rad2 = [rad, fliplr(rad)];
inBetween = [curve1, fliplr(curve2)]';

fill(rad2, inBetween, Color_std);
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
plot(rad, mean_bvr_r, '-k', 'LineWidth', 2);
yline(mean(mean_bvr_r),'--k', 'LineWidth', 2);
legend({'','','','',sprintf('mean = %f µL/min',mean(mean_bvr_r)),'',''});

axis tight;
aa = axis;
aa(3)=-5;
aa(4)=95;
axis(aa);
hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('radius in pixels')
title("Radial variations of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'meanbloodVolumeRatexradius.png')))

plot_bvr_r_variance = figure(1677);

hold on;
for i = 1:numCircles
    plot(fullTime, squeeze(bvr_r(i,:,:)), 'LineWidth', 2);
end
axis tight;
aa = axis;
aa(3)=-5;
aa(4)=95;
axis(aa);
hold off
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Time average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'bloodVolumeRatevariancextime.png')))

plot_bvr_t = figure(1579);

mean_bvr_t = squeeze(mean(bvr_r,1))';
mean_bvr_t_value = mean(mean_bvr_t(index_start:index_end));
mean_std_bvr_t = squeeze(rms(std_bvr_r,1))';


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
yline(mean_bvr_t_value,'--k', 'LineWidth', 2)

plot(fullTime(index_start), 1.7*mean_bvr_t_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_end), 1.7*mean_bvr_t_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_start:index_end),repmat(1.7*mean_bvr_t_value,index_end-index_start+1),'-k');
legend({'','','','',sprintf('mean = %f µL/min',mean_bvr_t_value),'',''});
axis tight;
aa = axis;
aa(3)=-5;
aa(4)=95;
axis(aa);
hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_%s", ToolBox.main_foldername,'bloodVolumeRateallradxtime.png')))

if flagBloodVelocityProfile
    for i = 1:numCircles
        plot_mean_velocity_profiles = figure(7579+i);
        for j=1:nb_sections_artery(i)
            plot(mean(velocity_profiles_r{i,j},2))
            hold on
        end
        colors = lines(nb_sections_artery(i));
        for j=1:nb_sections_artery(i)
            profile = mean(velocity_profiles_r{i,j},2);
            if any(profile<0) % edge case when there is negative velocities
                [~,locs] = findpeaks(-profile);
                % we find the minimums and set them as the borders of the
                % vessel profile
                if length(locs)>1
                    indx = locs(1):locs(end);
                else
                    if locs(1)>length(profile)/2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(profile);
                    end
                end
            else % main case
                indx = find(profile>0);
            end
            plot(indx,ones([1 length(indx)])*mean(mean(velocity_profiles_r{i,j},2)),'Color',colors(j,:))
            hold on
        end

        title(['Measured time-averaged velocity profiles at radius = ',num2str(rad(i)),' pix'])
        set(gca, 'Linewidth', 2)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i,'bloodVelocityProfiles.png')))


        title(['Mean velocity profiles at radius = ',num2str(rad(i)),' pix'])
        plot_inter_velocity_profile = figure(7503+i);
        Ninterp = 50;
        interp_profile = zeros([nb_sections_artery(i),Ninterp],'single');
        interp_profile_std = zeros([nb_sections_artery(i),Ninterp],'single');
        for j=1:nb_sections_artery(i)

            profile = mean(velocity_profiles_r{i,j},2); % mean velocity profile
            profile_std = mean(std_velocity_profiles_r{i,j},2);
            if any(profile<0) % edge case when there is negative velocities
                [~,locs] = findpeaks(-profile);
                % we find the minimums and set them as the borders of the
                % vessel profile
                if length(locs)>1
                    indx = locs(1):locs(end);
                else
                    if locs(1)>length(profile)/2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(profile);
                    end
                end
            else % main case
                indx = find(profile>0);
            end
            interp_profile(j,:) = interp1(1:length(indx),profile(indx),linspace(1,length(indx),Ninterp));
            interp_profile_std(j,:) = interp1(1:length(indx),profile_std(indx),linspace(1,length(indx),Ninterp));
        end
        mean_interp_profile = mean(interp_profile,1);
        std_interp_profile = mean(interp_profile_std,1);
        curve1 = mean_interp_profile + 0.5 * std_interp_profile;
        curve2 = mean_interp_profile - 0.5 * std_interp_profile;
        ft2 = [(1:Ninterp), fliplr(1:Ninterp)];
        inBetween = [curve1, fliplr(curve2)]';

        fill(ft2, inBetween, Color_std);
        hold on;
        plot(1:Ninterp, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(1:Ninterp, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(1:Ninterp, mean_interp_profile, '-k', 'LineWidth', 2);
        axis tight;


        % adding a poiseuille fiting (poly2)
        [~,centt] = max(mean_interp_profile);
        central_range = 1:Ninterp;%max(1,centt-round(Ninterp/6)):min(Ninterp,centt+round(Ninterp/6));
        r_range = (central_range - centt);
        f = fit(r_range',mean_interp_profile(central_range)','poly2');
        poiseuille_fit = f.p1*((1:Ninterp) -centt).^2+f.p2*((1:Ninterp) -centt)+f.p3;
        poiseuille_fit(poiseuille_fit<0)=0;
        plot(poiseuille_fit, '-r', 'LineWidth', 2);

        axis tight;
        aa = axis;
        aa(3)=-5;
        aa(4)=25;
        axis(aa);
        hold off

        title(['Interpolated time-averaged velocity profile at radius = ',num2str(rad(i)),' pix'])
        set(gca, 'Linewidth', 2)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodVolumeRate', sprintf("%s_circle_%d_%s", ToolBox.main_foldername,i,'interpolatedBloodVelocityProfile.png')))

    end
end
close all
fprintf("- Blood Volume Rate for all radii took : %ds", round(toc))

end
