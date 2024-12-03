function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, v_RMS, M0_ff_video, xy_barycenter, systolesIndexes, flagBloodVelocityProfile)

tic

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path);

veins_analysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;
force_width = [];

if ~isempty(PW_params.forcewidth)
    force_width = PW_params.forcewidth;
end

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

mkdir(ToolBox.PW_path_png, 'volumeRate')
mkdir(ToolBox.PW_path_eps, 'volumeRate')

[numX, numY, numFrames] = size(v_RMS);
[X, Y] = meshgrid(1:numY, 1:numX);
Color_std = [0.7 0.7 0.7];

fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));

L = (numY + numX) / 2;

%% 1. Mask Sectionning for all circles

% for the all circles output
tic
numCircles = PW_params.nbCircles;
maskSectionCircles = zeros(numX, numY, numCircles);

r1 = (PW_params.velocitySmallRadiusRatio) * L;
r2 = (PW_params.velocityBigRadiusRatio) * L;
dr = (PW_params.velocityBigRadiusRatio - PW_params.velocitySmallRadiusRatio) * L / numCircles; %PW_params.radius_gap
if veins_analysis
    maskAllSections = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'mask_artery_all_sections', maskArtery, maskVein);
else
    maskAllSections = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'mask_artery_all_sections', maskArtery);
end

parfor circleIdx = 1:numCircles
    rad_in = r1 + (circleIdx - 1) * dr;
    rad_out = rad_in + dr;
    c1 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= rad_in;
    c2 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= rad_out;
    maskSectionCircles(:, :, circleIdx) = xor(c1, c2);

    % save mask image
    if veins_analysis
        createMaskSection(ToolBox, M0_ff_img, rad_in, rad_out, xy_barycenter, sprintf('mask_vessel_section_circle_%d', circleIdx), maskArtery, maskVein);
    else
        createMaskSection(ToolBox, M0_ff_img, rad_in, rad_out, xy_barycenter, sprintf('mask_artery_section_circle_%d', circleIdx), maskArtery);
    end
end

fprintf("    1. Mask Sectionning for all circles output took %ds\n", round(toc))

%% 2. Properties of the sections for all circles output

tic

locs_A = cell(numCircles, 1);
widths_A = cell(numCircles, 1);
numSections_A = zeros(1, numCircles);

if veins_analysis
    locs_V = cell(numCircles, 1);
    widths_V = cell(numCircles, 1);
    numSections_V = zeros(1, numCircles);
end

parfor circleIdx = 1:numCircles

    maskSection_A = maskSectionCircles(:, :, circleIdx) .* maskArtery;
    [labels_A, numSection_A] = bwlabel(maskSection_A);
    row_A = zeros(numSection_A, 1);
    col_A = zeros(numSection_A, 1);
    width_A = zeros(numSection_A);

    for sectionIdx = 1:numSection_A
        [row, col] = find(labels_A == sectionIdx);
        row_A(sectionIdx) = round(mean(row));
        col_A(sectionIdx) = round(mean(col));
        width_A(sectionIdx) = 0.01 * size(maskArtery, 1);
    end

    widths_A{circleIdx} = width_A;
    locs_A{circleIdx} = [row_A col_A];
    numSections_A(circleIdx) = numSection_A;

    if veins_analysis
        maskSection_V = maskSectionCircles(:, :, circleIdx) .* maskVein;
        [labels_V, numSection_V] = bwlabel(maskSection_V);
        row_V = zeros(numSection_V);
        col_V = zeros(numSection_V);
        width_V = zeros(numSection_V);

        for sectionIdx = 1:numSection_V
            [row, col] = find(labels_V == sectionIdx);
            row_V(sectionIdx) = round(mean(row));
            col_V(sectionIdx) = round(mean(col));
            width_V(sectionIdx) = 0.01 * size(maskVein, 1);
        end

        widths_V{circleIdx} = width_V;
        locs_V{circleIdx} = [row_V col_V];
        numSections_V(circleIdx) = numSection_V;

    end

end

fprintf("    2. Initialisation of the sections for all circles output took %ds\n", round(toc))

%% 3. Cross-sections analysis for all circles output

tic

% Parameters 
flowRate_sliceHalfThickness = PW_params.flowRate_sliceHalfThickness;

% Arteries Cross-sections analysis
% Initialisation of the cells for arteries
vr_avg_A_r = cell(1, numCircles);
vr_std_A_r = cell(1, numCircles);
area_A_r = cell(1, numCircles);
width_avg_A_r = cell(1, numCircles);
width_std_A_r = cell(1, numCircles);
v_profiles_avg_A_r = cell(1, numCircles);
v_profiles_std_A_r = cell(1, numCircles);
sub_images_A_r = cell(1, numCircles);

mkdir(ToolBox.PW_path_png, 'crossSection')
mkdir(ToolBox.PW_path_png, 'projection')

% Cross-Section Analysis of the arteries
parfor circleIdx = 1:numCircles
    [vr_avg_A, vr_std_A, area_A, ~, ~, width_avg_A, v_profiles_avg_A, v_profiles_std_A, subImg_cell_A, ~, width_std_A] = crossSectionAnalysis2(ToolBox, locs_A{circleIdx}, widths_A{circleIdx}, maskArtery, v_RMS, flowRate_sliceHalfThickness, 'artery', flagBloodVelocityProfile, circleIdx, force_width, 1);

    if length(vr_avg_A) < 1
        continue
    end
    numSection = numSections_A(circleIdx);
    vr_avg_A_r{circleIdx} = reshape(vr_avg_A, 1, numSection, numFrames);
    vr_std_A_r{circleIdx} = reshape(vr_std_A, 1, numSection, numFrames);
    area_A_r{circleIdx} = reshape(area_A, 1, numSection);
    width_std_A_r{circleIdx} = reshape(width_std_A, 1, numSection);
    width_avg_A_r{circleIdx} = reshape(width_avg_A, 1, numX, numY);
    v_profiles_avg_A_r{circleIdx} = v_profiles_avg_A;
    v_profiles_std_A_r{circleIdx} = v_profiles_std_A;
    sub_images_A_r{circleIdx} = rescale(subImg_cell_A);

end

% Vein Cross-sections analysis
if veins_analysis
    % Initialisation of the cells for veins
    vr_avg_V_r = cell(1, numCircles);
    vr_std_V_r = cell(1, numCircles);
    area_V_r = cell(1, numCircles);
    width_avg_V_r = cell(1, numCircles);
    width_std_V_r = cell(1, numCircles);
    v_profiles_avg_V_r = cell(1, numCircles);
    v_profiles_std_V_r = cell(1, numCircles);
    sub_images_V_r = cell(1, numCircles);

    % Cross-Section Analysis of the veins
    parfor circleIdx = 1:numCircles
        [vr_avg_V, vr_std_V, area_V, ~, ~, width_avg_V, v_profiles_avg_V, v_profiles_std_V, subImg_cell_V, ~, width_std_V] = crossSectionAnalysis2(ToolBox, locs_V{circleIdx}, widths_V{circleIdx}, maskVein, v_RMS, flowRate_sliceHalfThickness, 'vein', flagBloodVelocityProfile, circleIdx, force_width, 1);

        if length(vr_avg_V) < 1
            continue
        end
        numSection = numSections_V(circleIdx);
        vr_avg_V_r{circleIdx} = reshape(vr_avg_V, 1, numSection, numFrames);
        vr_std_V_r{circleIdx} = reshape(vr_std_V, 1, numSection, numFrames);
        area_V_r{circleIdx} = reshape(area_V, 1, numSection);
        width_std_V_r{circleIdx} = reshape(width_std_V, 1, numSection);
        width_avg_V_r{circleIdx} = reshape(width_avg_V, 1, numX, numY);
        v_profiles_avg_V_r{circleIdx} = v_profiles_avg_V;
        v_profiles_std_V_r{circleIdx} = v_profiles_std_V;
        sub_images_V_r{circleIdx} = rescale(subImg_cell_V);

    end
end

fprintf("    3. Cross-sections analysis for all circles output took %ds\n", round(toc))

%% 4. Sections Image
tic

if isempty(PW_params.forcewidth)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

sectionImage(numCircles, M0_ff_img, width_avg_A_r, 'arteries')
sectionImage(numCircles, M0_ff_img, width_avg_V_r, 'veins')

fprintf("    4. Section image took %ds\n", round(toc))

%% 5. Width Image

% Arteries
figure("Visible","off")
% fill with zero images the zeros parts
subimage_cells = sub_images_A_r{1};
subimage_size = size(subimage_cells{1}, 1);

for circleIdx = 1:numCircles

    for j = 1:max(numSections_A)

        if isempty(sub_images_A_r{circleIdx}{j})
            sub_images_A_r{circleIdx}{j} = zeros(subimage_size, 'single');
        end

    end

end

montage(cell2mat(sub_images_A_r), "Size", [max(numSections_A), numCircles])
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, '5_all_arteries_sections_with_increasing_radius.png')))

% Veins
figure("Visible","off")
% fill with zero images the zeros parts
subimage_size = size(sub_images_V_r{1, 1}, 1);

for circleIdx = 1:numCircles

    for j = 1:max(numSections_V)

        if isempty(sub_images_V_r{circleIdx, j})
            sub_images_V_r{circleIdx, j} = zeros(subimage_size, 'single');
        end

    end

end

montage(sub_images_V_r(1:numCircles, 1:max(numSections_V)), "Size", [max(numSections_V), numCircles])
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, '5_all_veins_sections_with_increasing_radius.png')))

%% 6. Cross-Section Width Image

section_width_plot = figure("Visible","off");
mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate'), 'sectionsWidth')
mkdir(fullfile(ToolBox.PW_path_eps, 'volumeRate'), 'sectionsWidth')
x_center = x_barycenter;
y_center = y_barycenter;

for circleIdx = 1:numCircles
    section_width_plot.Position = [200 200 600 600];
    crossSectionWidthArtery = 2 * sqrt(area_A_r(circleIdx, 1:numSections_A(circleIdx)) / pi) * 1000;
    etiquettes_frame_values = append(string(round(crossSectionWidthArtery, 1)), "µm");
    graphMaskTags(section_width_plot, M0_ff_img, squeeze(width_avg_A_r(circleIdx, :, :)), locs_A{circleIdx}, etiquettes_frame_values, x_center, y_center, Fontsize = 12);
    title(sprintf("%s", 'Cross section width in arteries (µm)'));
    set(gca, 'FontSize', 14)
    vesselWidthsVideo(:, :, :, circleIdx) = frame2im(getframe(section_width_plot));

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'crossSectionWidthArteryImage.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', 'sectionsWidth', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'crossSectionWidthArteryImage.eps')))
end

writeGifOnDisc(vesselWidthsVideo, 'sectionsWidth.gif', 0.1);



figure("Visible","off")
cross_section_hist = histogram(2 * sqrt(area_A_r(area_A_r ~= 0) / pi) * 1000, 50, FaceColor = 'k');
aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);
title('Histogram of sections width (µm)');
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'histogram_of_section_width.png')))
writematrix(2 * sqrt(area_A_r / pi) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'section_widths.txt')));
writematrix(width_std_A_r * PW_params.cropSection_pixelSize / (2 ^ PW_params.k) * 1000, fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'standard_deviation_section_width.txt')));

plot_Bvr_full_field = figure("Visible","off");

Color_std = [0.7 0.7 0.7];
rad = ((PW_params.velocitySmallRadiusRatio * (numX + numY) / 2) + dr / 2:dr:(PW_params.velocityBigRadiusRatio * (numX + numY) / 2) - dr / 2)'';
BvrR = sum(vr_avg_A_r, 2);
std_BvrR = sqrt(sum(vr_std_A_r .^ 2, 2)); % sqrt of the sum of variances
mean_BvrR = squeeze(mean(BvrR(:, :, index_start:index_end), 3))';
mean_std_BvrR = squeeze(rms(std_BvrR(:, :, index_start:index_end), 3))'; % quadratic mean
curve1 = mean_BvrR + 0.5 * mean_std_BvrR;
curve2 = mean_BvrR - 0.5 * mean_std_BvrR;
rad2 = [rad, fliplr(rad)];
inBetween = [curve1, fliplr(curve2)]';

fill(rad2, inBetween, Color_std);
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
plot(rad, mean_BvrR, '-k', 'LineWidth', 2);
yline(mean(mean_BvrR), '--k', 'LineWidth', 2);
legend({'', '', '', '', sprintf('mean = %0.2f µL/min', mean(mean_BvrR)), '', ''});

axis tight;
aa = axis;
aa(3) = -5;
aa(4) = 95;
axis(aa);
hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('radius in pixels')
title("Time average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'LineWidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'meanvolumeRatexradius.png')))

plot_BvrR_variance = figure("Visible","off");

hold on;

for circleIdx = 1:numCircles
    plot(fullTime, squeeze(BvrR(circleIdx, :, :)), 'LineWidth', 2);
end

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial variations of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRatevariancextime.png')))

plot_BvrT = figure("Visible","off");

mean_BvrT = squeeze(mean(BvrR, 1))';
mean_BvrT_value = mean(mean_BvrT(index_start:index_end));
max_BvrT_value = max(mean_BvrT(index_start:index_end));
mean_std_BvrT = squeeze(rms(std_BvrR, 1))';

hold off
curve1 = mean_BvrT + 0.5 * mean_std_BvrT;
curve2 = mean_BvrT - 0.5 * mean_std_BvrT;
ft2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)]';

fill(ft2, inBetween, Color_std);
hold on;
yline(0, 'k-', 'LineWidth', 2)
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, mean_BvrT, '-k', 'LineWidth', 2);
yline(mean_BvrT_value, '--k', 'LineWidth', 2)

plot(fullTime(index_start), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_end), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_start:index_end), repmat(1.07 * max_BvrT_value, index_end - index_start + 1), '-k');
legend({'', '', '', '', '', sprintf('mean = %0.2f µL/min', mean_BvrT_value), '', ''});

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4) * 1.2])
box on

hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateallradxtime.png')))

ArterialResistivityIndex(fullTime, mean_BvrT, M0_ff_video, maskArtery)

%maskNeigbors = mat2gray(mean(imread(fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png'))), 3)) > 0; % import mask neigbors
graphCombined(M0_ff_video, imdilate(maskArtery, strel('disk', PW_params.local_background_width)) & maskAllSections, [], [], mean_BvrT, mean_std_BvrT, xy_barycenter, 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in arteries Full Field', 'µL/min', skip = ~PW_params.exportVideos);

if flagBloodVelocityProfile
    mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles'));

    for circleIdx = 1:numCircles
        plot_mean_velocity_profiles = figure("Visible","off");

        for j = 1:numSections_A(circleIdx)
            plot(mean(v_profiles_avg_A_r{circleIdx, j}, 2))
            hold on
        end

        colors = lines(numSections_A(circleIdx));

        for j = 1:numSections_A(circleIdx)
            profile = mean(v_profiles_avg_A_r{circleIdx, j}, 2);

            if any(profile < 0) % edge case when there is negative velocities
                [~, locs] = findpeaks(-profile);
                % we find the minimums and set them as the borders of the
                % vessel profile
                if length(locs) > 1
                    indx = locs(1):locs(end);
                else

                    if locs(1) > length(profile) / 2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(profile);
                    end

                end

            else % main case
                indx = find(profile > 0);
            end

            plot(indx, ones([1 length(indx)]) * mean(mean(v_profiles_avg_A_r{circleIdx, j}, 2)), 'Color', colors(j, :))
            hold on
        end

        title(['Measured time-averaged velocity profiles at radius = ', num2str(rad(circleIdx)), ' pix'])
        set(gca, 'Linewidth', 2)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'bloodVelocityProfiles.png')))

        title(['Mean velocity profiles at radius = ', num2str(rad(circleIdx)), ' pix'])
        plot_inter_velocity_profile = figure("Visible","off");
        Ninterp = 50;
        interp_profile = zeros([numSections_A(circleIdx), Ninterp], 'single');
        interp_profile_std = zeros([numSections_A(circleIdx), Ninterp], 'single');

        for j = 1:numSections_A(circleIdx)

            profile = mean(v_profiles_avg_A_r{circleIdx, j}, 2); % mean velocity profile
            profile_std = mean(v_profiles_std_A_r{circleIdx, j}, 2);

            if any(profile < 0) % edge case when there is negative velocities
                [~, locs] = findpeaks(-profile);
                % we find the minimums and set them as the borders of the
                % vessel profile
                if length(locs) > 1
                    indx = locs(1):locs(end);
                else

                    if locs(1) > length(profile) / 2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(profile);
                    end

                end

            else % main case
                indx = find(profile > 0);
            end

            interp_profile(j, :) = interp1(1:length(indx), profile(indx), linspace(1, length(indx), Ninterp));
            interp_profile_std(j, :) = interp1(1:length(indx), profile_std(indx), linspace(1, length(indx), Ninterp));
        end

        mean_interp_profile = mean(interp_profile, 1);
        std_interp_profile = mean(interp_profile_std, 1);
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
        [~, centt] = max(mean_interp_profile);
        central_range = 1:Ninterp; %max(1,centt-round(Ninterp/6)):min(Ninterp,centt+round(Ninterp/6));
        r_range = (central_range - centt);
        f = fit(r_range', mean_interp_profile(central_range)', 'poly2');
        poiseuille_fit = f.p1 * ((1:Ninterp) -centt) .^ 2 + f.p2 * ((1:Ninterp) -centt) + f.p3;
        poiseuille_fit(poiseuille_fit < 0) = 0;
        plot(poiseuille_fit, '-r', 'LineWidth', 2);

        axis tight;
        aa = axis;
        aa(3) = -10;
        aa(4) = 30;
        axis(aa);
        hold off

        title(['Interpolated time-averaged velocity profile at radius = ', num2str(rad(circleIdx)), ' pix'])
        set(gca, 'Linewidth', 2)
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles', sprintf("%s_circle_%d_%s", ToolBox.main_foldername, circleIdx, 'interpolatedBloodVelocityProfile.png')))

    end

end

plot_interp_pulse = figure("Visible","off");
Ninterp = 1000;

[interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(mean_BvrT, systolesIndexes, Ninterp, mean_std_BvrT);
dt = (fullTime(2) - fullTime(1));
pulseTime = dt * (1:Ninterp) * avgLength / Ninterp;

[~, amin] = min(interp_BvrT);
[~, amax] = max(interp_BvrT);
cshiftn = Ninterp - amin;

hold off

% Retinal Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime, fliplr(pulseTime)];
inBetween = [curve1, fliplr(curve2)]';
cRose = [254, 191, 210] / 255;
fill(ft2, inBetween, cRose, 'EdgeColor', 'none');
xline(pulseTime(end), 'k--', 'LineWidth', 2)

% Remaining Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
curve1 = curve1(1:amax + cshiftn);
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime(1:amax + cshiftn), fliplr(pulseTime(1:amax + cshiftn))];
inBetween = [curve1, fliplr(curve2)]';
cCrimson = [222, 49, 99] / 255;
xline(pulseTime(amax + cshiftn), 'k--', 'LineWidth', 2)
fill(ft2, inBetween, cCrimson, 'EdgeColor', 'none');

% Grey STD and Signal
interp_BvrT2 = repmat(interp_BvrT, 1, 3);
interp_std_BvrT2 = repmat(interp_std_BvrT, 1, 3);
pulseTime2 = dt * (-Ninterp + 1:Ninterp * 2) * avgLength / Ninterp;

hold on
curve1 = circshift(interp_BvrT2, cshiftn) + 0.5 * circshift(interp_std_BvrT2, cshiftn);
curve2 = circshift(interp_BvrT2, cshiftn) - 0.5 * circshift(interp_std_BvrT2, cshiftn);
ft2 = [pulseTime2, fliplr(pulseTime2)];
inBetween = [curve1, fliplr(curve2)]';
cSTD = [0.7 0.7 0.7];
fill(ft2, inBetween, cSTD, 'EdgeColor', 'none');
plot(pulseTime2, circshift(interp_BvrT2, cshiftn), '-k', 'LineWidth', 2);
%
% axis padded
% axP = axis;
% axis([pulseTime(1)-1/2 * pulseTime(end), 3/2 * pulseTime(end), 0, axP(4)])

yline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4)])
xlim([pulseTime(1)-1/2 * pulseTime(end), 3/2 * pulseTime(end)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('Time (s)')
ccinterpBvrT = circshift(interp_BvrT, cshiftn);
dt2 = pulseTime2(2) - pulseTime2(1);
stroke_volume_value = sum(ccinterpBvrT(1:amax + cshiftn)) * dt2 / 60 * 1000; % in nL
total_volume_value = sum(ccinterpBvrT) * dt2 / 60 * 1000;
title(sprintf("Retinal Stroke Volume : %02.0f nL and Total Volume : %02.0f nL", stroke_volume_value, total_volume_value));
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
box on
set(gca, 'LineWidth', 2)
box ('on', Clipping = 'off')

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'strokeAndTotalVolume.png')))

%close all
fprintf("- Blood Volume Rate for all radii took : %ds\n", round(toc))

end
