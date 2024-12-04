function [] = bloodVolumeRateForAllRadii(maskArtery, maskVein, v_RMS, M0_ff_video, xy_barycenter, systolesIndexes, flagBloodVelocityProfile)

tic

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path);

veins_analysis = PW_params.veins_analysis;
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

mkdir(ToolBox.PW_path_png, 'crossSection')
mkdir(ToolBox.PW_path_png, 'projection')

[vr_avg_A_r, vr_std_A_r, area_A_r, mask_A_r, v_profiles_avg_A_r, v_profiles_std_A_r, sub_images_A_r] = crossSectionAnalysisAllRad(numSections_A, locs_A, widths_A, maskArtery, v_RMS, 'artery', flagBloodVelocityProfile, force_width);
if veins_analysis
    [vr_avg_V_r, vr_std_V_r, area_V_r, mask_V_r, v_profiles_avg_V_r, v_profiles_std_V_r, sub_images_V_r] = crossSectionAnalysisAllRad(numSections_V, locs_V, widths_V, maskVein, v_RMS, 'vein', flagBloodVelocityProfile, force_width);
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

sectionImage(M0_ff_img, mask_A_r, 'Artery')
if veins_analysis
    sectionImage(M0_ff_img, mask_V_r, 'Vein')
end

widthImage(sub_images_A_r, 'artery')
if veins_analysis
    widthImage(sub_images_V_r, 'vein')
end

crossSectionWidthImage(M0_ff_img, xy_barycenter, area_A_r, numSections_A, mask_A_r, locs_A, 'Artery')
if veins_analysis
    crossSectionWidthImage(M0_ff_img, xy_barycenter, area_V_r, numSections_V, mask_V_r, locs_V, 'Vein')
end

widthHistogram(area_A_r, 'artery')
if veins_analysis
    widthHistogram(area_V_r, 'vein')
end

rad = ((PW_params.velocitySmallRadiusRatio * (numX + numY) / 2) + dr / 2:dr:(PW_params.velocityBigRadiusRatio * (numX + numY) / 2) - dr / 2);

[mean_BvrT_A, mean_std_BvrT_A] = plotRadius(vr_avg_A_r, vr_std_A_r, rad, index_start, index_end, 'Artery');
if veins_analysis
    [mean_BvrT_V, mean_std_BvrT_V] = plotRadius(vr_avg_V_r, vr_std_V_r, rad, index_start, index_end, 'Vein');
end

fprintf("    4. Sections Images Generation took %ds\n", round(toc))

%% 5. Blood Flow Profiles
tic

if flagBloodVelocityProfile
    mkdir(fullfile(ToolBox.PW_path_png, 'volumeRate', 'velocityProfiles'));

    interpolatedBloodVelocityProfile(v_profiles_avg_A_r, v_profiles_std_A_r, numSections_A, rad, 50)
    if veins_analysis
        interpolatedBloodVelocityProfile(v_profiles_avg_V_r, v_profiles_std_V_r, numSections_V, rad, 50)
    end

    fprintf("    5. Profiles Images Generation took %ds\n", round(toc))
end

%% 6. Arterial Indicators
tic

graphCombined(M0_ff_video, imdilate(mask, strel('disk', PW_params.local_background_width)) & maskAllSections, [], [], mean_BvrT_A, mean_std_BvrT_A, xy_barycenter, 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in arteries Full Field', 'µL/min', skip = ~PW_params.exportVideos);
if veins_analysis
    graphCombined(M0_ff_video, imdilate(mask, strel('disk', PW_params.local_background_width)) & maskAllSections, [], [], mean_BvrT_V, mean_std_BvrT_V, xy_barycenter, 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in veins Full Field', 'µL/min', skip = ~PW_params.exportVideos);
end

ArterialResistivityIndex(fullTime, mean_BvrT_A, M0_ff_video, maskArtery)

strokeAndTotalVolume(mean_BvrT_A, mean_std_BvrT_A, systolesIndexes, fullTime, 1000)

fprintf("    6. Arterial Indicators Images Generation took %ds\n", round(toc))
%close all
fprintf("- Blood Volume Rate for all radii took : %ds\n", round(toc))

end
