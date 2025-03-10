function [] = bloodVolumeRate(maskArtery, maskVein, v_RMS, M0_ff_video, xy_barycenter, systolesIndexes)

TB = getGlobalToolBox;
params = TB.getParams;

veins_analysis = params.veins_analysis;
force_width = [];

if ~isempty(params.forcewidth)
    force_width = params.forcewidth;
end

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

mkdir(TB.path_png, 'volumeRate')
mkdir(TB.path_eps, 'volumeRate')

[numX, numY, numFrames] = size(v_RMS);
[X, Y] = meshgrid(1:numY, 1:numX);

t = linspace(0, numFrames * TB.stride / TB.fs / 1000, numFrames);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));

L = (numY + numX) / 2;

%% 1. Mask Sectionning for all circles

% for the all circles output
tic
numCircles = params.nbCircles;
maskSectionCircles = zeros(numX, numY, numCircles);

r1 = (params.velocitySmallRadiusRatio) * L;
r2 = (params.velocityBigRadiusRatio) * L;
dr = (params.velocityBigRadiusRatio - params.velocitySmallRadiusRatio) * L / numCircles;
if veins_analysis
    createMaskSection(TB, M0_ff_img, r1, r2, xy_barycenter, 'mask_all_sections', maskArtery, maskVein);
else
    createMaskSection(TB, M0_ff_img, r1, r2, xy_barycenter, 'mask_artery_all_sections', maskArtery);
end

parfor circleIdx = 1:numCircles
    rad_in = r1 + (circleIdx - 1) * dr;
    rad_out = rad_in + dr;
    c1 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= rad_in;
    c2 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= rad_out;
    maskSectionCircles(:, :, circleIdx) = xor(c1, c2);

    % save mask image
    if veins_analysis
        createMaskSection(TB, M0_ff_img, rad_in, rad_out, xy_barycenter, sprintf('mask_vessel_section_circle_%d', circleIdx), maskArtery, maskVein);
    else
        createMaskSection(TB, M0_ff_img, rad_in, rad_out, xy_barycenter, sprintf('mask_artery_section_circle_%d', circleIdx), maskArtery);
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
    width_A = zeros(numSection_A, 1);

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
        row_V = zeros(numSection_V, 1);
        col_V = zeros(numSection_V, 1);
        width_V = zeros(numSection_V, 1);

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

mkdir(fullfile(TB.path_png, 'volumeRate'), 'crossSection')
mkdir(fullfile(TB.path_png, 'volumeRate'), 'projection')

[vr_avg_A_r, vr_std_A_r, area_A_r, mask_A_r, v_profiles_avg_A_r, v_profiles_std_A_r, sub_images_A_r, width_avg_A_r, width_std_A_r, vtop_avg_A_r, vtop_std_A_r] = crossSectionAnalysisAllRad(numSections_A, locs_A, widths_A, maskArtery, v_RMS, 'artery', force_width);
if veins_analysis
    [vr_avg_V_r, vr_std_V_r, area_V_r, mask_V_r, v_profiles_avg_V_r, v_profiles_std_V_r, sub_images_V_r, width_avg_V_r, width_std_V_r, vtop_avg_V_r, vtop_std_V_r] = crossSectionAnalysisAllRad(numSections_V, locs_V, widths_V, maskVein, v_RMS, 'vein', force_width);
end

[area_A_mat, ~, vr_avg_A_mat, vr_std_A_mat] = reshapeSections(numFrames, numSections_A, area_A_r, width_std_A_r, vr_avg_A_r, vr_std_A_r);
if veins_analysis
    [area_V_mat, ~, vr_avg_V_mat, vr_std_V_mat] = reshapeSections(numFrames, numSections_V, area_V_r, width_std_V_r, vr_avg_V_r, vr_std_V_r);
end

plot2csvForAllRadSection(t, vr_avg_A_r, vr_std_A_r, vr_avg_A_mat, vr_std_A_mat, 'A')
if veins_analysis
    plot2csvForAllRadSection(t, vr_avg_V_r, vr_std_V_r, vr_avg_V_mat, vr_std_V_mat, 'V')
end

topvel2csv(t, vtop_avg_A_r, vtop_std_A_r, 'A');
if veins_analysis
    topvel2csv(t, vtop_avg_V_r, vtop_std_V_r, 'V');
end

fprintf("    3. Cross-sections analysis for all circles output took %ds\n", round(toc))

%% 4. Sections Image
tic

if isempty(params.forcewidth) && ~isempty(systolesIndexes)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

if params.json.BloodVolumeRate.sectionImage
    sectionImage(M0_ff_img, mask_A_r, 'A')
    if veins_analysis
        sectionImage(M0_ff_img, mask_V_r, 'V')
    end
end

subImageSize = checkSubImgSize(sub_images_A_r);

if params.json.BloodVolumeRate.widthImage
    widthImage(subImageSize, sub_images_A_r, numSections_A, 'artery')
    if veins_analysis
        widthImage(subImageSize, sub_images_V_r, numSections_V, 'vein')
    end
end

if params.json.BloodVolumeRate.crossSectionImages
    crossSectionImages(M0_ff_img, xy_barycenter, area_A_mat, vr_avg_A_mat, v_profiles_avg_A_r,  mask_A_r, locs_A, 'Artery')
    if veins_analysis
        crossSectionImages(M0_ff_img, xy_barycenter, area_V_mat, vr_avg_V_mat, v_profiles_avg_V_r, mask_V_r, locs_V, 'Vein')
    end
end

if params.json.BloodVolumeRate.widthHistogram
    widthHistogram(width_avg_A_r, width_std_A_r,area_A_mat, 'artery');
    if veins_analysis
        widthHistogram(width_avg_V_r, width_std_V_r,area_V_mat, 'artery');
    end
end

rad = ((params.velocitySmallRadiusRatio * (numX + numY) / 2) + dr / 2:dr:(params.velocityBigRadiusRatio * (numX + numY) / 2) - dr / 2);

if params.json.BloodVolumeRate.plotRadius
    [mean_BvrT_A, mean_std_BvrT_A] = plotRadius(vr_avg_A_mat, vr_std_A_mat, t, rad, index_start, index_end, 'Artery');
    if veins_analysis
        [mean_BvrT_V, mean_std_BvrT_V] = plotRadius(vr_avg_V_mat, vr_std_V_mat, t, rad, index_start, index_end, 'Vein');
    end
end

if params.json.BloodVolumeRate.BloodFlowProfiles
    mkdir(fullfile(TB.path_png, 'volumeRate', 'velocityProfiles'));
    interpolatedBloodVelocityProfile(v_profiles_avg_A_r, v_profiles_std_A_r, numSections_A, 'A', rad, 50)
    if veins_analysis
        interpolatedBloodVelocityProfile(v_profiles_avg_V_r, v_profiles_std_V_r, numSections_V, 'V', rad, 50)
    end
end

fprintf("    4. Sections Images Generation took %ds\n", round(toc))

%% 5. Arterial Indicators
tic

% Call for arterial analysis
graphCombined(M0_ff_video, imdilate(maskArtery, strel('disk', params.local_background_width)), ...
    mean_BvrT_A, mean_std_BvrT_A, xy_barycenter, 'arterial_bvr', ...
    'etiquettes_locs', [], ...
    'etiquettes_values', [], ...
    'ylabl', 'Volume Rate (µL/min)', ...
    'xlabl', 'Time (s)', ...
    'fig_title', 'Blood Volume Rate', ...
    'unit', 'µL/min', ...
    'skip', ~params.exportVideos, ...
    'Color', 'Artery', ...
    'Visible', false);

% Call for venous analysis (if veins_analysis is true)
if veins_analysis
    graphCombined(M0_ff_video, imdilate(maskVein, strel('disk', params.local_background_width)), ...
        mean_BvrT_V, mean_std_BvrT_V, xy_barycenter, 'venous_bvr', ...
        'etiquettes_locs', [], ...
        'etiquettes_values', [], ...
        'ylabl', 'Volume Rate (µL/min)', ...
        'xlabl', 'Time (s)', ...
        'fig_title', 'Blood Volume Rate', ...
        'unit', 'µL/min', ...
        'skip', ~params.exportVideos, ...
        'Color', 'Vein', ...
        'Visible', false);
end

if params.json.BloodVolumeRate.ARIBVR
    ArterialResistivityIndex(t, mean_BvrT_A, maskArtery, 'BVR', 'volumeRate');
end

if params.json.BloodVolumeRate.strokeAndTotalVolume && ~isempty(systolesIndexes)
strokeAndTotalVolume(mean_BvrT_A, mean_std_BvrT_A, systolesIndexes, t, 1000);
end

fprintf("    5. Arterial Indicators Images Generation took %ds\n", round(toc))

end
