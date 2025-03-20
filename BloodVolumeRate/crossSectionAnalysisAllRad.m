function [vr_avg_r, vr_std_r, area_r, mask_r, v_profiles_avg_r, v_profiles_std_r, sub_images_r, width_avg_r, width_std_r, vtop_avg_r, vtop_std_r] = crossSectionAnalysisAllRad(numSections, locs, mask, v_RMS, vesselName)

% Parameters
TB = getGlobalToolBox;
[numX, numY, ~] = size(v_RMS);
numCircles = size(numSections, 2);

% Arteries Cross-sections analysis
% Initialisation of the cells for arteries
vr_avg_r = cell(1, numCircles); % Average volume rate
vr_std_r = cell(1, numCircles); % Standard deviation of volume rate
vtop_avg_r = cell(1, numCircles); % Top velocity
vtop_std_r = cell(1, numCircles); % Standard deviation of velocity
area_r = cell(1, numCircles); % Cross-sectional area
mask_r = zeros(numX, numY, numCircles); % Cross-section mask
rejected_mask_all_rad = zeros(numX, numY, 3); % Cross-section mask
v_profiles_avg_r = cell(1, numCircles); % Velocity profiles
v_profiles_std_r = cell(1, numCircles); % Standard deviation of velocity profiles
sub_images_r = cell(1, numCircles); % Sub-images of vessels
width_avg_r = cell(1, numCircles); % Cross-section width
width_std_r = cell(1, numCircles); % Standard deviation of cross-section width

% Cross-Section Analysis of the arteries
parfor c_idx = 1:numCircles
    % Call crossSectionAnalysis2
    circleName = sprintf('C%d_%s', c_idx, vesselName);
    [results] = crossSectionAnalysis2(TB, locs{c_idx}, mask, v_RMS, circleName);

    % Map outputs to variables
    vr_avg_r{c_idx} = results.Q;
    vr_std_r{c_idx} = results.Q_std;
    vtop_avg_r{c_idx} = results.v;
    vtop_std_r{c_idx} = results.v_std;
    area_r{c_idx} = results.A;
    mask_r(:, :, c_idx) = results.crossSectionMask;
    rejected_mask_all_rad = results.rejected_masks + rejected_mask_all_rad;
    v_profiles_avg_r{c_idx} = results.v_profiles;
    v_profiles_std_r{c_idx} = results.v_profiles_std;
    sub_images_r{c_idx} = results.subImg_cell;
    width_avg_r{c_idx} = results.D;
    width_std_r{c_idx} = results.D_std;
end

% Rescale sub-images
for c_idx = 1:numCircles
    sub_images = sub_images_r{c_idx};
    for sectionIdx = 1:size(sub_images, 2)
        sub_images_r{c_idx}{sectionIdx} = rescale(sub_images{sectionIdx});
    end
end

imwrite(rejected_mask_all_rad, fullfile(TB.path_png, 'volumeRate', sprintf("%s_rejected_masks_%s.png", TB.main_foldername, vesselName)))

end