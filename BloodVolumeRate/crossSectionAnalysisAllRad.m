function [vr_avg_r, vr_std_r, area_r, mask_r, v_profiles_avg_r, v_profiles_std_r, sub_images_r, width_avg_r, width_std_r, vtop_avg_r, vtop_std_r] = crossSectionAnalysisAllRad(numSections, locs, widths, mask, v_RMS, name)

% Parameters 
TB = getGlobalToolBox;
params = TB.getParams;
flowRate_sliceHalfThickness = params.flowRate_sliceHalfThickness;
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
parfor circleIdx = 1:numCircles
    % Call crossSectionAnalysis2
    [vr_avg, vr_std, area, top_velocity, std_velocity, maskCross, ...
        v_profiles_avg, v_profiles_std, subImg_cell, ...
        width_avg, width_std, rejected_masks] = ...
        crossSectionAnalysis2(TB, locs{circleIdx}, widths{circleIdx}, ...
        mask, v_RMS, flowRate_sliceHalfThickness, name, circleIdx);

    % Map outputs to variables
    vr_avg_r{circleIdx} = vr_avg;
    vr_std_r{circleIdx} = vr_std;
    vtop_avg_r{circleIdx} = top_velocity;
    vtop_std_r{circleIdx} = std_velocity;
    area_r{circleIdx} = area;
    mask_r(:, :, circleIdx) = maskCross;
    rejected_mask_all_rad = rejected_masks + rejected_mask_all_rad;
    v_profiles_avg_r{circleIdx} = v_profiles_avg;
    v_profiles_std_r{circleIdx} = v_profiles_std;
    sub_images_r{circleIdx} = subImg_cell;
    width_avg_r{circleIdx} = width_avg;
    width_std_r{circleIdx} = width_std;
end

% Rescale sub-images
for circleIdx = 1:numCircles
    sub_images = sub_images_r{circleIdx};
    for sectionIdx = 1:size(sub_images, 2)
        sub_images_r{circleIdx}{sectionIdx} = rescale(sub_images{sectionIdx});
    end
end

imwrite(rejected_mask_all_rad, fullfile(TB.path_png, 'volumeRate', sprintf("%s_rejected_masks_%s.png", TB.main_foldername, name)))

end