function [vr_avg_r, vr_std_r, area_r, mask_r, v_profiles_avg_r, v_profiles_std_r, sub_images_r, width_std_r] = crossSectionAnalysisAllRad(numSections, locs, widths, mask, v_RMS, name, flagBloodVelocityProfile, force_width)

% Parameters 
ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
flowRate_sliceHalfThickness = PW_params.flowRate_sliceHalfThickness;
[numX, numY, ~] = size(v_RMS);
numCircles = size(numSections, 2);

% Arteries Cross-sections analysis
% Initialisation of the cells for arteries
vr_avg_r = cell(1, numCircles);
vr_std_r = cell(1, numCircles);
area_r = cell(1, numCircles);
mask_r = zeros(numX, numY, numCircles);
v_profiles_avg_r = cell(1, numCircles);
v_profiles_std_r = cell(1, numCircles);
sub_images_r = cell(1, numCircles);
width_std_r = cell(1, numCircles);

% Cross-Section Analysis of the arteries
parfor circleIdx = 1:numCircles
    [vr_avg, vr_std, area, ~, ~, maskCross, v_profiles_avg, v_profiles_std, subImg_cell, ~, width_std] = crossSectionAnalysis2(ToolBox, locs{circleIdx}, widths{circleIdx}, mask, v_RMS, flowRate_sliceHalfThickness, name, flagBloodVelocityProfile, circleIdx, force_width, 1);

    if length(vr_avg) < 1
        continue
    end
    vr_avg_r{circleIdx} = vr_avg;
    vr_std_r{circleIdx} = vr_std;
    area_r{circleIdx} = area;
    mask_r(:, :, circleIdx) = maskCross;
    v_profiles_avg_r{circleIdx} = v_profiles_avg;
    v_profiles_std_r{circleIdx} = v_profiles_std;
    sub_images_r{circleIdx} = subImg_cell;
    width_std_r{circleIdx} = width_std;

end

for circleIdx = 1:numCircles
    sub_images = sub_images_r{circleIdx};
    for sectionIdx = 1:size(sub_images, 2)
        sub_images_r{circleIdx}{sectionIdx} = rescale(sub_images{sectionIdx});
    end
end

end