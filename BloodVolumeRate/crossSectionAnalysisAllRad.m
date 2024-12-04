function [vr_avg_r, vr_std_r, area_r, mask_r, v_profiles_avg_r, v_profiles_std_r, sub_images_r] = crossSectionAnalysisAllRad(numSections, locs, widths, mask, v_RMS, name, flagBloodVelocityProfile, force_width)

% Parameters 
ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox);
flowRate_sliceHalfThickness = PW_params.flowRate_sliceHalfThickness;
numSectionMax = max(numSections);

% Arteries Cross-sections analysis
% Initialisation of the cells for arteries
vr_avg_r = zeros(numSectionMax, numFrames, numCircles);
vr_std_r = zeros(numSectionMax, numFrames, numCircles);
area_r = zeros(numSectionMax, numCircles);
mask_r = zeros(numX, numY, numCircles);
v_profiles_avg_r = zeros(numSectionMax, numCircles);
v_profiles_std_r = zeros(numSectionMax, numCircles);
sub_images_r = cell(numSectionMax, numCircles);

% Cross-Section Analysis of the arteries
parfor circleIdx = 1:numCircles
    [vr_avg, vr_std, area, ~, ~, maskCross, v_profiles_avg, v_profiles_std, subImg_cell, ~, ~] = crossSectionAnalysis2(ToolBox, locs{circleIdx}, widths{circleIdx}, mask, v_RMS, flowRate_sliceHalfThickness, name, flagBloodVelocityProfile, circleIdx, force_width, 1);

    if length(vr_avg) < 1
        continue
    end
    numSection = numSections_A(circleIdx);
    vr_avg_r(:, :, circleIdx) = reshape(vr_avg, 1, numSection, numFrames);
    vr_std_r(:, :, circleIdx) = reshape(vr_std, 1, numSection, numFrames);
    area_r(:, circleIdx) = reshape(area, 1, numSection);
    mask_r(:, :, circleIdx) = reshape(maskCross, 1, numX, numY);
    v_profiles_avg_r{:, circleIdx} = v_profiles_avg;
    v_profiles_std_r{:, circleIdx} = v_profiles_std;
    sub_images_r{:, circleIdx} = rescale(subImg_cell);

end

for circleIdx = 1:numCircles
    sub_images = sub_images_r{circleIdx};
    for sectionIdx = 1:size(sub_images, 2)
        sub_images_r{sectionIdx, circleIdx} = rescale(sub_images{sectionIdx});
    end
end

end