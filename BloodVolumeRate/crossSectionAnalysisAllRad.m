function [Q_cell, dQ_cell, v_cell, dv_cell, v_profiles_cell, dv_profiles_cell, ...
              A_cell, D_cell, dD_cell, mask_mat, subImg_cell] = ...
    crossSectionAnalysisAllRad(numSections, locs, mask, v_RMS, vesselName)

% Parameters
ToolBox = getGlobalToolBox;
[numX, numY, ~] = size(v_RMS);
numCircles = size(numSections, 2);

% Arteries Cross-sections analysis
% Initialisation of the cells for arteries
Q_cell = cell(1, numCircles); % Average volume rate
dQ_cell = cell(1, numCircles); % Standard deviation of volume rate
v_cell = cell(1, numCircles); % Top velocity
dv_cell = cell(1, numCircles); % Standard deviation of velocity
v_profiles_cell = cell(1, numCircles); % Top velocity
dv_profiles_cell = cell(1, numCircles); % Standard deviation of velocity
A_cell = cell(1, numCircles); % Cross-sectional area
D_cell = cell(1, numCircles); % Cross-section width
dD_cell = cell(1, numCircles); % Standard deviation of cross-section width
mask_mat = zeros(numX, numY, numCircles); % Cross-section mask
rejected_mask = zeros(numX, numY, 3); % Cross-section mask
subImg_cell = cell(1, numCircles); % Sub-images of vessels

% Cross-Section Analysis of the arteries
parfor c_idx = 1:numCircles
    % Call crossSectionAnalysis2
    circleName = sprintf('C%d_%s', c_idx, vesselName);
    [results] = crossSectionAnalysis2(ToolBox, locs{c_idx}, mask, v_RMS, circleName);

    % Map outputs to variables
    v_cell{c_idx} = results.v;
    dv_cell{c_idx} = results.dv;
    v_profiles_cell{c_idx} = results.v_profiles;
    dv_profiles_cell{c_idx} = results.dv_profiles;
    Q_cell{c_idx} = results.Q;
    dQ_cell{c_idx} = results.dQ;

    A_cell{c_idx} = results.A;
    D_cell{c_idx} = results.D;
    dD_cell{c_idx} = results.dD;

    mask_mat(:, :, c_idx) = results.crossSectionMask;
    rejected_mask = results.rejected_masks + rejected_mask;
    subImg_cell{c_idx} = results.subImg_cell;
end

imwrite(rejected_mask, fullfile(ToolBox.path_png, 'volumeRate', sprintf("%s_rejected_masks_%s.png", ToolBox.main_foldername, vesselName)))

end
