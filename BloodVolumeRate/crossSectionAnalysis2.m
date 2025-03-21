function [results] = crossSectionAnalysis2(ToolBox, locs, mask, v_RMS, circleName)

% Perform cross-section analysis on blood vessels.
%
% Inputs:
%   ToolBox     - Struct, contains parameters and paths.
%   locs        - Nx2 array, locations of vessel centers.
%   mask        - 2D array, mask for the region of interest.
%   v_RMS       - 3D array, velocity data over time.
%   circleName  - String, name of the circle (for saving results).
%
% Outputs:
%   results     - Struct containing analysis results.

% Initialize parameters
params = ToolBox.getParams;
numSections = size(locs, 1);
[numX, numY, numFrames] = size(v_RMS);

% Initialize the results struct with preallocated fields.
results = struct();

% Sub-images
results.subImg_cell = cell(1, numSections);

% Velocity and Flow Rate
results.v = zeros(numSections, numFrames); % Average velocity
results.dv = zeros(numSections, numFrames);
results.Q = zeros(numSections, numFrames); % Volumetric flow rate
results.dQ = zeros(numSections, numFrames);
results.v_profiles = cell(numSections, numFrames); % Velocity profiles
results.dv_profiles = cell(numSections, numFrames);

% Vessel Dimensions
results.D = zeros(numSections, 1);
results.dD = zeros(numSections, 1);
results.A = zeros(numSections, 1);
results.dA = zeros(numSections, 1);

% Masks
results.mask_sections = zeros(numX, numY, numSections);

% Compute mean velocity over time
v_RMS_mean_masked = squeeze(mean(v_RMS, 3)) .* mask;

% Define sub-image dimensions
subImgHW = round(0.01 * size(v_RMS_mean_masked, 1) * params.json.BloodVolumeRateAnalysis.ScaleFactorWidth);

% Initialize rejected masks
rejected_masks = zeros(numX, numY, 3);
crossSectionMask = zeros(numX, numY);

for n = 1:numSections
    % Define sub-image dimensions
    xRange = max(round(-subImgHW / 2) + locs(n, 1), 1):min(round(subImgHW / 2) + locs(n, 1), numX);
    yRange = max(round(-subImgHW / 2) + locs(n, 2), 1):min(round(subImgHW / 2) + locs(n, 2), numY);
    subImg = v_RMS_mean_masked(yRange, xRange);

    % Crop and rotate sub-image
    subImg = cropCircle(subImg);
    [subImg, tilt_angle] = rotateSubImage(subImg);
    results.subImg_cell{n} = rescale(subImg);

    % Update cross-section mask
    [crossSectionMask, maskCurrentSlice] = updateCrossSectionMask(crossSectionMask, mask, subImg, locs, n, tilt_angle, params);
    results.mask_sections(:, :, n) = maskCurrentSlice;

    % Compute the Vessel Cross Section
    figName = sprintf('%s%d', circleName, n);
    [D, dD, A, dA, c1, c2, rsquare] = computeVesselCrossSection(subImg, figName, ToolBox);
    results.D(n) = D;
    results.dD(n) = dD;
    results.A(n) = A;
    results.dA(n) = dA;

    % Generate figures
    saveCrossSectionFigure(subImg, D, ToolBox, figName);

    % Update rejected masks
    if rsquare < 0.6 || isnan(D) || D > mean(sum(subImg ~= 0, 2))
        rejected_masks(:, :, 1) = rejected_masks(:, :, 1) + maskCurrentSlice;
    else
        rejected_masks(:, :, 2) = rejected_masks(:, :, 2) + maskCurrentSlice;
    end

    % Compute blood volume rate and average velocity

    for t = 1:numFrames
        tmp = v_RMS(:, :, t) .* crossSectionMask;
        subFrame = tmp(yRange, xRange);
        subFrame = cropCircle(subFrame);
        subFrame = imrotate(subFrame, tilt_angle, 'bilinear', 'crop');

        v_profile = mean(subFrame, 1);
        v_cross = mean(subFrame(c1:c2, :), 2);

        % Compute average velocity
        v = mean(v_profile(c1:c2));

        % Compute standard deviation of velocity
        dv = std(v_cross);

        % Compute volumetric flow rate
        Q = v * A * 60; % microL/min

        % Uncertainty in volumetric flow rate
        if v ~= 0 && A ~= 0
            dQ = Q * sqrt((dv / v) ^ 2 + (dA / A) ^ 2 + (dA * dv / (A * v)) ^ 2);
        else
            dQ = 0; % Handle division by zero
        end

        % Handle NaN values
        if isnan(v)
            v = 0;
        end

        if isnan(Q)
            Q = 0;
        end

        if isnan(dv)
            dv = 0;
        end

        if isnan(dQ)
            dQ = 0;
        end

        % Store results
        results.v(n, t) = v;
        results.dv(n, t) = dv;
        results.Q(n, t) = Q;
        results.dQ(n, t) = dQ;
        results.v_profiles{n, t} = mean(subFrame, 1);
        results.dv_profiles{n, t} = std(subFrame, [], 1);
    end

end

results.crossSectionMask = crossSectionMask;
results.rejected_masks = rejected_masks;

close all;
end
