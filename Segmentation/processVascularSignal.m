function [mask, R_Vessel] = processVascularSignal(M0_ff_video, maskClean, maskVesselnessClean, corrParams, cmap, prefix, ToolBox)
% processVascularSignal - Processes vascular signal and computes correlation maps.
%
% Inputs:
%   M0_ff_video         - 3D video data (e.g., MRI or CT).
%   maskClean           - Clean mask for the vessel type (artery or vein).
%   maskVesselnessClean - Vesselness mask.
%   maskDiaphragm       - Diaphragm mask for filtering.
%   params              - Parameters for thresholding and processing.
%   cmap                - Colormap for visualization.
%   prefix              - Prefix for saving files (e.g., 'artery' or 'vein').
%   folder_steps        - Folder to save intermediate results.
%   ToolBox             - Toolbox structure containing paths and parameters.
%
% Outputs:
%   mask       - Final thresholded mask for the vessel type.
%   R_Vessel   - Correlation map for the vessel type.

params = ToolBox.getParams;

[numX, numY, ~] = size(M0_ff_video);
diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);

% Step 1: Compute the vascular signal
vascularSignal = sum(M0_ff_video .* maskClean, [1 2]);
vascularSignal = vascularSignal ./ nnz(maskClean);
vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);

% Step 2: Compute the centered M0

M0_ff_video_centered = M0_ff_video .* maskDiaphragm - (sum(M0_ff_video .* maskDiaphragm, [1 2]) ./ nnz(maskDiaphragm));

% Step 2: Compute the correlation map
R_Vessel = mean((M0_ff_video_centered .* vascularSignal_centered), 3) ./ ...
    (std(M0_ff_video_centered, [], 3) .* std(vascularSignal_centered, [], 3));

% Save the correlation map
saveImage(rescale(R_Vessel), ToolBox, sprintf("%s_2_1_CorrelMatrix.png", prefix), isStep = true);

% Step 3: Apply vesselness mask to the correlation map
R_Vessel = R_Vessel .* maskVesselnessClean;
R_Vessel = R_Vessel .* (R_Vessel > 0);

% Save the masked correlation map
saveImage(rescale(R_Vessel), ToolBox, sprintf("%s_2_2_R_Vessel.png", prefix), isStep = true);

% Step 4: Apply thresholding (manual or automatic)
if corrParams.threshold >= -1 && corrParams.threshold <= 1
    % Manual threshold
    mask = R_Vessel >= corrParams.threshold;
    graphThreshHistogram(R_Vessel, corrParams.threshold, mask, cmap, [prefix '_2_3']);
else
    % Automatic Otsu thresholding
    mask = autoOtsuThresholding(R_Vessel, maskVesselnessClean, corrParams.classes, [prefix '_2_3']);
    mask = mask | maskClean; % Include the initial clean mask
end

% Save the final thresholded mask
saveImage(mask, ToolBox, sprintf("%s_2_4_Thresh.png", prefix), isStep = true, cmap = cmap);
end
