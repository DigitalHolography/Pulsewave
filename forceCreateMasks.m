function [maskArtery, maskVein, maskVessel, maskBackground, maskCRA, maskCRV, maskSection] = forceCreateMasks(videoM0, videoM1M0, path, ToolBox)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PW_params = Parameters_json(path);
mkdir(ToolBox.PW_path_png, 'mask')
[Nx, Ny, N_frame] = size(videoM0);

%% Manual Mask Import

maskArtery = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;
maskVein = mat2gray(mean(imread(fullfile(ToolBox.PW_path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;

%% Compute correlation to segment veins and arteries

meanIm = squeeze(mean(videoM0, 3));
meanM1M0 = squeeze(mean(videoM1M0, 3));
videoM0_zero = videoM0 - meanIm;

% compute pulse in 3 dimentions for correlation in main arteries
pulse = squeeze(mean(videoM0 .* maskArtery, [1 2]));
pulseInit = pulse - mean(pulse, "all");
pulseInit3d = zeros(Nx, Ny, N_frame);

for xx = 1:Nx

    for yy = 1:Ny
        pulseInit3d(xx, yy, :) = pulseInit;
    end
end

correlationMatrix = squeeze(mean((videoM0_zero .* pulseInit3d), 3));
correlationMatrixArtery = correlationMatrix ./ max(correlationMatrix, [], 'all');

%% Circles Sectioning

meanIm = squeeze(mean(videoM0, 3)); % Because highest intensities in CRA usually
meanIm_M1M0 = squeeze(mean(videoM1M0, 3)); % Because velocities coming from the CRA are out-of-plane
blurred_mask = imgaussfilt(double(meanIm .* meanIm_M1M0) .* correlationMatrixArtery, PW_params.gauss_filt_size_for_barycentre * size(meanIm .* meanIm_M1M0, 1), 'Padding', 0);
[ToolBox.y_barycentre, ToolBox.x_barycentre] = find(blurred_mask == max(blurred_mask, [], 'all'));

%% Create Vessel Mask

maskVessel = maskArtery | maskVein;

%% Create Background Mask

maskBackground = not(maskVessel);

%% Create Mask Section
[x, y] = meshgrid(1:Ny, 1:Nx);

radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (Ny + Nx) / 2;
radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (Ny + Nx) / 2;

circleMask1 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius1;
circleMask2 = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= radius2;

maskSection = xor(circleMask1, circleMask2);

%% Create CRA and CRV Mask

stdM1M0 = std2(meanM1M0);
maskCRA = meanM1M0 > (PW_params.CRACRV_Threshold * stdM1M0);
maskCRV = meanM1M0 < (-PW_params.CRACRV_Threshold * stdM1M0);

%% Saving masks as PNG
foldername = ToolBox.main_foldername;

imwrite(mat2gray(single(maskArtery)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskArtery_New.png')), 'png');
imwrite(mat2gray(single(maskVein)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVein_New.png')), 'png');
imwrite(mat2gray(single(maskVessel)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskVessel_New.png')), 'png');
imwrite(mat2gray(single(maskBackground)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskBackground_New.png')), 'png');
imwrite(mat2gray(single(maskSection)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSection_New.png')), 'png');
imwrite(mat2gray(single(maskCRA)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRA_New.png')), 'png');
imwrite(mat2gray(single(maskCRV)), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskCRV_New.png')), 'png');

fprintf("Manually made Masks have been used\n");
end