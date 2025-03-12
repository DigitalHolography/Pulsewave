function [avgVolumeRate, stdVolumeRate, crossSectionArea, topVelocity, stdVelocity, crossSectionMask, velocityProfiles, stdVelocityProfiles, subImg_cell, crossSectionWidth, stdCrossSectionWidth, rejected_MasksRGB] = crossSectionAnalysis2(TB, locs, width, mask, v_RMS, slice_half_thickness, type_of_vessel, circleIdx, force_width)
% Perform cross-section analysis on blood vessels.
%
% Inputs:
%   ToolBox                 - Struct, contains parameters and paths.
%   locs                    - Nx2 array, locations of vessel centers.
%   width                   - Nx1 array, widths of vessels.
%   mask                    - 2D array, mask for the region of interest.
%   v_RMS                   - 3D array, velocity data over time.
%   slice_half_thickness    - Scalar, half-thickness of the slice.
%   type_of_vessel          - String, type of vessel ('artery' or 'vein').
%   circleIdx               - Scalar, index of the circle (optional).
%   force_width             - Scalar, force a specific width (optional).
%
% Outputs:
%   avgVolumeRate           - NxF array, average volume rate over time.
%   stdVolumeRate           - NxF array, standard deviation of volume rate.
%   crossSectionArea        - Nx1 array, cross-sectional area of vessels.
%   topVelocity             - NxF array, top velocity over time.
%   stdVelocity             - NxF array, standard deviation of velocity.
%   crossSectionMask        - 2D array, mask for cross-sections.
%   velocityProfiles        - Cell array, velocity profiles for each vessel.
%   stdVelocityProfiles     - Cell array, std of velocity profiles.
%   subImg_cell             - Cell array, sub-images of vessels.
%   crossSectionWidth       - Nx1 array, width of cross-sections.
%   stdCrossSectionWidth    - Nx1 array, std of cross-section width.

% Initialize parameters
numSections = size(locs, 1);
[numX, numY, numFrames] = size(v_RMS);
params = TB.getParams;
k = params.k;
circleName = sprintf('circle_%d', circleIdx);

avgVolumeRate = zeros(numSections, numFrames);
stdVolumeRate = zeros(numSections, numFrames);
crossSectionArea = zeros(numSections, 1);
stdCrossSectionArea = zeros(numSections, 1);
avgVelocity = zeros(numSections, numFrames);
topVelocity = zeros(numSections, numFrames);
stdVelocity = zeros(numSections, numFrames);
crossSectionMask = zeros(numX, numY);
velocityProfiles = cell([1 numSections]);
stdVelocityProfiles = cell([1 numSections]);
subImg_cell = cell([1 numSections]);
crossSectionWidth = zeros(numSections, 1);
stdCrossSectionWidth = zeros(numSections, 1);
mask_sections = zeros(numX, numY, numSections);
tilt_angle_list = zeros(1, length(locs));
rejected_MasksRGB = zeros(numX, numY, 3);

v_RMS_mean_masked = squeeze(mean(v_RMS, 3)) .* mask;

for sectionIdx = 1:numSections % sectionIdx: vessel_number

    if strcmp(type_of_vessel, 'artery')
        name_section = sprintf('A%d', sectionIdx);
    else
        name_section = sprintf('V%d', sectionIdx);
    end

    if width(sectionIdx) > 2

        subImgHW = round(width(sectionIdx) * params.cropSection_scaleFactorWidth);
        %FIXME bords d IMG,

        xRange = max(round(-subImgHW / 2) + locs(sectionIdx, 2), 1):min(round(subImgHW / 2) + locs(sectionIdx, 2), numX);
        yRange = max(round(-subImgHW / 2) + locs(sectionIdx, 1), 1):min(round(subImgHW / 2) + locs(sectionIdx, 1), numY);
        subImg = v_RMS_mean_masked(yRange, xRange);

        %make disk mask
        %FIXME img anamorphique
        subImg = cropCircle(subImg);

        % Rotate the sub-image to align the blood vessel vertically
        [subImg, tilt_angle] = rotateSubImage(subImg);

        subImg_cell{sectionIdx} = subImg;
        tilt_angle_list(sectionIdx) = tilt_angle;

        profile = mean(subImg, 1);
        if mean(profile) > 0
            central_range = find(profile > 0.1 * max(profile));
        else % case of fully negative vessel taken in the choroid; this makes the vessel detection independent of sign
            central_range = find(profile < 0.1 * min(profile));
        end
        centt = mean(central_range);
        r_range = (central_range - centt) * params.cropSection_pixelSize / 2 ^ k;
        [p1, p2, p3, rsquare, p1_err, p2_err, p3_err] = customPoly2Fit(r_range', profile(central_range)');
        [r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err);

        if rsquare < 0.6 % if bad fit : reject the fit
            crossSectionWidth(sectionIdx) = mean(sum(subImg ~= 0, 2));
            stdCrossSectionWidth(sectionIdx) = std(sum(subImg ~= 0, 2));
        else
            crossSectionWidth(sectionIdx) = abs(r1 - r2) / (params.cropSection_pixelSize / 2 ^ k);
            stdCrossSectionWidth(sectionIdx) = sqrt(r1_err ^ 2 + r2_err ^ 2) / (params.cropSection_pixelSize / 2 ^ k);
        end

        if isnan(crossSectionWidth(sectionIdx)) || crossSectionWidth(sectionIdx) > mean(sum(subImg ~= 0, 2))
            crossSectionWidth(sectionIdx) = mean(sum(subImg ~= 0, 2));
            stdCrossSectionWidth(sectionIdx) = std(sum(subImg ~= 0, 2));
        end
        
        crossSectionMask = updateCrossSectionMask(crossSectionMask, mask, subImg, locs, sectionIdx, tilt_angle, slice_half_thickness, params);
        mask_sections(:, :, sectionIdx) = updateCrossSectionMask(crossSectionMask, mask, subImg, locs, sectionIdx, tilt_angle, slice_half_thickness, params);

        % Create Figures
        poiseuilleProfileFigure(subImg, profile, centt, central_range, p1, p2, p3, r1, r2, rsquare, circleName, name_section, TB)
        saveCrossSectionFigure(subImg, crossSectionWidth(sectionIdx), TB, circleName, name_section)

        if rsquare < 0.6 || isnan(crossSectionWidth(sectionIdx)) || crossSectionWidth(sectionIdx) > mean(sum(subImg ~= 0, 2))
            rejected_MasksRGB(:, :, 1) = rejected_MasksRGB(:, :, 1) + mask_sections(:, :, sectionIdx);
        else
            rejected_MasksRGB(:, :, 2) = rejected_MasksRGB(:, :, 2) + mask_sections(:, :, sectionIdx);
        end
        
    end

    if ~isempty(force_width)
        crossSectionWidth(sectionIdx) = force_width;
    end

    crossSectionArea(sectionIdx) = pi * ((crossSectionWidth(sectionIdx) * (params.cropSection_pixelSize / 2 ^ k) / 2)) ^ 2; % /2 because radius=d/2 - 0.0102/2^k mm = size pixel with k coef interpolation
    stdCrossSectionArea(sectionIdx) = pi * (1/2 * (params.cropSection_pixelSize / 2 ^ k)) ^ 2 * sqrt(stdCrossSectionWidth(sectionIdx) ^ 4 + 2 * stdCrossSectionWidth(sectionIdx) ^ 2 * crossSectionWidth(sectionIdx) ^ 2);
end

%% Blood Volume Rate computation

for sectionIdx = 1:numSections
    stdProfils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)), numFrames], 'single');
    profils = zeros([length(round(-subImgHW / 2):round(subImgHW / 2)), numFrames], 'single');

    for tt = 1:numFrames

        current_frame = v_RMS(:, :, tt);

        %FIXME mean(v_RMS,3)
        tmp = current_frame .* mask_sections(:, :, sectionIdx);

        %tmp_velocity = zeros(1,size(nnz(tmp(:))));
        xRange = round(-subImgHW / 2) + locs(sectionIdx, 2):round(subImgHW / 2) + locs(sectionIdx, 2);
        yRange = round(-subImgHW / 2) + locs(sectionIdx, 1):round(subImgHW / 2) + locs(sectionIdx, 1);
        subFrame = tmp(yRange, xRange);
        subFrame = cropCircle(subFrame);
        subFrame = imrotate(subFrame, tilt_angle_list(sectionIdx), 'bilinear', 'crop');
        avg_profil = mean(subFrame, 1);
        profils(:, tt) = avg_profil;

        % for ll = 1:size(subFrame, 1)
        %     subFrame(ll, :) = subFrame(ll, :) - avg_profil;
        % 
        % end

        %FIXME calcul std avg avec des v = 0
        %avgVelocity(sectionIdx,tt) = sum(tmp(:))/nnz(tmp(:));
        avgVelocity(sectionIdx, tt) = mean(tmp(tmp ~= 0));
        topVelocity(sectionIdx, tt) = mean(max(subFrame, 2)) + mean(min(subFrame, 2));

        if isnan(avgVelocity(sectionIdx, tt))
            avgVelocity(sectionIdx, tt) = 0;
        end
        if isnan(topVelocity(sectionIdx, tt))
            topVelocity(sectionIdx, tt) = 0;
        end

        %stdVelocity(sectionIdx,tt) = std(tmp(tmp~=0));
        stdProfils(:, tt) = std(subFrame, [], 1);


        stdVelocity(sectionIdx, tt) = std(max(subFrame, [], 2)); % mean of std along first dimension (columns)

        if isnan(stdVelocity(sectionIdx, tt))
            stdVelocity(sectionIdx, tt) = 0;
        end

        avgVolumeRate(sectionIdx, tt) = topVelocity(sectionIdx, tt)/2 * crossSectionArea(sectionIdx) * 60; % microL/min
        stdVolumeRate(sectionIdx, tt) = stdVelocity(sectionIdx, tt) * crossSectionArea(sectionIdx) * 60; % microL/min

        %     figure(101)
        %     plot(plot_values);
        %     findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
        %     % text(locs,pks,num2str((1:numel(pks))'))
        %     text(locs,pks,string(round(avg_blood_rate,3)))
        %     title("Peaks of luminosity")
        %     pbaspect([1.618 1 1]);

        %print(['-f' num2str(70+sectionIdx)],'-dpng',fullfile(ToolBox.path_png,strcat(ToolBox.main_foldername,['_Artery_Section_' num2str(sectionIdx) '.png']))) ;
        %print(['-f' num2str(1000+sectionIdx)],'-dpng',fullfile(ToolBox.path_png,strcat(ToolBox.main_foldername,['_Proj_Artery_Section_' num2str(sectionIdx) '.png']))) ;

    end

    velocityProfiles{sectionIdx} = profils;
    stdVelocityProfiles{sectionIdx} = stdProfils;

    avgVolumeRate(sectionIdx, :) = filloutliers(avgVolumeRate(sectionIdx, :), 'linear');
    stdVolumeRate(sectionIdx, tt) = sqrt(stdVelocity(sectionIdx, tt) ^ 2 * stdCrossSectionArea(sectionIdx) ^ 2 + stdVelocity(sectionIdx, tt) ^ 2 * crossSectionArea(sectionIdx) ^ 2 + stdCrossSectionArea(sectionIdx) ^ 2 * avgVelocity(sectionIdx, tt) ^ 2) * 60; % microL/min

end % sectionIdx

close all
end
