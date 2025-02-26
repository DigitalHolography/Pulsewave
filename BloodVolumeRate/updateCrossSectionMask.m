function crossSectionMask = updateCrossSectionMask(crossSectionMask, mask, subImg, locs, sectionIdx, tilt_angle, slice_half_thickness, PW_params)
% Update the cross-section mask for a vessel section.
maskSlice_subImg = false(size(subImg, 1), size(subImg, 2));
slice_center = round(size(subImg, 1) / 2);
slice_half_thickness_tmp = min(slice_half_thickness, floor(size(subImg, 1) / 2));
maskSlice_subImg((slice_center - slice_half_thickness_tmp):(slice_center + slice_half_thickness_tmp), :) = true;
maskSlice_subImg = imrotate(double(maskSlice_subImg), -tilt_angle, 'bilinear', 'crop');
maskSlice_subImg = maskSlice_subImg > PW_params.cropSection_maskThreshold;

maskCurrentSlice = false(size(mask));
maskCurrentSlice(1:size(maskSlice_subImg, 1), 1:size(maskSlice_subImg, 2)) = maskSlice_subImg;
shift_x = locs(sectionIdx, 2) - round(size(maskSlice_subImg, 1) / 2);
shift_y = locs(sectionIdx, 1) - round(size(maskSlice_subImg, 2) / 2);
maskCurrentSlice = circshift(maskCurrentSlice, [shift_y shift_x]);
maskCurrentSlice = maskCurrentSlice .* mask;
crossSectionMask = crossSectionMask + maskCurrentSlice;
end