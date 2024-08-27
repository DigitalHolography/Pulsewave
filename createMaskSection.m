function [maskSection, VesselImageRGB] = createMaskSection(meanIm, maskArtery,radius1,radius2,name_, ToolBox, path)

    PW_params = Parameters_json(path);

    [N, M] = size(maskArtery);
    [x, y] = meshgrid(1:M, 1:N);

    xx = ToolBox.x_barycentre;
    yy = ToolBox.y_barycentre;

    % skel = bwskel(logical(maskArtery), "MinBranchLength", 100); % create skeleton
    %
    % bp = imdilate(bwmorph(skel,'branchpoints'),strel('square',3));
    % ep =  bwmorph(skel,'endpoints');
    %
    % skel_branches = skel - (bp & skel) - ep;
    %
    % skel_branches = bwareaopen(skel_branches,100);
    %
    % cercle_mask = sqrt((x - xx).^2 + (y - yy).^2) <= M/10;
    %
    % skel_dilate = imdilate(skel_branches, strel('disk',2));
    %
    % maskdisk = skel_dilate | cercle_mask;
    % maskdisk = bwareafilt(maskdisk, 1);
    %
    % prim_branches = maskdisk&skel_branches;
    %
    %
    % CC = bwconncomp(prim_branches);
    % k = CC.NumObjects;
    %
    % end_points_rad = zeros(k,2);
    %
    % mask_copy = prim_branches;
    %
    % for i = 1:k
    %     zone_k = bwareafilt(mask_copy, 1);
    %     mask_copy(zone_k==1) = 0;
    %
    %     ep =  bwmorph(zone_k,'endpoints');
    %     [row,col]=find(ep);
    %
    %     if size(row,1) ~=2
    %         end_points_rad(i,1) = 0;
    %         end_points_rad(i,2) = Inf;
    %     else
    %         rad1 =  sqrt((col(1) - xx).^2 + (row(1) - yy).^2);
    %         rad2 =  sqrt((col(2) - xx).^2 + (row(2) - yy).^2);
    %
    %         end_points_rad(i,1) = min(rad1, rad2);
    %         end_points_rad(i,2) = max(rad1, rad2);
    %     end
    % end
    %
    % radius1 = max(end_points_rad(:,1));
    % radius2 = min(end_points_rad(:,2));
    %
    % radius = (radius2+ radius1)/2;
    %
    % % cercle_mask1 = sqrt((x - xx).^2 + (y - yy).^2) <= radius*1.1;
    % % cercle_mask2 = sqrt((x - xx).^2 + (y - yy).^2) <= radius*0.9;

    

    cercle_mask1 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= radius1;
    cercle_mask2 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= radius2;

    maskSection = xor(cercle_mask1, cercle_mask2);

    %% Create Colormap Artery/Vein

    mask_artery = maskArtery;

    [hue_artery, sat_artery, val] = createHSVmap(meanIm, mask_artery - mask_artery .* maskSection, 0, 0);
    [hue_sectionA, sat_sectionA, ~] = createHSVmap(meanIm, maskSection .* mask_artery, 0.15, 0.15);
    sat_section_artery = sat_sectionA;
    val = val .* (~maskSection) + val .* maskSection + maskSection .* (~(mask_artery));
    VesselImageRGB = hsv2rgb(hue_artery + hue_sectionA, sat_artery + sat_section_artery, val);
    mask_all = mask_artery | maskSection;
    VesselImageRGB = VesselImageRGB .* mask_all + ones(size(VesselImageRGB)) .* meanIm .* ~mask_all;

    figure(156)
    imshow(VesselImageRGB);
    title('Mask Artery Section');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    axis off
    axis image

    imwrite(VesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', strcat(ToolBox.main_foldername, name_)), 'png');

end
