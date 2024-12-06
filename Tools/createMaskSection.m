function [maskSection, VesselImageRGB] = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, figname, maskArtery, maskVein)

arguments
    ToolBox
    M0_ff_img
    r1
    r2
    xy_barycenter
    figname
    maskArtery
    maskVein = []
end

[numX, numY] = size(maskArtery);
[x, y] = meshgrid(1:numY, 1:numX);

xx = xy_barycenter(1);
yy = xy_barycenter(2);

% skel = bwskel(logical(maskArtery), "MinBranchLength", 100); % create skeleton
%
% bp = imdilate(bwmorph(skel,'branchpoints'),strel('square',3));
% ep =  bwmorph(skel,'endpoints');
%
% skel_branches = skel - (bp & skel) - ep;
%
% skel_branches = bwareaopen(skel_branches,100);
%
% cercle_mask = sqrt((x - xx).^2 + (y - yy).^2) <= numY/10;
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

cercle_mask1 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r1;
cercle_mask2 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r2;

maskSection = xor(cercle_mask1, cercle_mask2);

%% Create Colormap Artery/Vein
M0_ff_img = mat2gray(M0_ff_img);
M0_ff_img_RGB = cat(3, M0_ff_img, M0_ff_img, M0_ff_img);
maskSectionRGB = cat(3, maskSection, maskSection, maskSection);

if size(maskVein, 1) == 0
    [hueA, satA, valA] = createHSVmap(M0_ff_img, ~maskSection & maskArtery, 0, 0);
    [hueSectionA, satSectionA, valSectionA] = createHSVmap(M0_ff_img, maskSection & maskArtery, 0.15, 0.15);
    
    hueArtery = hueA .* (~maskSection & maskArtery) + hueSectionA .* (maskSection & maskArtery);
    satArtery = satA .* (~maskSection & maskArtery) + satSectionA .* (maskSection & maskArtery);
    valArtery = valA .* (~maskSection & maskArtery) + valSectionA .* (maskSection & maskArtery);
    
    VesselImageRGB = hsv2rgb(hueArtery, satArtery, valArtery) + maskSectionRGB .* ~maskArtery + (M0_ff_img_RGB .* ~maskArtery);
    imwrite(VesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s.png", ToolBox.main_foldername, figname)), 'png');
else
    maskVessel = maskArtery | maskVein;
    [hueA, satA, valA] = createHSVmap(M0_ff_img, ~maskSection & maskArtery, 0, 0);
    [hueV, satV, valV] = createHSVmap(M0_ff_img, ~maskSection & maskVein, 0.7, 0.7);
    [hueSectionA, satSectionA, valSectionA] = createHSVmap(M0_ff_img, maskSection & maskArtery, 0.15, 0.15);
    [hueSectionV, satSectionV, valSectionV] = createHSVmap(M0_ff_img, maskSection & maskVein, 0.5, 0.5);
    
    hueVessel = hueA .* (~maskSection & maskArtery) + hueV .* (~maskSection & maskVein) + hueSectionA .* (maskSection & maskArtery) + hueSectionV .* (maskSection & maskVein);
    satVessel = satA .* (~maskSection & maskArtery) + satV .* (~maskSection & maskVein) + satSectionA .* (maskSection & maskArtery) + satSectionV .* (maskSection & maskVein);
    valVessel = valA .* (~maskSection & maskArtery) + valV .* (~maskSection & maskVein) + valSectionA .* (maskSection & maskArtery) + valSectionV .* (maskSection & maskVein);
    
    VesselImageRGB = hsv2rgb(hueVessel, satVessel, valVessel) + maskSectionRGB .* ~maskVessel + (M0_ff_img_RGB .* ~maskVessel);
    imwrite(VesselImageRGB, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s.png", ToolBox.main_foldername, figname)), 'png');
end

end
