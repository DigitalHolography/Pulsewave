function [maskSection, VesselImageRGB] = createMaskSection(TB, img, r1, r2, xy_barycenter, figname, maskArtery, maskVein,option)

arguments
    TB
    img
    r1
    r2
    xy_barycenter
    figname
    maskArtery
    maskVein = []
    option.thin = 0
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

if option.thin>0
    cercle_mask1 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r1 ;
    cercle_mask11 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r1 - option.thin ;
    maskSection1 = xor(cercle_mask1, cercle_mask11);
    cercle_mask2 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r2;
    cercle_mask21 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r2 - option.thin;
    maskSection2 = xor(cercle_mask2, cercle_mask21);
    maskSection = xor(maskSection1, maskSection2);
else
    cercle_mask1 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r1 ;
    cercle_mask2 = sqrt((x - xx) .^ 2 + (y - yy) .^ 2) <= r2;
    maskSection = xor(cercle_mask1, cercle_mask2);
end

%% Create Colormap Artery/Vein
img = mat2gray(img);
img_RGB = cat(3, img, img, img);
maskSectionRGB = cat(3, maskSection, maskSection, maskSection);

maskArterySection = maskArtery & maskSection;

cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapArtery_Section = cmapLAB(256, [1 1 1], 0, [1 1 0], 1/3, [1 0 0], 2/3, [0 0 0], 1);

M0_Artery = setcmap(img, maskArtery, cmapArtery);
M0_Artery_Section = setcmap(img, maskArterySection, cmapArtery_Section);

if isempty(maskVein)

    VesselImageRGB = ...
        M0_Artery .* ~maskSection + ...
        M0_Artery_Section + ...
        maskSectionRGB .* ~maskArtery + ...
        img_RGB .* ~maskArtery .* ~maskSection;

    if ~isfolder(fullfile(TB.path_png, 'mask'))
        mkdir(fullfile(TB.path_png, 'mask'));
    end

    imwrite(VesselImageRGB, fullfile(TB.path_png, 'mask', sprintf("%s_%s.png", TB.main_foldername, figname)), 'png');

else

    maskVein_Section = maskVein & maskSection;
    maskAV = maskArtery & maskVein;
    maskAV_Section = maskAV & maskSection;

    cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
    cmapVein_Section = cmapLAB(256, [1 1 1], 0, [0 1 1], 1/3, [0 0 1], 2/3, [0 0 0], 1);
    cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);
    cmapAV_Section = cmapLAB(256, [0 0 0], 0, [0 1 0], 1/3, [1 1 1], 1);

    M0_Vein = setcmap(img, maskVein, cmapVein);
    M0_Vein_Section = setcmap(img, maskVein_Section, cmapVein_Section);
    M0_AV = setcmap(img, maskAV, cmapAV);
    M0_AV_Section = setcmap(img, maskAV_Section, cmapAV_Section);

    maskVessel = maskArtery | maskVein;

    VesselImageRGB = ...
        M0_Artery .* ~maskAV .* ~maskSection + ...
        M0_Artery_Section .* ~maskAV_Section + ...
        M0_Vein .* ~maskAV .* ~maskSection + ...
        M0_Vein_Section .* ~maskAV_Section + ...
        M0_AV .* ~maskSection + ...
        M0_AV_Section + ...
        maskSectionRGB .* ~maskVessel + ...
        img_RGB .* ~maskArtery .* ~maskVein .* ~maskSection;

    if ~isfolder(fullfile(TB.path_png, 'mask'))
        mkdir(fullfile(TB.path_png, 'mask'));
    end
    imwrite(VesselImageRGB, fullfile(TB.path_png, 'mask', sprintf("%s_%s.png", TB.main_foldername, figname)), 'png');
end

end