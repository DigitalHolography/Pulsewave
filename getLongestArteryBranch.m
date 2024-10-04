function [maskLongArtery] = getLongestArteryBranch(maskArtery,ToolBox,path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
PW_params = Parameters_json(path);
[numX, numY] = size(maskArtery);

skel = bwskel(maskArtery);
[x, y] = meshgrid(1:numY, 1:numX);
cercleMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius/1.5 * (numY + numX) / 2;

skel = skel & ~cercleMask; % takes out arteries near the center 
[L,n] = bwlabel(skel&~imdilate(bwmorph(skel, 'branchpoints'), strel('disk', 2)));

% FIX ME we need a graph based approach to get full successive branches.

% simply getting the longuest individual branch for now
mx =0;
index=1;
for i=1:n % for each individual branch
    s = sum(L==i,'all'); % number of pixel of this branch
    if s>mx
        mx = s;
        index = i;
    end
end

for i=1:n % for each individual branch
    sk_mask = L==i;
    L(imdilate(sk_mask, strel('disk', floor(numX * PW_params.masks_radius/5)))&maskArtery) = i;
end

figure()
imagesc(L);
foldername = ToolBox.main_foldername;
imwrite(L, jet(max(L,[], 'all') + 1), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSeparateLabeled.png')), 'png');



maskLongArtery = (L==index);





end