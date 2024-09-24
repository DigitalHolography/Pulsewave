function [maskLongArtery,L,adjMatrix] = getLongestArteryBranch(maskArtery,ToolBox,path)
% Returns the longest artery branch mask
PW_params = Parameters_json(path);
[Nx, Ny] = size(maskArtery);

%% Skeletonize and label the individual branches
skel = bwskel(maskArtery);
[x, y] = meshgrid(1:Ny, 1:Nx);
cercleMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius/6 * (Ny + Nx) / 2; % making a mask to cut the center

skel = skel & ~cercleMask; % takes out arteries near the center 
[L,n] = bwlabel(skel&~imdilate(bwmorph(skel, 'branchpoints'), strel('disk', 2))); % labeling the skeleton with branch points off to get individual branches

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
%% Get the label mask back to initial size
for i=1:n % for each individual branch
    sk_mask = L==i;
    L(imdilate(sk_mask, strel('disk', floor(Nx * PW_params.masks_radius/5)))&maskArtery) = i;
end

figure(72)
imagesc(L);
foldername = ToolBox.main_foldername;
imwrite(L, jet(max(L,[], 'all') + 1), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSeparateBranchesLabeled.png')), 'png');


%% compute the adjency matrix of the artery tree
adjMatrix = logical(size([n,n]));
% skel2 = bwskel(maskArtery);
% bp = find(bwmorph(skel2, 'branchpoints'));
% 
% for i=1:length(bp) % iterates over all branch points
%     circle = logical(size(maskArtery));
%     circle(bp(i)) = true;
%     circle = imdilate(circle , strel('disk', 2));
%     neighboors = find(L&circle);
%     for j=1:length(neighboors) % iterates over all the neigbors of that branch point 
%         for k=1:length(neighboors)
%             adjMatrix(L(j),L(k))=true;
%         end
%     end
% end




maskLongArtery = (L==index);





end