function [maskLongArtery, L, adjMatrix] = getLongestArteryBranch(maskArtery, U)
% Returns the longest artery branch mask
PW_params = Parameters_json(ToolBox.PW_path);
[numX, numY] = size(maskArtery);

%% Skeletonize and label the individual branches
skel = bwskel(maskArtery);
[x, y] = meshgrid(1:numY, 1:numX);
cercleMask = sqrt((x - ToolBox.x_barycentre) .^ 2 + (y - ToolBox.y_barycentre) .^ 2) <= PW_params.masks_radius / 6 * (numY + numX) / 2; % making a mask to cut the center

skel = skel & ~cercleMask; % takes out arteries near the center
[L, n] = bwlabel(skel & ~imdilate(bwmorph(skel, 'branchpoints'), strel('disk', 2))); % labeling the skeleton with branch points off to get individual branches

% FIX ME we need a graph based approach to get full successive branches.

% simply getting the longuest individual branch for now
mx = 0;
index = 1;

for i = 1:n % for each individual branch
    s = sum(L == i, 'all'); % number of pixel of this branch
    
    if s > mx
        mx = s;
        index = i;
    end
    
end

%% Get the label mask back to initial size
for i = 1:n % for each individual branch
    sk_mask = L == i;
    L(bwareafilt(imdilate(sk_mask, strel('disk', floor(numX * PW_params.masks_radius / 5))) & maskArtery, 1)) = i;
end

figure(72)
imagesc(L);
foldername = ToolBox.main_foldername;
imwrite(L, jet(max(L, [], 'all') + 1), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskSeparateBranchesLabeled.png')), 'png');

%% compute the adjency matrix of the artery tree and make the graph
adjMatrix = false(n);
skel2 = bwskel(maskArtery);
bp = find(bwmorph(skel2, 'branchpoints'));

all_circles = false(size(maskArtery)); %display only

for i = 1:length(bp) % iterates over all branch points
    circle = false(size(maskArtery));
    circle(bp(i)) = true;
    circle = imdilate(circle, strel('disk', 10));
    all_circles = all_circles | circle;
    neighboors = find(L & circle);
    
    for j = 1:length(neighboors) % iterates over all the neigbors of that branch point
        
        for k = 1:length(neighboors)
            adjMatrix(L(neighboors(j)), L(neighboors(k))) = true;
        end
        
    end
    
end

figure(402); imshow(all_circles | skel2); saveas(402, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'maskTree.png')));
g = graph(adjMatrix); figure(403); plot(g); saveas(403, fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", foldername, 'graphTree.png')));
%% get the graph's biggest connected component
conn = conncomp(g);
idx_most_frequent_conn = mode(conn);

%% get the closest branch to the CRA
closest_branch = 0;
dist_closest_branch = Inf;

for i = 1:length(bp) % iterates over all branch points
    [xx, yy] = ind2sub(size(maskArtery), bp(i));
    
    if (xx - ToolBox.x_barycentre) ^ 2 + (yy - ToolBox.y_barycentre) ^ 2 < dist_closest_branch && conn(L(bp(i))) == idx_most_frequent_conn
        dist_closest_branch = (xx - ToolBox.x_barycentre) ^ 2 + (yy - ToolBox.y_barycentre) ^ 2;
        closest_branch = L(bp(i));
    end
    
end

%% get the shortest path between branches from the closest branch to the CRA gathering the most pixels

% d = distances(g); % get the distance between each branch
%
% [~,i] = max(d(closest_branch,d(closest_branch,:)<Inf));
% P = shortestpath(g,closest_branch,i);
%
% maskLongArtery = false(size(maskArtery));
% for idx= 1:n
%     if conn(idx) == idx_most_frequent_conn && ismember(idx,P)
%         maskLongArtery = maskLongArtery | (L==idx);
%     end
% end
max_pixels_intensity = 0;
end_path = [];

for idx = 1:n
    
    if conn(idx) == idx_most_frequent_conn
        P = shortestpath(g, closest_branch, idx);
        pixels_intensity = 0;
        
        for j = 1:length(P)
            pixels_intensity = pixels_intensity + sum((L == P(j)), [1, 2]);
        end
        
        if pixels_intensity > max_pixels_intensity
            max_pixels_intensity = pixels_intensity;
            end_path = P;
        end
        
    end
    
end

maskLongArtery = false(size(maskArtery));

for idx = 1:n
    
    if conn(idx) == idx_most_frequent_conn && ismember(idx, end_path)
        maskLongArtery = maskLongArtery | (L == idx);
    end
    
end

%%FIX ME going back to simplest option for reliability
m = 0;

for idx = 1:n
    
    if sum((L == idx), [1, 2]) > m
        maskLongArtery = (L == idx);
        m = sum((L == idx), [1, 2]);
    end
    
end

figure(408); imshow(maskLongArtery);

end
