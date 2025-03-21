function [mask, RG_video] = region_growing_for_vessel(img, seed_map, conditionMask)
% Input
%   img (double, numX x numY) : filtered image to segment
%   seed_map (logical, numX x numY)  : initialisation of the region growing
%   alpha (double) : adjust the region growing
%   conditionMask : pixels were we follow positive gradient
% Output :
%   mask (logical, numX x numY) : result of the region growing segmentation

%% INITIALISATION
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

[numX, numY] = size(img);
mask = seed_map; % will be update at each step ;  represents the pixels already visited
max_iter = round(numX / 4);
floor = params.json.RegionGrowing.FloorThreshold * mean2(img(conditionMask));
alpha = params.json.RegionGrowing.Alpha;

%Represents the seeds points as a list of coordinates
[row, col] = find(seed_map);
seeds_position = [row col];

idx = 1; %index of the first seed not seen yet
last = nnz(mask); %index of the last seed not seen yet

%Chose which pixel are considered as neighbours of a seed
moves = [1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];

%% REGION GROWING
RG_video = zeros(numX, numY, max_iter);

for n = 1:max_iter

    ctr = 0; %To count number of points added to the seeds matrix

    for i = idx:last

        sp = seeds_position(i, :); %Seed Point
        mask(sp(1), sp(2)) = 1; %Marking the seed as visited
        seed_value = img(sp(1), sp(2));

        for m = 1:size(moves, 1)

            new_coord = sp + moves(m, :); %Neighbour Pixel
            x = new_coord(1, 1);
            y = new_coord(1, 2);

            if (x >= 1 && x <= numX && y >= 1 && y <= numY && mask(x, y) == 0) %Checks if seed point is valid or not

                neighbour_value = img(x, y); %value of neighbour
                d = seed_value - neighbour_value;

                if (d > 0 && neighbour_value > floor) || (d > alpha && conditionMask(x, y) == 1) %Serves as Predicate for growing the region
                    %if  (abs(d) > abs(alpha/8) && conditionMask(x,y) == 1 ) %Serves as Predicate for growing the region

                    mask(x, y) = 1; %Marking the pixel as visited
                    seeds_position = cat(1, seeds_position, new_coord(1, :)); %Adding the pixels to seeds
                    ctr = ctr + 1;

                end

            end

        end

    end

    idx = last + 1; %Index of first seed not visted yet
    last = last + ctr; %Updated size of seeds matrix

    if (idx >= last) %Indicates that no pixels can be merged anymore
        break;
    end

    RG_video(:, :, n) = mask;
end

end
