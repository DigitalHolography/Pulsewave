function spectrum_video(SH, maskArtery, maskNeighbors)

ToolBox = getGlobalToolBox;

[~, ~, batch_size, numFrames] = size(SH);
[numX, numY] = size(maskArtery);

if sz(1) ~= numX || sz(2) ~= numY
    maskArtery = logical(imresize(maskArtery, [sz(1), sz(2)]));
    maskNeighbors = logical(imresize(maskNeighbors, [sz(1), sz(2)])) & ~maskArtery;
end

%% save new mask image
SH_neighbors_rgb = zeros(sz(1), sz(2), 3);

SH_neighbors_rgb(:, :, 1) = maskArtery;
SH_neighbors_rgb(:, :, 2) = maskNeighbors;
SH_neighbors_rgb(:, :, 3) = zeros(sz(1), sz(2));

if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
    mkdir(fullfile(ToolBox.path_png, 'mask'));
end

imwrite(SH_neighbors_rgb, fullfile(ToolBox.path_png, 'mask', strcat(ToolBox.main_foldername, '_SH_neighbors_rgb.png')));

spectrum_video = zeros(numX, numY, batch_size, numFrames);
%% make video
for fri = 1:numFrames
    clf;
    spectrum_ploting(SH(:, :, :, fri), maskArtery, maskNeighbors, ToolBox.fs, ToolBox.f1, ToolBox.f2);
    fi = figure(33533);
    frame = getframe(fi);
    spectrum_video(:, :, :, fri) = frame.cdata;
end

writeVideoOnDisc(spectrum_video, fullfile(ToolBox.path_avi, strcat(ToolBox.main_foldername, '_spectrum_video')));

end
