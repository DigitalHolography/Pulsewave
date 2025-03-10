function spectrum_video(SH, maskArtery, maskNeighbors)

TB = getGlobalToolBox;

sz = size(SH);
batch_size = sz(3);
numFrames = sz(4);
[numX, numY] = size(maskArtery);
if sz(1)~=numX ||sz(2)~=numY
    maskArtery = logical(imresize(maskArtery, [sz(1), sz(2)]));
    maskNeighbors = logical(imresize(maskNeighbors, [sz(1), sz(2)])) & ~maskArtery;
end

%% save new mask image
SH_neighbors_rgb = zeros(sz(1),sz(2),3);

SH_neighbors_rgb(:,:,1) = maskArtery;
SH_neighbors_rgb(:,:,2) = maskNeighbors;
SH_neighbors_rgb(:,:,3) = zeros(sz(1),sz(2));
if ~isfolder(fullfile(TB.path_png, 'mask'))
    mkdir(fullfile(TB.path_png, 'mask'));
end
imwrite(SH_neighbors_rgb, fullfile(TB.path_png, 'mask', strcat(TB.main_foldername, '_SH_neighbors_rgb.png')));


%% make video
for fri =1:numFrames
    clf;
    spectrum_ploting(SH(:,:,:,fri), maskArtery, maskNeighbors,TB.fs,TB.f1,TB.f2);
    fi=figure(33533);
    frame = getframe(fi);
    spectrum_video(:,:,:,fri) = frame.cdata;
end

writeVideoOnDisc(spectrum_video, fullfile(TB.path_avi, strcat(TB.main_foldername, '_spectrum_video')));

end