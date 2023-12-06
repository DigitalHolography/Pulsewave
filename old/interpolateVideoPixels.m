function videoOut = interpolateVideoPixels(video)

videoOut = zeros(2*size(videoOut, 1)-1, 2*size(videoOut, 2)-1, size(videoOut, 3));
% frame interpolation
for mm = 1:size(video, 3)
        videoOut(:,:,pp) = squeeze(interp2(squeeze(video(:,:,pp)), 1));
end

end


