function [heatmap,lags] = zmXCorrVideos(vid,vidRef)

% zero-mean reference video
vidRef = vidRef - mean(vidRef,3);

% zero-mean video
vid = vid - mean(vid, 3); 

%compute lagged local pulse-to-average pulse xcorrelations
C = zeros(size(vid));
R = zeros(size(vid));
for kk = 1:size(vid, 3)
    R = vid .* circshift(vidRef, kk, 3);
    C(:,:,kk) = squeeze(mean(R, 3));
end

%outputs
[heatmap, lags] = max(C, [], 3);

end

