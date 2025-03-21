function registered = register_video_from_shifts(video, shifts)
% Registers a redered video from precomputed frame shifts.
% Works with color images as well.
%
% video: a 4D array representing a video
% shifts: a 2 x num_frames array containing translations required to
%         register the video
registered = zeros(size(video), 'single');

for i = 1:size(video, 3)
    registered(:, :, i) = circshift(video(:, :, i), floor(shifts(:, i)));
end

end
