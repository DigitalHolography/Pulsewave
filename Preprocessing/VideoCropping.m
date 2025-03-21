function obj = VideoCropping(obj)
%Crop a video (matrix dim 3)
params = Parameters_json(obj.directory, obj.param_name);
firstFrame = params.videoStartFrameIndex;
lastFrame = params.videoEndFrameIndex;
[~, ~, numFrames] = size(obj.M0_ff_video);

if firstFrame > 0 && firstFrame < numFrames || lastFrame > 1 && lastFrame <= numFrames

    if lastFrame == -1
        lastFrame = numFrames;
    end

    if firstFrame == -1
        firstFrame = 1;
    end

    obj.M0_ff_video = obj.M0_ff_video(:, :, firstFrame:lastFrame);
    obj.M0_data_video = obj.M0_data_video(:, :, firstFrame:lastFrame);
    obj.M1_data_video = obj.M1_data_video(:, :, firstFrame:lastFrame);
    obj.M2_data_video = obj.M2_data_video(:, :, firstFrame:lastFrame);

    disp(['Data cube frame: ', num2str(firstFrame), '/', num2str(numFrames), ' to ', num2str(lastFrame), '/', num2str(numFrames)])
else
    disp('Wrong value for the first frame. Set as 1.')
    disp('Wrong value for the last frame. Set as the end.')
    disp(['Data cube frame: 1/', num2str(numFrames), ' to ', num2str(numFrames), '/', num2str(numFrames)])
end

end
