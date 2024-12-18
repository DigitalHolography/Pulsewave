function obj = VideoCropping(obj)
%Crop a video (matrix dim 3)
PW_params = Parameters_json(obj.directory, obj.PW_param_name);
firstFrame = PW_params.videoStartFrameIndex;
lastFrame = PW_params.videoEndFrameIndex;
[~, ~, numFrames] = size(obj.M0_ff_video);
logs = obj.load_logs;

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

    str_tosave = sprintf('Data cube frame: %s/%s to %s/%s', num2str(firstFrame), num2str(numFrames), num2str(lastFrame), num2str(numFrames));
    logs = strcat(logs, '\r', str_tosave);
else
    disp('Wrong value for the first frame. Set as 1.')
    disp('Wrong value for the last frame. Set as the end.')
    disp(['Data cube frame: 1/', num2str(numFrames), ' to ', num2str(numFrames), '/', num2str(numFrames)])

    str_tosave = sprintf('Wrong value for the first frame. Set as 1. \rWrong value for the last frame. Set as the end. \rData cube frame: 1/%s to %s/%s', num2str(numFrames), num2str(numFrames), num2str(numFrames));
    logs = strcat(logs, '\r\n\n', str_tosave, '\n');
end

obj.load_logs = logs;
end