function obj = VideoRemoveOutliers(obj)
%% Outlier Cleaning

params = Parameters_json(obj.directory, obj.param_name);

if ~params.json.Preprocess.RemoveOutliersOption
    return
end

[numX, numY, numFrames] = size(obj.f_RMS_video);
window_size = params.json.Preprocess.WindowSize;

tmp_M0_ff_cleaned = zeros(numX, numY, numFrames);
tmp_f_RMS_cleaned = zeros(numX, numY, numFrames);
tmp_f_AVG_cleaned = zeros(numX, numY, numFrames);
tmp_M0_ff = obj.M0_ff_video;
tmp_f_RMS = obj.f_RMS_video;
tmp_f_AVG = obj.f_AVG_video;

% Compute the frame-wise mean (spatial average per frame)
frame_means = squeeze(mean(obj.f_RMS_video, [1 2])); % 1D array: numFrames x 1

% Identify outlier frames
outlier_frames_mask = isoutlier(frame_means);

% If no outliers detected, return early
if ~any(outlier_frames_mask)
    return;
end

parfor xx = 1:numX

    for yy = 1:numY
        tmp_M0_ff_cleaned(xx, yy, :) = filloutliers(tmp_M0_ff(xx, yy, :), 'linear', 'movmedian', window_size);
    end

end

obj.M0_ff_video = tmp_M0_ff_cleaned;

clear tmp_M0_ff_cleaned tmp_M0_ff

parfor xx = 1:numX

    for yy = 1:numY
        tmp_f_RMS_cleaned(xx, yy, :) = filloutliers(tmp_f_RMS(xx, yy, :), 'linear', 'movmedian', window_size);
    end

end

obj.f_RMS_video = tmp_f_RMS_cleaned;

clear tmp_f_RMS_cleaned tmp_f_RMS

parfor xx = 1:numX

    for yy = 1:numY
        tmp_f_AVG_cleaned(xx, yy, :) = filloutliers(tmp_f_AVG(xx, yy, :), 'linear', 'movmedian', window_size);
    end

end

obj.f_AVG_video = tmp_f_AVG_cleaned;

clear tmp_f_AVG_cleaned tmp_f_AVG
end
