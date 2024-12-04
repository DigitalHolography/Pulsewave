function obj = VideoRemoveOutliers(obj)
%% Outlier Cleaning

PW_params = Parameters_json(obj.directory,obj.PW_param_name);

if ~PW_params.removeOutliers
    return
end

[numX, numY, numFrames] = size(obj.f_RMS_video);
window_size = ceil(numFrames / 50);

tmp_M0_disp_cleaned = zeros(numX, numY, numFrames);
tmp_f_RMS_cleaned = zeros(numX, numY, numFrames);
tmp_f_AVG_cleaned = zeros(numX, numY, numFrames);
tmp_M0_disp = obj.M0_disp_video;
tmp_f_RMS = obj.f_RMS_video;
tmp_f_AVG = obj.f_AVG_video;

parfor xx = 1:numX
    
    for yy = 1:numY
        tmp_M0_disp_cleaned(xx, yy, :) = filloutliers(tmp_M0_disp(xx, yy, :), 'linear', 'movmedian', window_size);
    end
    
end

obj.M0_disp_video = tmp_M0_disp_cleaned;

clear tmp_M0_disp_cleaned tmp_M0_disp

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
