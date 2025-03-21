function obj = VideoInterpolating(obj) %ref = TRUE indicates the object is the reference
[numX, numY, numFrames] = size(obj.M0_ff_video);
params = Parameters_json(obj.directory, obj.param_name);
kInterp = params.json.Preprocess.InterpolationFactor;
numX = (numX - 1) * (2 ^ kInterp - 1) + numX;
numY = (numY - 1) * (2 ^ kInterp - 1) + numY;

if kInterp == 0
    return
end

% Reference M0
tmp_M0_ff = zeros(numX, numY, numFrames);
tmp_Calc_M0_ff = obj.M0_ff_video;

parfor frameIdx = 1:numFrames
    tmp_M0_ff(:, :, frameIdx) = interp2(tmp_Calc_M0_ff(:, :, frameIdx), kInterp);
end

obj.M0_ff_video = tmp_M0_ff;
clear tmp_M0_ff tmp_Calc_M0_ff

% M1M0
tmpM1M0 = zeros(numX, numY, numFrames);
tmpCalcM1M0 = obj.f_AVG_video;

parfor frameIdx = 1:numFrames
    tmpM1M0(:, :, frameIdx) = interp2(tmpCalcM1M0(:, :, frameIdx), kInterp);
end

obj.f_AVG_video = tmpM1M0;
clear tmpM1M0 tmpCalcM1M0

% M2M0
tmpM2M0 = zeros(numX, numY, numFrames);
tmpCalcM2M0 = obj.f_RMS_video;

parfor frameIdx = 1:numFrames
    tmpM2M0(:, :, frameIdx) = interp2(tmpCalcM2M0(:, :, frameIdx), kInterp);
end

obj.f_RMS_video = single(tmpM2M0);
clear tmpM2M0 tmpCalcM2M0

end
