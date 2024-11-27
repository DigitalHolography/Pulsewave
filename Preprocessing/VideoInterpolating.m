function obj = VideoInterpolating(obj) %ref = TRUE indicates the object is the reference
    [numX, numY, numFrames] = size(obj.M0_disp_video);
    PW_params = Parameters_json(obj.directory);
    kInterp = PW_params.k;
    numX = (numX - 1) * (2 ^ kInterp - 1) + numX;
    numY = (numY - 1) * (2 ^ kInterp - 1) + numY;

    if kInterp == 0
        return
    end

    % Reference M0
    tmpReferenceM0 = zeros(numX, numY, numFrames);
    tmpCalcRef = obj.M0_disp_video;

    parfor frameIdx = 1:numFrames
        tmpReferenceM0(:, :, frameIdx) = interp2(tmpCalcRef(:, :, frameIdx), kInterp);
    end

    obj.M0_disp_video = tmpReferenceM0;
    clear tmpReferenceM0 tmpCalcRef

    % M0
    tmpM0 = zeros(numX, numY, numFrames);
    tmpCalcM0 = obj.M0_data_video;

    parfor frameIdx = 1:numFrames % loop over frames
        tmpM0(:, :, frameIdx) = interp2(tmpCalcM0(:, :, frameIdx), kInterp);
    end

    obj.M0_data_video = tmpM0;
    clear tmpM0 tmpCalcM0

    % M1
    tmpM1 = zeros(numX, numY, numFrames);
    tmpCalcM1 = obj.M1_data_video;

    parfor frameIdx = 1:numFrames % loop over frames
        tmpM1(:, :, frameIdx) = interp2(tmpCalcM1(:, :, frameIdx), kInterp);
    end

    obj.M1_data_video = tmpM1;
    clear tmpM1 tmpCalcM1

    % M2
    tmpM2 = zeros(numX, numY, numFrames);
    tmpCalcM2 = obj.M2_data_video;

    parfor frameIdx = 1:numFrames % loop over frames
        tmpM2(:, :, frameIdx) = interp2(tmpCalcM2(:, :, frameIdx), kInterp);
    end

    obj.M2_data_video = tmpM2;
    clear tmpM2 tmpCalcM2

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