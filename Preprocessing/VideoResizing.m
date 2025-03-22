function obj = VideoResizing(obj)

params = Parameters_json(obj.directory, obj.param_name);

out_height = params.json.Preprocess.Resize.FrameHeight;
out_width = params.json.Preprocess.Resize.FrameWidth;
out_numFrames = params.json.Preprocess.Resize.VideoLength;

[numX, numY, numFrames] = size(obj.M0_ff_video);

isHeightResized = out_height > 0;
isWidthResized = out_width > 0;
isLengthResized = out_numFrames > 0;

if isHeightResized && ~isWidthResized
    out_width = numX;

elseif ~isHeightResized && isWidthResized
    out_height = numY;

elseif ~isHeightResized && ~isWidthResized
    out_width = max(numX, numY);
    out_height = max(numX, numY);

    if (numX == numY) && ~isLengthResized
        % If the video is spatially isomorphic and no Heigth/Width/Length
        % is inputed then the function has no more job to do
        return
    end

end

if ~isLengthResized
    out_numFrames = numFrames;
end

[Xq, Yq, Zq] = meshgrid(linspace(1, numY, out_width), linspace(1, numX, out_height), linspace(1, numFrames, out_numFrames));
% tmp_ref = zeros(numX, numY, numFrames);

tmp_calc = obj.M0_ff_video;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.M0_ff_video = single(tmp);

tmp_calc = obj.f_AVG_video;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.f_AVG_video = single(tmp);

tmp_calc = obj.f_RMS_video;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.f_RMS_video = single(tmp);

disp(['Resized data cube : ', num2str(out_width), 'x', num2str(out_height), 'x', num2str(out_numFrames)])

end
