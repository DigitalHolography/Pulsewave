function obj = VideoRegistering(obj)
    tic
    % Registers the video using intensity based registration
    PW_params = Parameters_json(obj.directory,obj.PW_param_name);

    if ~PW_params.registerVideoFlag
        return % do nothing if not required
    end

    video = obj.M0_ff_video;
    numX = size(video, 1);
    numY = size(video, 2);
    x = linspace(-numX / 2, numX / 2, numY);
    y = linspace(-numY / 2, numY / 2, numX);
    [X, Y] = meshgrid(x, y);

    disc_ratio = 0.9; % parametrize this coef if needed
    disc = X .^ 2 + Y .^ 2 < (disc_ratio * min(numX, numY) / 2) ^ 2;
    video_reg = video .* disc - disc .* sum(video .* disc, [1, 2]) / nnz(disc); % minus the mean in the disc of each frame
    video_reg = reshape(video_reg, size(video, 1), size(video, 2), 1, size(video, 3)); % insert a dimension to match reegistration functions

    video_reg = video_reg ./ (max(abs(video_reg), [], [1, 2])); % rescaling each frame but keeps mean at zero

    image_ref = mean(video_reg(:, :, PW_params.refAvgStart:PW_params.refAvgEnd), 3); % ref image is from 10 to 20
    [~, shifts] = register_video_from_reference(video_reg, image_ref);

    obj.M0_ff_video = register_video_from_shifts(video, shifts);

    obj.M0_data_video = register_video_from_shifts(obj.M0_data_video, shifts);
    obj.M1_data_video = register_video_from_shifts(obj.M1_data_video, shifts);
    obj.M2_data_video = register_video_from_shifts(obj.M2_data_video, shifts);

    toc

end