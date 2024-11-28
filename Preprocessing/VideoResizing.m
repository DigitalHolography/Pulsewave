function obj = VideoResizing(obj)
    tic
    PW_params = Parameters_json(obj.directory);
    out_height = PW_params.frameHeight;
    out_width = PW_params.frameWidth;
    out_numFrames = PW_params.videoLength;
    [numX, numY, numFrames] = size(obj.M0_disp_video);

    if out_numFrames < numFrames && out_numFrames > 0
        % we average the input images to get out_numFrames frames

        % batch = floor(numFrames / out_numFrames);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.M0_disp_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.M0_disp_video = single(tmp_ref);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.M0_data_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.M0_data_video = single(tmp_ref);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.M1_data_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.M1_data_video = single(tmp_ref);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.M2_data_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.M2_data_video = single(tmp_ref);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.f_AVG_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.f_AVG_video = single(tmp_ref);

        tmp_ref = zeros([numX, numY, out_numFrames]);
        tmp_calc = obj.f_RMS_video;

        for i = 1:out_numFrames
            tmp_ref(:, :, i) = mean(tmp_calc(:, :, floor(i / out_numFrames * numFrames):floor((i + 1) / out_numFrames * numFrames)), 3);
        end

        obj.f_RMS_video = single(tmp_ref);

    end

    if out_height < 0 && out_width < 0
        if numX ~= numY
            out_width = max(numX,numY);
            out_height = max(numX,numY);
        elseif out_numFrames < 0
            return % do nothing if not required
        end
    end

    if out_height < 0
        out_height = numX;
    end

    if out_width < 0
        out_width = numY;
    end

    if out_numFrames < 0
        out_numFrames = numFrames;
    end

    [Xq, Yq, Zq] = meshgrid(linspace(1, numY, out_width), linspace(1, numX, out_height), linspace(1, numFrames, out_numFrames));
    % tmp_ref = zeros(numX, numY, numFrames);

    tmp_calc_ref = obj.M0_disp_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.M0_disp_video = single(tmp_ref);

    tmp_calc_ref = obj.M0_data_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.M0_data_video = single(tmp_ref);

    tmp_calc_ref = obj.M1_data_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.M1_data_video = single(tmp_ref);

    tmp_calc_ref = obj.M2_data_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.M2_data_video = single(tmp_ref);

    tmp_calc_ref = obj.f_AVG_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.f_AVG_video = single(tmp_ref);

    tmp_calc_ref = obj.f_RMS_video;
    tmp_ref = interp3(tmp_calc_ref, Xq, Yq, Zq);
    obj.f_RMS_video = single(tmp_ref);

    disp(['Resized data cube : ', num2str(out_width), 'x', num2str(out_height), 'x', num2str(out_numFrames)])
    % logs = obj.load_logs;
    % str_tosave = sprintf("Resized data cube : %s x %s x %s", num2str(out_width), num2str(out_height), num2str(out_numFrames));
    % logs = strcat(logs, '\r\n\n', str_tosave, '\n');
    toc
end