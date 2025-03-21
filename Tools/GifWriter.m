classdef GifWriter < handle
%GifWriter Handles the creation and management of Gifs
%   Detailed explanation goes here

properties
    name
    filename_gif
    filename_avi
    filename_mp4
    images
    timePeriod
    timePeriodMin
    numX
    numY
    numFrames
    numFramesFixed
    isRGB
    t
end

methods

    function obj = GifWriter(name, gifLength, timePeriodMin, numFramesFixed, opt)
        %GifWriter Construct an instance of this class
        %   filename: where want your Gif to be built
        %   time_period_min: minimal time between each frame of your GIF

        arguments
            name
            gifLength
            timePeriodMin = NaN
            numFramesFixed = NaN
            opt.ToolBox = [];
        end

        if isempty(opt.ToolBox)
            ToolBox = getGlobalToolBox;
        else
            ToolBox = opt.ToolBox;
        end

        params = ToolBox.getParams;
        obj.name = name;
        obj.filename_gif = fullfile(ToolBox.path_gif, sprintf("%s_%s.gif", ToolBox.folder_name, name));
        obj.filename_avi = fullfile(ToolBox.path_avi, sprintf("%s_%s.avi", ToolBox.folder_name, name));
        obj.filename_mp4 = fullfile(ToolBox.path_mp4, sprintf("%s_%s.mp4", ToolBox.folder_name, name));

        if isnan(timePeriodMin)
            obj.timePeriodMin = params.timePeriodMin;
        else
            obj.timePeriodMin = timePeriodMin;
        end

        obj.numFramesFixed = numFramesFixed;

        obj.timePeriod = ToolBox.stride / ToolBox.fs / 1000;
        obj.numFrames = gifLength;
        obj.t = tic;

    end

    function obj = write(obj, frame, frameIdx)
        % Sets the frame to the gif

        % Checks if it is a frame obj or an image
        if isa(frame, 'struct')
            image = frame2im(frame);
        else
            image = frame;
        end

        if isempty(obj.images) % allocate on first frame
            obj.numX = size(image, 1);
            obj.numY = size(image, 2);

            if size(image, 3) == 3
                obj.isRGB = true;
                obj.images = zeros(obj.numX, obj.numY, 3, obj.numFrames, 'like', image);
            else

                obj.images = zeros(obj.numX, obj.numY, 1, obj.numFrames, 'like', image);
            end

        end

        obj.images(:, :, :, frameIdx) = image;

    end

    function obj = generate(obj)
        % Generate the gif from the current array of frames
        h = waitbar(0, 'Generate GIF file...');

        if obj.numY > 1000
            num_X = round(size(obj.images, 1) * 1000 / size(obj.images, 2));
            num_Y = 1000;
        else
            num_X = obj.numX;
            num_Y = obj.numY;
        end

        if obj.timePeriod < obj.timePeriodMin % in case you ask too fast gif

            if isnan(obj.numFramesFixed)
                num_T = floor(obj.numFrames * obj.timePeriod / obj.timePeriodMin);
            else
                num_T = obj.numFramesFixed;
            end

            if obj.isRGB
                % time-interp
                images_interp_t(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.numX, obj.numY, num_T], "nearest");
                images_interp_t(:, :, 2, :) = imresize3(squeeze(obj.images(:, :, 2, :)), [obj.numX, obj.numY, num_T], "nearest");
                images_interp_t(:, :, 3, :) = imresize3(squeeze(obj.images(:, :, 3, :)), [obj.numX, obj.numY, num_T], "nearest");
                % xy-interp
                images_interp(:, :, 1, :) = imresize3(squeeze(images_interp_t(:, :, 1, :)), [num_X, num_Y, num_T], "cubic");
                images_interp(:, :, 2, :) = imresize3(squeeze(images_interp_t(:, :, 2, :)), [num_X, num_Y, num_T], "cubic");
                images_interp(:, :, 3, :) = imresize3(squeeze(images_interp_t(:, :, 3, :)), [num_X, num_Y, num_T], "cubic");
            else
                % time-interp
                images_interp_t(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.numX, obj.numY, num_T], "nearest");
                % xy-interp
                images_interp(:, :, 1, :) = imresize3(squeeze(images_interp_t(:, :, 1, :)), [num_X, num_Y, num_T], "cubic");
            end

            images_interp(images_interp < 0) = 0;
            images_interp(images_interp > 256) = 256;

            for tt = 1:num_T
                waitbar((tt - 1) / num_T, h);

                if obj.isRGB
                    [A, map] = rgb2ind(images_interp(:, :, :, tt), 256);
                else
                    [A, map] = gray2ind(images_interp(:, :, :, tt), 256);
                end

                if tt == 1
                    imwrite(A, map, obj.filename_gif, "gif", "LoopCount", Inf, "DelayTime", obj.timePeriodMin);
                else
                    imwrite(A, map, obj.filename_gif, "gif", "WriteMode", "append", "DelayTime", obj.timePeriodMin);
                end

            end

            % avi
            parfeval(backgroundPool, @writeVideoOnDisc, 0, images_interp, obj.filename_avi);
            % mp4
            parfeval(backgroundPool, @writeVideoOnDisc, 0, images_interp, obj.filename_mp4, 'MPEG-4');

        else
            num_T = obj.numFrames;

            if obj.isRGB
                % xy-interp
                images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [num_X, num_Y, num_T], "cubic");
                images_interp(:, :, 2, :) = imresize3(squeeze(obj.images(:, :, 2, :)), [num_X, num_Y, num_T], "cubic");
                images_interp(:, :, 3, :) = imresize3(squeeze(obj.images(:, :, 3, :)), [num_X, num_Y, num_T], "cubic");
            else
                % xy-interp
                images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [num_X, num_Y, num_T], "cubic");
            end

            images_interp(images_interp < 0) = 0;
            images_interp(images_interp > 256) = 256;

            for tt = 1:num_T

                if obj.isRGB
                    [A, map] = rgb2ind(images_interp(:, :, :, tt), 256);
                else
                    [A, map] = gray2ind(images_interp(:, :, :, tt), 256);
                end

                waitbar((tt - 1) / num_T, h);

                if tt == 1
                    imwrite(A, map, obj.filename_gif, "gif", "LoopCount", Inf, "DelayTime", obj.timePeriod);
                else
                    imwrite(A, map, obj.filename_gif, "gif", "WriteMode", "append", "DelayTime", obj.timePeriod);
                end

            end

            % avi
            parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(images_interp), obj.filename_avi);
            % mp4
            parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(obj.images), obj.filename_mp4, 'MPEG-4');

        end

        close(h)

        fprintf("    - %s.gif took %ds\n", obj.name, round(toc(obj.t)));

    end

end

end
