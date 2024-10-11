classdef GifWriter < handle
    %GifWriter Handles the creation and management of Gifs
    %   Detailed explanation goes here

    properties
        filename
        images
        timePeriod
        timePeriodMin
        numX
        numY
        numFrames
        isRGB
    end

    methods

        function obj = GifWriter(filename, timePeriod, timePeriodMin, gifLength)
            %GifWriter Construct an instance of this class
            %   filename: where want your Gif to be built
            %   time_period_min: minimal time between each frame of your GIF

            obj.filename = filename;
            obj.timePeriodMin = timePeriodMin;
            obj.timePeriod = timePeriod; % ToolBox.stride / ToolBox.fs / 1000
            obj.numFrames = gifLength;

        end

        function obj = write(obj, frame, frameIdx)
            % Sets the frame to the gif

            % Checks if it is a frame obj or an image
            if isa(frame, 'struct')
                image = frame2im(frame);
            else
                image = frame;
            end

            if isempty(obj.images)
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

            if obj.timePeriod < obj.timePeriodMin

                num_T = floor(obj.numFrames * obj.timePeriod / obj.timePeriodMin);

                if obj.isRGB
                    images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.numX obj.numY num_T], "nearest");
                    images_interp(:, :, 2, :) = imresize3(squeeze(obj.images(:, :, 2, :)), [obj.numX obj.numY num_T], "nearest");
                    images_interp(:, :, 3, :) = imresize3(squeeze(obj.images(:, :, 3, :)), [obj.numX obj.numY num_T], "nearest");
                else
                    images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.numX obj.numY num_T], "nearest");
                end

                for tt = 1:num_T
                    waitbar((tt - 1) / num_T, h);

                    if obj.isRGB
                        [A, map] = rgb2ind(images_interp(:, :, :, tt), 256);
                    else
                        [A, map] = gray2ind(images_interp(:, :, :, tt), 256);
                    end

                    if tt == 1
                        imwrite(A, map, obj.filename, "gif", "LoopCount", Inf, "DelayTime", obj.timePeriodMin);
                    else
                        imwrite(A, map, obj.filename, "gif", "WriteMode", "append", "DelayTime", obj.timePeriodMin);
                    end

                end

            else

                for tt = 1:obj.numFrames
                    if obj.isRGB
                        [A, map] = rgb2ind(obj.images(:, :, :, tt), 256);
                    else
                        [A, map] = gray2ind(obj.images(:, :, :, tt), 256);
                    end

                    waitbar((tt - 1) / obj.numFrames, h);

                    if tt == 1
                        imwrite(A, map, obj.filename, "gif", "LoopCount", Inf, "DelayTime", obj.timePeriod);
                    else
                        imwrite(A, map, obj.filename, "gif", "WriteMode", "append", "DelayTime", obj.timePeriod);
                    end

                end

            end

            close(h)
        end

    end

end
