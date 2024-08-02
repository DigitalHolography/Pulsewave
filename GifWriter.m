classdef GifWriter
    %GifWriter Handles the creation and management of Gifs
    %   Detailed explanation goes here

    properties
        filename_gif
        gifLength
        images
        time_period
        time_period_min
        Nx
        Ny
        isRGB
    end

    methods

        function obj = GifWriter(filename, time_period_min, ToolBox)
            %GifWriter Construct an instance of this class
            %   filename: where want your Gif to be built
            %   time_period_min: minimal time between each frame of your GIF

            obj.filename_gif = fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, filename));
            obj.time_period_min = time_period_min;
            obj.time_period = ToolBox.stride / ToolBox.fs / 1000;

        end

        function obj = write(obj, frame)
            % Appends a frame to the gif
            if isa(frame, 'struct')
                image = frame2im(frame);
            else
                image = frame;
            end

            if size(image, 3) == 3
                obj.isRGB = true;
            end

            obj.images = cat(4, obj.images, image);

            obj.gifLength = size(obj.images, 4);
            obj.Nx = size(image, 1);
            obj.Ny = size(image, 2);

        end

        function obj = generate(obj)
            % Generate the gif from the current array of frames            
            h = waitbar(0, 'Generate GIF file...');
            
            if obj.time_period < obj.time_period_min

                num_T = floor(obj.gifLength * obj.time_period / obj.time_period_min);

                if obj.isRGB
                    images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.Nx obj.Ny num_T], "nearest");
                    images_interp(:, :, 2, :) = imresize3(squeeze(obj.images(:, :, 2, :)), [obj.Nx obj.Ny num_T], "nearest");
                    images_interp(:, :, 3, :) = imresize3(squeeze(obj.images(:, :, 3, :)), [obj.Nx obj.Ny num_T], "nearest");
                else
                    images_interp(:, :, 1, :) = imresize3(squeeze(obj.images(:, :, 1, :)), [obj.Nx obj.Ny num_T], "nearest");
                end

                for tt = 1:num_T
                    waitbar((tt - 1) / num_T, h);

                    if obj.isRGB
                        [A, map] = rgb2ind(images_interp(:, :, :, tt), 256);
                    else
                        [A, map] = gray2ind(images_interp(:, :, :, tt), 256);
                    end

                    if tt == 1
                        imwrite(A, map, obj.filename_gif, "gif", "LoopCount", Inf, "DelayTime", obj.time_period_min);
                    else
                        imwrite(A, map, obj.filename_gif, "gif", "WriteMode", "append", "DelayTime", obj.time_period_min);
                    end

                end

            else

                for tt = 1:obj.gifLength
                    if obj.isRGB
                        [A, map] = rgb2ind(obj.images(:, :, :, tt), 256);
                    else
                        [A, map] = gray2ind(obj.images(:, :, :, tt), 256);
                    end

                    waitbar((tt - 1) / obj.gifLength, h);

                    if tt == 1
                        imwrite(A, map, obj.filename_gif, "gif", "LoopCount", Inf, "DelayTime", obj.time_period);
                    else
                        imwrite(A, map, obj.filename_gif, "gif", "WriteMode", "append", "DelayTime", obj.time_period);
                    end

                end

            end

            close(h)
        end

    end

end
