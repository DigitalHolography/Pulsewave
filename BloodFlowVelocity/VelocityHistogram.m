function histoVideo = VelocityHistogram(v_video, mask, name, n)
% VelocityHistogram.m: Creates an Histogram of the velocities from a video
% and passes the figure as a mat
%
% Inputs:
%   v_video     :   video of the velocities
%   mask        :   a logical mask to locate the velocities
%   name        :   'Arteries' or 'Veins' to change the name of the figs
%               and the color maps
%   n           :   Resolution of the Histrogram
%
% Outputs:
%   histoVideo  :   Mat of the histogram figure

arguments
    v_video
    mask
    name
    n = 256
end

tVelocityVideo = tic;

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);
exportVideos = PW_params.exportVideos;

[numX, numY, numFrames] = size(v_video);

v_histo = v_video .* mask;
v_min = min(v_histo, [], 'all');
v_max = max(v_histo, [], 'all');

if strcmp(name, 'Arteries')
    cmap = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
elseif strcmp(name, 'Veins')
    cmap = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
else
    cmap = cmapLAB(256, [0 0 0], 0, [0 1 0], 1/3, [1/2 1 1/2], 2/3, [1 1 1], 1);
end

yAx = [v_min v_max];

%% Velocity Histogram

X = linspace(v_min, v_max, n);
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
histo = zeros(n, numFrames);
D = (v_max - v_min) / (n - 1);

fDistrib = figure;
fDistrib.Position(3:4) = [600 275];

indexMin = find(X == v_min);
indexMax = find(X == v_max);
imagesc(xAx, yAx, histo(indexMin:indexMax, :))
set(gca, 'YDir', 'normal')
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
colormap(cmap)
f = getframe(gcf);
[numX_fig, numY_fig, ~] = size(f.cdata);

if exportVideos
    histoVideo = zeros(numX_fig, numY_fig, 3, numFrames);
    gifWriter = GifWriter(sprintf("histogramVelocity%s", name), numFrames);

    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if mask(xx, yy) ~= 0
                    i = find(and(X >= v_histo(xx, yy, frameIdx), X < v_histo(xx, yy, frameIdx) + D));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
                end

            end

        end
        figure(fDistrib)
        imagesc(xAx, yAx, histo(indexMin:indexMax, :))
        set(gca, 'YDir', 'normal')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title(sprintf("Velocity distribution in %s", name))
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        f = getframe(gcf);
        histoVideo(:, :, :, frameIdx) = f.cdata;
        gifWriter.write(f, frameIdx);

    end

    gifWriter.generate();
    gifWriter.delete();

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))

    % AVI

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_histogramVelocity%s.mp4", ToolBox.main_foldername, name)));

    tmp = mat2gray(histoVideo);
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);

    % MP4

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, sprintf("%s_histogramVelocity%s.mp4", ToolBox.main_foldername, name)), 'MPEG-4');
    tmp = mat2gray(histoVideo);
    open(w)

    for frameIdx = 1:numFrames
        writeVideo(w, tmp(:, :, :, frameIdx));
    end

    close(w);
else
    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if mask(xx, yy) ~= 0
                    i = find( and(X >= v_histo(xx, yy, frameIdx), X < v_histo(xx, yy, frameIdx) + D));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
                end

            end

        end

    end

    figure(fDistrib)
    imagesc(xAx, yAx, histo(indexMin:indexMax, :))
    set(gca, 'YDir', 'normal')
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    title(sprintf("Velocity distribution in %s", name))
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    f = getframe(gcf);
    histoVideo = f.cdata;
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))

end

fprintf("- Velocity Histogram %s Timing : %ds\n", name, round(toc(tVelocityVideo)))

end