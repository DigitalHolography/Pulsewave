function graphCombined(Videofield, mask, signal, stdsignal, xy_barycenter, dirname, opt)
    % Combines a video and a signal into an animated GIF, overlaying them frame by frame.
    %
    % Inputs:
    %   Videofield: Video to be displayed (grayscale).
    %   mask: Mask to overlay in red.
    %   signal: Signal data to plot.
    %   stdsignal: Standard deviation of the signal.
    %   xy_barycenter: Barycenter coordinates [x, y].
    %   dirname: Directory name for saving outputs.
    %   opt: Optional parameters (see below).
    %
    % Optional Parameters (opt):
    %   etiquettes_locs: Locations of tags (empty if not needed).
    %   etiquettes_values: Values for tags (empty if not needed).
    %   ylabl: Y-axis label.
    %   xlabl: X-axis label.
    %   fig_title: Figure title.
    %   unit: Unit for the signal.
    %   skip: Skip frames (logical, default false).
    %   Color: Color for mask overlay.
    %   Visible: Visibility of figures (logical, default false).

    arguments
        Videofield
        mask
        signal
        stdsignal
        xy_barycenter
        dirname
        opt.etiquettes_locs = []
        opt.etiquettes_values = []
        opt.ylabl = ''
        opt.xlabl = ''
        opt.fig_title = ''
        opt.unit = ''
        opt.skip logical = false
        opt.Color = 'none'
        opt.Visible logical = false
    end

    % Get global toolbox settings
    ToolBox = getGlobalToolBox;

    % Rescale video and get dimensions
    Videofield_rescaled = rescale(Videofield);
    [~, ~, numFrames] = size(Videofield_rescaled);
    x_barycenter = xy_barycenter(1);
    y_barycenter = xy_barycenter(2);

    % Set y-axis limits for the signal plot
    ylimm = [min(-1, min(signal)), max(signal) * 1.3];

    % Determine starting frame based on skip option
    if opt.skip
        startingFrame = numFrames;
    else 
        startingFrame = 1;
    end

    % Set parfor argument based on visibility
    parforArg = Inf; % Default (no parallelization if visible)
    if opt.Visible
        parforArg = 0;
    end

    % Precompute etiquettes frame values if provided
    if ~isempty(opt.etiquettes_locs)
        etiquettes_frame_values = round(opt.etiquettes_values(:, 1), 1);
    else
        etiquettes_frame_values = [];
    end

    % Initialize video plot
    etiquette_locs = opt.etiquettes_locs;
    color = opt.Color;

    videoPlot = figure(410);
    videoPlot.Visible = opt.Visible;
    graphMaskTags(videoPlot, Videofield_rescaled(:, :, end), mask, etiquette_locs, etiquettes_frame_values, x_barycenter, y_barycenter, Color = color);
        
    videoPlot.Position = [200 200 600 600];

    % Preallocate video data array
    videoPlotFrames = zeros([size(getframe(gca).cdata), numFrames], 'single');

    % Generate video frames
    parfor (frameIdx = startingFrame:numFrames, parforArg)
        graphMaskTags(videoPlot, Videofield_rescaled(:, :, frameIdx), mask, etiquette_locs, etiquettes_frame_values, x_barycenter, y_barycenter, Color = color);
        set(gca, 'FontSize', 14);
        videoPlotFrames(:, :, :, frameIdx) = frame2im(getframe(gca));
    end

    % Initialize signal plot
    signalPlot = figure(530);
    signalPlot.Visible = opt.Visible;
    signalPlot.Position = [200 200 600 300];

    % Preallocate signal plot data array
    signalPlotFrames = zeros([size(getframe(signalPlot).cdata), numFrames], 'single');

    % Generate signal plot frames
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
    unit = opt.unit;
    fig_title = opt.fig_title;
    x_label = opt.xlabl;
    y_label = opt.ylabl;
    parfor (frameIdx = startingFrame:numFrames, parforArg)
        graphSignalStd(signalPlot, signal, stdsignal, numFrames, y_label, x_label, sprintf('Average %s', fig_title), unit, ylimm = ylimm, cropIndx = frameIdx, fullTime = fullTime);
        signalPlotFrames(:, :, :, frameIdx) = frame2im(getframe(signalPlot));
    end

    % Combine video and signal frames
    videoInterp = imresize(mat2gray(videoPlotFrames), [size(signalPlotFrames, 2), size(signalPlotFrames, 2)]);
    videoInterp = max(0, min(videoInterp, 1)); % Ensure values are within [0, 1]
    combinedFrames = cat(1, videoInterp, mat2gray(signalPlotFrames));

    % Save final frames as PNGs
    imwrite(mat2gray(signalPlotFrames(:, :, :, end)), fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s_plot.png", ToolBox.PW_folder_name, dirname)));
    imwrite(combinedFrames(:, :, :, end), fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s_combined.png", ToolBox.PW_folder_name, dirname)));

    % Save as GIF if not skipping frames
    if ~opt.skip
        writeGifOnDisc(mat2gray(signalPlotFrames), sprintf("%s_plot", dirname), 0.04);
        writeGifOnDisc(combinedFrames, sprintf("%s_combined", dirname), 0.04);
    end

    % Close figures
    close(videoPlot, signalPlot);
end