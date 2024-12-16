function graphCombined(Videofield, mask, etiquettes_locs, etiquettes_values, signal, stdsignal, xy_barycenter, ylabl, xlabl, titl, unit, NameValueArgs)
% Combines a video and a signal to an animated gif combined on top of each
% other frame by frame
% Videofield : Video to be displayed on top grayscale
% mask : mask in red over it
% etiquettes_locs, localization of etiquettes (tags) patch (nbetiquettes 2)
% (leave empty if no etiquettes required)
% etiquettes_values, values inside etiquettes (tags) (nbetiquettes numFrames)
arguments
    Videofield
    mask
    etiquettes_locs
    etiquettes_values
    signal
    stdsignal
    xy_barycenter
    ylabl
    xlabl
    titl
    unit
    NameValueArgs.skip logical = false
    NameValueArgs.Color (1, 3) double = [1 0 0]
end

ToolBox = getGlobalToolBox;
Videofield_rescaled = rescale(Videofield);
[numX, numY, numFrames] = size(Videofield_rescaled);
x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);
ylimm = [min(-1, min(signal)), max(signal) * 1.3];

if NameValueArgs.skip
    startingvalue = numFrames - 1;
else
    startingvalue = 1;
end

for frameIdx = startingvalue:numFrames
    video_plot = figure(410);
    video_plot.Position = [200 200 600 600];
    
    if ~isempty(etiquettes_locs)
        etiquettes_frame_values = round(etiquettes_values(:, frameIdx), 1);
    else
        etiquettes_frame_values = [];
    end
    
    graphMaskTags(video_plot, Videofield_rescaled(:, :, frameIdx), mask, etiquettes_locs, etiquettes_frame_values, x_barycenter, y_barycenter, Color = NameValueArgs.Color);
    title(sprintf("%s : %02.0f %s", titl, round(signal(frameIdx))), unit);
    set(gca, 'FontSize', 14)
    video_plot_frame = getframe(video_plot);
    video_plot_video(:, :, :, frameIdx) = frame2im(video_plot_frame);
end

for frameIdx = startingvalue:numFrames
    signal_plot = figure(530);
    signal_plot.Position = [200 200 600 300];
    
    graphSignalStd(signal_plot, signal, stdsignal, numFrames, ylabl, xlabl, sprintf('Average %s', titl), unit, ylimm = ylimm, cropIndx = frameIdx);
    
    signal_plot_frame = getframe(signal_plot);
    signal_plot_video(:, :, :, frameIdx) = signal_plot_frame.cdata;
    
end

combined_plot_video = cat(1, mat2gray(video_plot_video), mat2gray(signal_plot_video));

if ~NameValueArgs.skip
    writeGifOnDisc(combined_plot_video, sprintf("%s_combined", titl), 0.04);
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(combined_plot_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, sprintf('%s_combined.avi', titl))));
else
    imwrite(combined_plot_video(:, :, :, end), fullfile(ToolBox.PW_path_gif, sprintf("%s_%s_combined.png", ToolBox.PW_folder_name, titl)));
end

close(410, 530)

end
