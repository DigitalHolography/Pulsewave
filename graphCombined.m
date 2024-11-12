function graphCombined(Videofield,mask,etiquettes_locs,etiquettes_values,signal,stdsignal,ToolBox,path,ylabl,xlabl,titl,unit,NameValueArgs)
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
    ToolBox
    path
    ylabl
    xlabl
    titl
    unit
    NameValueArgs.skip logical = false
end

Videofield_rescaled = rescale(Videofield);
[numX, numY, numFrames] = size(Videofield_rescaled);
x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;
ylimm = [min(signal),max(signal)];
if NameValueArgs.skip
    startingvalue = numFrames-1;
else
    startingvalue = 1;
end
for frameIdx = startingvalue:numFrames
    video_plot = figure(410);
    video_plot.Position = [200 200 600 600];
    if ~isempty(etiquettes_locs)
        etiquettes_frame_values = round(etiquettes_values(:,frameIdx),1);
    else 
        etiquettes_frame_values = [];
    end
    graphMaskTags(video_plot, Videofield_rescaled(:,:,frameIdx),mask, etiquettes_locs, etiquettes_frame_values,x_center,y_center);
    title(sprintf("%s : %02.0f %s",titl, round(signal(frameIdx))),unit);
    set(gca, 'FontSize', 14)
    video_plot_frame = getframe(video_plot);
    video_plot_video(:, :, :, frameIdx) = frame2im(video_plot_frame);
end

for frameIdx = startingvalue:numFrames
    signal_plot = figure(530);
    signal_plot.Position = [200 200 600 300];
    signal_crop = signal(1:frameIdx);
    stdsignal_crop = stdsignal(1:frameIdx);

    graphSignalStd(signal_plot, signal_crop, stdsignal_crop,ToolBox,numFrames,ylabl,xlabl,sprintf('Average %s',titl),unit,ylimm=ylimm)

    signal_plot_frame = getframe(signal_plot);
    signal_plot_video(:, :, :, frameIdx) = signal_plot_frame.cdata;

end

combined_plot_video = cat(1, mat2gray(video_plot_video), mat2gray(signal_plot_video));

timePeriod = ToolBox.stride / ToolBox.fs / 1000;
if ~NameValueArgs.skip
    writeGifOnDisc(combined_plot_video, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s_combined.gif", ToolBox.PW_folder_name, titl)), timePeriod);
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(combined_plot_video), fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, sprintf('%s_combined.avi',titl))));
else
    imwrite(combined_plot_video(:,:,:,end),fullfile(ToolBox.PW_path_gif, sprintf("%s_%s_combined.png", ToolBox.PW_folder_name, titl)));
end
close(410,530)

end