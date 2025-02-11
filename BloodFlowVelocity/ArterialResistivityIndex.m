function [] = ArterialResistivityIndex(t, v_video, maskArtery, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps
cArtery = [255 22 18] / 255;

if size(v_video, 3) > 1
    arterial_signal = squeeze(sum(v_video .* maskArtery, [1 2])) / nnz(maskArtery);
else
    arterial_signal = v_video;
end
arterial_signal = filloutliers(arterial_signal, 'center');
arterial_signal_smooth = smoothdata(arterial_signal, 'rlowess');

vMin = min(arterial_signal_smooth);
vMax = max(arterial_signal_smooth);
vMean = mean(arterial_signal_smooth);
vStd = std(arterial_signal);

ARI = (vMax - vMin) / vMax;
API = (vMax - vMin) / vMean;

% ARI Graph
graphSignal(sprintf('ARI_%s', name), folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('ARI %s = %0.2f', name, ARI), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMax], yLineLabels = {'', '', ''});

% API Graph
graphSignal(sprintf('API_%s', name), folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('API %s = %0.2f', name, API), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMean, vMax], yLineLabels = {'', '', '', ''});

fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_', 'PW_main_outputs', '.txt')), 'a');
if strcmp(name,'velocity')
    fprintf(fileID, 'Mean Velocity artery : %f (mm/s) \r\n',vMean);
    fprintf(fileID, 'Std Velocity artery : %f (mm/s) \r\n',vStd);
    fprintf(fileID, 'Max Velocity artery : %f (mm/s) \r\n',vMax);
    fprintf(fileID, 'Min Velocity artery : %f (mm/s) \r\n',vMin);
end
fprintf(fileID, 'Arterial Resistivity Index (%s) : %f  \r\n', name, ARI);
fprintf(fileID, 'Arterial Pulsatility Index (%s) : %f  \r\n', name, API);
fclose(fileID);

end
