function [] = ArterialResistivityIndex(t, v_video, maskArtery, name, folder)

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

ARI = (vMax - vMin) / vMax;
API = (vMax - vMin) / vMean;

% ARI Graph
graphSignal(sprintf('ARI_%s', name), folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('ARI = %0.2f', ARI), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMax], yLineLabels = {'', '', ''});

% API Graph
graphSignal(sprintf('API_%s', name), folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('API = %0.2f', API), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMean, vMax], yLineLabels = {'', '', '', ''});

end
