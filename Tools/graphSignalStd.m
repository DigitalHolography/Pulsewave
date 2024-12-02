function graphSignalStd(ToolBox, figId, signal, stdsignal, numFrames, ylabl, xlabl, titl, unit, NameValueArgs)
% Plots on an existing graph the signal and its std

arguments
    ToolBox
    figId
    signal
    stdsignal
    numFrames
    ylabl
    xlabl
    titl
    unit
    NameValueArgs.ylimm double = [min(signal) max(signal)]
    NameValueArgs.cropIndx double = 0
end

mean_signal = mean(signal);

if NameValueArgs.cropIndx > 0
    signal = signal(1:NameValueArgs.cropIndx);
    stdsignal = stdsignal(1:NameValueArgs.cropIndx);
end

Color_std = [0.7, 0.7, 0.7];
figure(figId)
fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
axss = [fullTime(1), fullTime(end), NameValueArgs.ylimm];

if length(signal) ~= numFrames % for a variable length of the signal
    fullTime = fullTime(1:length(signal));
end

curve1 = signal + 0.5 * stdsignal;
curve2 = signal - 0.5 * stdsignal;
tmp_fullTime = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

fill(tmp_fullTime, inBetween, Color_std);
hold on;
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, signal, '-k', 'LineWidth', 2);
yline(mean_signal, '--k', 'LineWidth', 2)
hold off;

ylabel(ylabl)
xlabel(xlabl)
title(sprintf("%s : %02.0f %s", titl, round(mean_signal)), unit)

axis tight;

if length(axss) == 4
    axis(axss);
end

fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])

end
