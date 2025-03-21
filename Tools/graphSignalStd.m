function graphSignalStd(figId, U, dU, numFrames, ylabl, xlabl, fig_title, unit, NameValueArgs)
% Plots on an existing graph the signal and its std

arguments
    figId
    U % Signal
    dU % Uncertainty of the Signal
    numFrames
    ylabl
    xlabl
    fig_title
    unit
    NameValueArgs.ylimm double = [min(U) max(U)]
    NameValueArgs.cropIndx double = 0
    NameValueArgs.fullTime
end

ToolBox = getGlobalToolBox;
mean_signal = mean(U);

if NameValueArgs.cropIndx > 0
    U = U(1:NameValueArgs.cropIndx);
    dU = dU(1:NameValueArgs.cropIndx);
end

Color_std = [0.7, 0.7, 0.7];
figure(figId);

if ~isempty(ToolBox)
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
else % in a parfor no ToolBox
    fullTime = NameValueArgs.fullTime;
end

axss = [fullTime(1), fullTime(end), NameValueArgs.ylimm];

if length(U) ~= numFrames % for a variable length of the signal
    fullTime = fullTime(1:length(U));
end

curve1 = U + dU;
curve2 = U - dU;
tmp_fullTime = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

fill(tmp_fullTime, inBetween, Color_std);
hold on;
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, U, '-k', 'LineWidth', 2);
yline(mean_signal, '--k', 'LineWidth', 2)
hold off;

ylabel(ylabl)
xlabel(xlabl)
title(sprintf("%s : %02.0f %s", fig_title, round(mean_signal), unit))

axis tight;

if length(axss) == 4
    axis(axss);
end

fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])

end
