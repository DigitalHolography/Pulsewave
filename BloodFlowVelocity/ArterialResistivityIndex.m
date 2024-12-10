function [] = ArterialResistivityIndex(t, v_video, M0_ff_video, maskArtery, name, folder)

ToolBox = getGlobalToolBox;

M0_ff_video = rescale(M0_ff_video);
M0_ff_image = rescale(mean(M0_ff_video, 3));

arterial_signal = squeeze(sum(v_video .* maskArtery, [1 2])) / nnz(maskArtery);
arterial_signal = filloutliers(arterial_signal, 'center');
arterial_signal_smooth = smoothdata(arterial_signal, 'rlowess');

%% ARI CALC

[vMin, minIdx] = min(arterial_signal_smooth);
[vMax, maxIdx] = max(arterial_signal_smooth);

ARI = (vMax - vMin) / vMax;

ARImap = squeeze((v_video(:, :, maxIdx) - v_video(:, :, minIdx)) ./ v_video(:, :, maxIdx));
ARImap(ARImap > 1) = 1;
ARImap = ARImap .* (ARImap .* maskArtery > 0);
ARImap(isnan(ARImap)) = 0;

%% ARI FIG

% Color Maps
cArtery = [255 22 18] / 255;
cmapARI = ones(256, 3);
cmapARI(193:256, 2:3) = [linspace(1, 0, 64)' linspace(1, 0, 64)'];
cmapAPI = cmapPerception('rocket');

% ARI Graph
graphSignal('6_ARI', folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('ARI = %0.2f', ARI), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMax], yLineLabels = {'', '', ''});

% ARI Map
ARImapRGB = setcmap(ARImap, maskArtery, cmapARI) .* maskArtery + M0_ff_image .* ~maskArtery;

figure("Visible", "off")
imagesc(ARImapRGB);
title(strcat('Arterial Resistivity Index avg. : ', sprintf(" %3.2f", ARI)));
axis image
axis off
set(gca, 'LineWidth', 2);
fontsize(gca, 14, "points");
c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
c.Label.String = 'Arterial resistivity index';
c.Label.FontSize = 14;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s_ARImapFig.png", ToolBox.main_foldername, name)))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s_ARImapFig.eps", ToolBox.main_foldername, name)))

%% API CALC

[vMean, meanIdx] = mean(arterial_signal_smooth);

API = (vMax - vMin) / vMean;

APImap = squeeze((v_video(:, :, maxIdx) - v_video(:, :, minIdx)) ./ v_video(:, :, meanIdx));
APImap = APImap .* (APImap .* maskArtery > 0);
APImap(isnan(APImap)) = 0;

%% API FIG

% API Graph
graphSignal('6_API', folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('API = %0.2f', API), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, vMin, vMean, vMax], yLineLabels = {'', '', '', ''});

% API Map

maxAPI = max(APImap);
APImapRGB = rescale(setcmap(APImap, maskArtery, cmapAPI)) .* maskArtery + M0_ff_image .* ~maskArtery;

figure("Visible", "off")
imagesc(APImapRGB);
title(strcat('Arterial Pulsatility Index avg. : ', sprintf(" %3.2f", API)));
axis image
axis off
set(gca, 'LineWidth', 2);
fontsize(gca, 14, "points");
c = colorbar('southoutside', 'Ticks', linspace(0, 5, 6), TickLabels, linspace(0, maxAPI, 6));
c.Label.String = 'Arterial pulsatility index';
c.Label.FontSize = 14;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s_APImapFig.png", ToolBox.main_foldername, name)))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s_APImapFig.eps", ToolBox.main_foldername, name)))

end
