function extendedPulseAnalysis(M0_ff_video, f_RMS_video, f_AVG_mean, v_RMS, maskArtery, maskVein, maskSection, sysIdxList)
% extendedPulseAnalysis - Performs extended pulse analysis on Doppler data.
% Inputs:
%   M0_ff_video: M0 flat-field corrected video.
%   f_RMS_video: RMS frequency video.
%   f_AVG_mean: Average frequency map.
%   v_RMS: RMS velocity map.
%   maskArtery: Mask for arteries.
%   maskVein: Mask for veins.
%   maskSection: Mask for the region of interest.
%   sysIdxList: List of systolic indices.

% Check if sysIdxList is empty
if isempty(sysIdxList)
    warning('sysIdxList is empty. Skipping extended pulse analysis.');
    return;
end

tic;

% Set parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numFramesInterp = params.PulseAnalysis.OneCycleInterpolationPoints;
[~, ~, numFrames] = size(f_RMS_video);
strXlabel = 'Time (s)';
strYlabel = 'Frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;
folder = 'pulseAnalysis';

% Apply masks to isolate arteries, veins, and background
maskArtery = maskArtery & maskSection;
maskVein = maskVein & maskSection;
maskVessel = maskArtery | maskVein;
maskBackground = ~maskVessel & maskSection;

% Define colors for plotting
cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

%% 3) Display and save raw heatmaps
t3 = tic;

% Display and save AVG frequency heatmap
figure("Visible", "off");
imagesc(f_AVG_mean);
colormap gray;
title('AVG Frequency Map RAW');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'AVG Doppler Frequency (kHz)';
c.Label.FontSize = 12;
axis off;
axis image;
range = clim;
imwrite(rescale(f_AVG_mean), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG.png')), 'png');

clear f_AVG_mean;

% Create and save colorbar for AVG image
colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap gray;
f_AVG_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false);
set(gca, 'LineWidth', 3);
f_AVG_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_AVG_colorbar, 'Title');
titleString = 'AVG Doppler Frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG_colorbar.png')));
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG_colorbar.eps')));

fprintf("    3. Raw heatmaps generation took %ds\n", round(toc(t3)));

%% 4) Calculate raw signals of arteries, background, and veins
tic;

% Compute background signal
background_signal = f_RMS_video .* maskBackground;
background_signal = squeeze(sum(background_signal, [1 2])) / nnz(maskBackground);

% Compute arterial signal
arterial_signal = f_RMS_video .* maskArtery;
arterial_signal = squeeze(sum(arterial_signal, [1 2])) / nnz(maskArtery);

% Compute venous signal if veins analysis is enabled
if veinsAnalysis
    venous_signal = f_RMS_video .* maskVein;
    venous_signal = squeeze(sum(venous_signal, [1 2])) / nnz(maskVein);
end

% Plot raw signals
if veinsAnalysis
    graphSignal('4_signalsRaw', folder, ...
        t, arterial_signal, '-', cArtery, ...
        t, background_signal, ':', cBlack, ...
        t, venous_signal, '-', cVein, ...
        Title = 'Arterial Pulse Waveform and Background Signal', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Artery', 'Background', 'Vein'}, TxtName = {'FullArterialSignal', 'FullBackgroundSignal', 'FullVenousSignal'});
else
    graphSignal('4_signalsRaw', folder, ...
        t, arterial_signal, '-', cArtery, ...
        t, background_signal, ':', cBlack, ...
        Title = 'Arterial Pulse Waveform and Background Signal', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Artery', 'Background'}, TxtName = {'FullArterialSignal', 'FullBackgroundSignal'});
end

fprintf("    4. Calculation of raw signals took %ds\n", round(toc));

%% 5) Smoothing signals
tic;

% Smooth arterial signal
delta_f_arterial_signal = arterial_signal - background_signal;
delta_f_arterial_smooth = smoothdata(delta_f_arterial_signal, 'lowess');

% Plot smoothed arterial signal
graphSignal('5_arterialSignalSmoothed', folder, ...
    t, delta_f_arterial_signal, ':', cArtery, ...
    t, delta_f_arterial_smooth, '-', cArtery, ...
    Title = 'Arterial Signal Smoothing', xlabel = strXlabel, ylabel = strYlabel, ...
    Legends = {'Noisy', 'Robust Linear Regression'});

% Smooth venous signal if veins analysis is enabled
if veinsAnalysis
    delta_f_venous_signal = venous_signal - background_signal;
    delta_f_venous_smooth = smoothdata(delta_f_venous_signal, 'lowess');

    % Plot smoothed venous signal
    graphSignal('5_venousSignalSmoothed', folder, ...
        t, delta_f_venous_signal, ':', cVein, ...
        t, delta_f_venous_smooth, '-', cVein, ...
        Title = 'Smoothed Venous Signal', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Noisy', 'Robust Linear Regression'});
end

% Compute noise and data reliability index
outOfNoiseThreshold = params.PulseAnalysis.OneCycleOutOfNoiseThreshold;
noise = sqrt(abs(delta_f_arterial_signal .^ 2 - delta_f_arterial_smooth .^ 2));
idxOutNoise = find(noise > outOfNoiseThreshold * std(noise));
dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(delta_f_arterial_smooth))));
disp(['        Data Reliability Index 2: ' num2str(dataReliabilityIndex2) ' %']);

% Plot filtered pulse vs. residual
graphSignal('5_filteredPulseVsResidual', folder, ...
    t, delta_f_arterial_smooth, ':', cArtery, ...
    t, noise, '-', cBlack, ...
    Title = 'Signal vs. Noise', xlabel = strXlabel, ...
    Legend = {'Filtered Arterial Pulse', 'Residual'}, ...
    yLines = [0 std(noise) 2 * std(noise) 3 * std(noise)], yLineLabels = {'', '1 std', '2 std', '3 std'}, ...
    TxtName = {'FilteredArterialPulse', 'ResidualArterialPulse'});

fprintf("    5. Smoothing signals took %ds\n", round(toc));

%% 6) Calculation of pulse signal derivative and finding/smoothing pulses
tic;

% Compute derivatives of arterial and venous signals
fullArterialPulseDerivative = gradient(arterial_signal);
fullArterialPulseSmoothDerivative = gradient(delta_f_arterial_smooth);

if veinsAnalysis
    fullVenousPulseDerivative = gradient(venous_signal);
    fullVenousPulseSmoothDerivative = gradient(delta_f_venous_smooth);
end

% Plot derivatives
if veinsAnalysis
    textsX = [t(sysIdxList) - 0.3, t(sysIdxList) - 0.3];
    textsY = [fullArterialPulseSmoothDerivative(sysIdxList)' + 0.03, fullVenousPulseSmoothDerivative(sysIdxList)' + 0.03];
    texts = arrayfun(@num2str, 1:length(sysIdxList) * 2, 'UniformOutput', false);

    graphSignal('6_derivative', folder, ...
        t, fullArterialPulseDerivative, ':', cArtery, ...
        t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
        t, fullVenousPulseDerivative, ':', cVein, ...
        t, fullVenousPulseSmoothDerivative, '-', cVein, ...
        Title = 'Derivative of Arterial and Venous Pulse Waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
        Legend = {'\delta <p(t)> - <b(t)>', 'From Smoothed Data'}, ...
        TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative', 'FullVenousPulseDerivative', 'FullVenousPulseSmoothDerivative'}, ...
        TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts);
else
    textsX = t(sysIdxList) - 0.3;
    textsY = fullArterialPulseSmoothDerivative(sysIdxList)' + 0.03;
    texts = arrayfun(@num2str, 1:length(sysIdxList), 'UniformOutput', false);

    graphSignal('6_derivative', folder, ...
        t, fullArterialPulseDerivative, ':', cArtery, ...
        t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
        Title = 'Derivative of Arterial Waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
        Legend = {'\delta <p(t)> - <b(t)>', 'From Smoothed Data'}, ...
        TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative'}, ...
        TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts);
end

fprintf("    6. Calculation of pulse signal derivative took %ds\n", round(toc));

%% 7) Creation of the average pulse for in-plane arteries
tic;

interp_t = 1:numFramesInterp;

% Create one-cycle average pulse
fprintf("    Average Pulse\n");
[onePulseVideo, ~, ~, onePulseVideoM0] = createOneCycle(f_RMS_video, M0_ff_video, maskArtery, sysIdxList, numFramesInterp);

clear f_RMS_video;

% Create one-cycle average pulse minus background
fprintf("    Average Pulse minus Background\n");

scalingFactor = 1000 * 1000 * 2 * params.json.PulseAnalysis.Lambda / sin(params.json.PulseAnalysis.Phi);
delta_f_RMS = v_RMS / scalingFactor;
[onePulseVideominusBKG, selectedPulseIdx, cycles_signal, ~] = createOneCycle(delta_f_RMS, M0_ff_video, maskArtery, sysIdxList, numFramesInterp);

clear delta_f_RMS;

% Compute average arterial pulse
avgArterialPulseHz = squeeze(sum(onePulseVideominusBKG .* maskArtery, [1 2])) / nnz(maskArtery);
avgArterialPulseVelocityInPlane = avgArterialPulseHz * scalingFactor;

% Export videos if enabled
if exportVideos
    writeGifOnDisc(mat2gray(onePulseVideo), "one_cycle");
    writeGifOnDisc(mat2gray(onePulseVideoM0), "one_cycleM0");
end

% Save average arterial pulse wave velocity to a text file
tmp = [interp_t(1:length(avgArterialPulseHz))', avgArterialPulseVelocityInPlane];
fileID = fopen(fullfile(ToolBox.path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'avgPulse.txt')), 'w');
fprintf(fileID, '%f %f \r\n', tmp');
fclose(fileID);

% Plot RMS Doppler frequency for different cycles
figure("Visible", "off");

for ii = 1:size(cycles_signal, 1)

    if ismember(ii, selectedPulseIdx)
        plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 1);
    else
        plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
    end

    hold on;
end

plot(interp_t, movavgvar(squeeze(mean(cycles_signal(:, :), 1)), 5), 'k-', 'LineWidth', 2);
fontsize(gca, 12, "points");
xlabel('Average Cardiac Cycle Duration (s)');
ylabel('RMS Doppler Frequency (kHz)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '7_RMS_Doppler_frequency_for_different_cycles.png')));
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '7_RMS_Doppler_frequency_for_different_cycles.eps')));

%% Arterial pulse wave analysis
[~, idx_sys] = max(avgArterialPulseHz);

if ~isnan(onePulseVideoM0)
    %% Diastolic and systolic Doppler frequency heatmaps

    gwRatio = params.json.FlatFieldCorrection.GWRatio;
    border = params.json.FlatFieldCorrection.Border;

    heatmap_dia_raw = squeeze(mean(onePulseVideominusBKG(:, :, floor(0.9 * numFramesInterp):numFramesInterp), 3));
    heatmap_dia = flat_field_correction(heatmap_dia_raw, ceil(gwRatio * size(heatmap_dia_raw, 1)), border);

    heatmap_sys_raw = squeeze(mean(onePulseVideominusBKG(:, :, max(ceil(idx_sys - 0.05 * numFramesInterp), 1):min(ceil(idx_sys + 0.05 * numFramesInterp), numFramesInterp)), 3));
    heatmap_sys = flat_field_correction(heatmap_sys_raw, ceil(gwRatio * size(heatmap_sys_raw, 1)), border);

    % Plot and save diastolic and systolic heatmaps
    figure("Visible", "off");
    imagesc(heatmap_dia_raw);
    colormap gray;
    title('Bottom Diastole RMS Frequency Map');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler Frequency (kHz)';
    c.Label.FontSize = 12;
    axis off;
    axis image;
    range = clim;

    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFig.png')));
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFig.eps')));

    figure("Visible", "off");
    imagesc(heatmap_dia);
    colormap gray;
    title('Bottom Diastole RMS Frequency Map (Flatfield)');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler Frequency (kHz)';
    c.Label.FontSize = 12;
    axis off;
    axis image;

    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFlatfieldFig.png')));
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFlatfieldFig.eps')));

    figure("Visible", "off");
    imagesc(heatmap_sys_raw);
    colormap gray;
    title('Peak Systole RMS Frequency Map');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler Frequency (kHz)';
    c.Label.FontSize = 12;
    axis off;
    axis image;
    clim([min(range), max(range)]);

    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFig.png')));
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFig.eps')));

    figure("Visible", "off");
    imagesc(heatmap_sys);
    colormap gray;
    title('Peak Systole RMS Frequency Map (Flatfield)');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler Frequency (kHz)';
    c.Label.FontSize = 12;
    axis off;
    axis image;
    clim([min(range), max(range)]);

    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFlatfieldFig.png')));
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFlatfieldFig.eps')));

    % Save heatmaps as images
    imwrite(rescale(heatmap_sys), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMap.png')), 'png');
    imwrite(rescale(heatmap_dia), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMap.png')), 'png');
end

%% Cleanup and close figures
clear idx_sys max_diff_pulse acc pulse_arteries_blurred_sys diff_pulse_sys pulse_arteries_blurred_dia diff_pulse_dia diff_avgPulse;
close all;
end
