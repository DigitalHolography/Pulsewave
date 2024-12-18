function extendedPulseAnalysis(M0_ff_video, f_RMS_video, f_AVG_mean, v_RMS, maskArtery, maskVein, maskSection, sysIdxList)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path, ToolBox.PW_param_name);
numFramesInterp = PW_params.oneCycleNinterp;

tic

[~, ~, numFrames] = size(f_RMS_video);
strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

veinsAnalysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;
folder = 'pulseAnalysis';

maskArtery = maskArtery & maskSection;
maskVein = maskVein & maskSection;
maskVessel = maskArtery | maskVein;
maskBackground = ~maskVessel & maskSection;

cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

%% 3) Display and save raw heatmaps

t3 = tic;
%% 3) 1) Doppler AVG frequency heatmap

%  Doppler AVG frequency heatmap
figure("Visible", "off")
imagesc(f_AVG_mean);
colormap gray
title('AVG frequency map RAW');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'AVG Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(f_AVG_mean), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG.png')), 'png');

clear f_AVG_mean

% Colorbar for AVG image
colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
f_AVG_colorbar = colorbar('north');
clim(range)
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
f_AVG_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_AVG_colorbar, 'Title');
titleString = 'AVG Doppler frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG_colorbar.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_AVG_colorbar.eps')))

%% 3) 2) Doppler RMS frequency heatmap

tic

f_RMS = squeeze(mean(f_RMS_video, 3));

%  Doppler AVG frequency heatmap
figure("Visible", "off")
imagesc(f_RMS);
colormap gray
title('RMS frequency map RAW');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(f_RMS), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_RMS.png')), 'png');

clear f_RMS

% Colorbar for AVG image
colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
f_RMS_colorbar = colorbar('north');
clim(range)
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
f_RMS_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_RMS_colorbar, 'Title');
titleString = 'RMS Doppler frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_colorbarRMSFrequency.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '3_colorbarRMSFrequency.eps')))

fprintf("    3. Raw heatmaps generation took %ds\n", round(toc(t3)))

%% 4 ) Calculate raw signals of arteries, background and veins

tic

background_signal = f_RMS_video .* maskBackground;
background_signal = squeeze(sum(background_signal, [1 2])) / nnz(maskBackground);

arterial_signal = f_RMS_video .* maskArtery;
arterial_signal = squeeze(sum(arterial_signal, [1 2])) / nnz(maskArtery);

if veinsAnalysis
    venous_signal = f_RMS_video .* maskVein;
    venous_signal = squeeze(sum(venous_signal, [1 2])) / nnz(maskVein);
end

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

fprintf("    4. Calculation of raw signals of arteries, background and veins took %ds\n", round(toc))

%% 5) Smoothing signals

tic

delta_f_arterial_signal = arterial_signal - background_signal;
delta_f_arterial_smooth = smoothdata(delta_f_arterial_signal, 'lowess');

graphSignal('5_arterialSignalSmoothed', folder, ...
    t, delta_f_arterial_signal, ':', cArtery, ...
    t, delta_f_arterial_smooth, '-', cArtery, ...
    Title = 'Arterial Signal Smoothing', xlabel = strXlabel, ylabel = strYlabel, ...
    Legends = {'Noisy', 'Robust linear regression'});

if veinsAnalysis
    delta_f_venous_signal = venous_signal - background_signal;
    delta_f_venous_smooth = smoothdata(delta_f_venous_signal, 'lowess');

    graphSignal('5_signalsSmoothed', folder, ...
        t, delta_f_arterial_smooth, '-', cArtery, ...
        t, delta_f_venous_smooth, '-.', cVein, ...
        Title = 'Smoothed Arterial Signal and Venous Signal', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Artery', 'Vein'}, TxtName = {'FullArterialSignalSmoothed', 'FullVenousSignalSmoothed'});

    graphSignal('5_venousSignalSmoothed', folder, ...
        t, delta_f_venous_signal, ':', cVein, ...
        t, delta_f_venous_smooth, '-', cVein, ...
        Title = 'Smoothed venous signal', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Noisy', 'Robust linear regression'});
end

noise = sqrt(abs(delta_f_arterial_signal .^ 2 - delta_f_arterial_smooth .^ 2));
idxOutNoise = find(noise > PW_params.pulseAnal_outNoiseThreshold * std(noise));

dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(delta_f_arterial_smooth))));
disp(['        data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);

if veinsAnalysis
    noise = sqrt(abs(delta_f_venous_signal .^ 2 - delta_f_venous_smooth .^ 2));
    idxOutNoise = find(noise > PW_params.pulseAnal_outNoiseThreshold * std(noise));

    dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(delta_f_arterial_smooth))));
    disp(['        data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);
end

graphSignal('5_filteredPulseVsResidual', folder, ...
    t, delta_f_arterial_smooth, ':', cArtery, ...
    t, noise, '-', cBlack, ...
    Title = 'signal vs. noise', xlabel = strXlabel, ...
    Legend = {'filtered arterial pulse', 'residual'}, ...
    yLines = [0 std(noise) 2 * std(noise) 3 * std(noise)], yLineLabels = {'', '1 std', '2 std', '3 std'}, ...
    TxtName = {'FilteredArterialPulse', 'ResidualArterialPulse'});

fprintf("    5. Smoothing signals took %ds\n", round(toc))

%% 6) Calculation of pulse signal derivative and finding/smoothing pulses

tic

fullArterialPulseDerivative = gradient(arterial_signal);
fullArterialPulseSmoothDerivative = gradient(delta_f_arterial_smooth);

if veinsAnalysis
    fullVenousPulseDerivative = gradient(venous_signal);
    fullVenousPulseSmoothDerivative = gradient(delta_f_venous_smooth);
end

% now cleanup dataCube to create_one_cycle()
% strategy : use average pulse profiles to detect and
% find  noisy frames @ >3 std from zero-mean
% replace them with cleaned data

if veinsAnalysis
    textsX = [t(sysIdxList) - 0.3, t(sysIdxList) - 0.3];
    textsY = [fullArterialPulseSmoothDerivative(sysIdxList)' + 0.03, fullVenousPulseSmoothDerivative(sysIdxList)' + 0.03];
    texts = cell(1, length(sysIdxList) * 2);

    for n = 1:length(sysIdxList)
        texts{n} = num2str(n);
        texts{n + length(sysIdxList)} = num2str(n);
    end

    graphSignal('6_derivative', folder, ...
        t, fullArterialPulseDerivative, ':', cArtery, ...
        t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
        t, fullVenousPulseDerivative, ':', cVein, ...
        t, fullVenousPulseSmoothDerivative, '-', cVein, ...
        Title = 'derivative of the arterial and venous pulse waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
        Legend = {'\delta <p(t)> - <b(t)>', 'from smoothed data'}, ...
        TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative', 'fullVenousPulseDerivative', 'FullVenousPulseSmoothDerivative'}, ...
        TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts);

else
    textsX = t(sysIdxList) - 0.3;
    textsY = fullArterialPulseSmoothDerivative(sysIdxList)' + 0.03;
    texts = cell(1, length(sysIdxList));

    for n = 1:numel(sysIdxList)
        texts{n} = num2str(n);
    end

    graphSignal('6_derivative', folder, ...
        t, fullArterialPulseDerivative, ':', cArtery, ...
        t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
        Title = 'derivative of the arterial waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
        Legend = {'\delta <p(t)> - <b(t)>', 'from smoothed data'}, ...
        TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative'}, ...
        TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts);
end

fprintf("    6.  Calculation of pulse signal derivative and finding/smoothing pulses signals took %ds\n", round(toc))

%% 7) Creation of the avg pulse for In-plane arteries ~5min

tic

interp_t = 1:numFramesInterp;

fprintf("    Average Pulse\n")
[onePulseVideo, ~, ~, onePulseVideoM0] = createOneCycle(f_RMS_video, M0_ff_video, maskArtery, sysIdxList, numFramesInterp);

clear f_RMS_video f_RMS_video

fprintf("    Average Pulse minus Background\n")
delta_f_RMS = v_RMS / ToolBox.ScalingFactorVelocityInPlane;

[onePulseVideominusBKG, selectedPulseIdx, cycles_signal, ~] = createOneCycle(delta_f_RMS, M0_ff_video, maskArtery, sysIdxList, numFramesInterp);

clear delta_f_RMS

avgArterialPulseHz = squeeze(sum(onePulseVideominusBKG .* maskArtery, [1 2])) / nnz(maskArtery);
avgArterialPulseVelocityInPlane = avgArterialPulseHz * ToolBox.ScalingFactorVelocityInPlane;

if exportVideos
    % avi
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideo), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.avi')));
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideoM0), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.avi')));

    % mp4

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideo), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.mp4')), 'MPEG-4');
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideoM0), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.mp4')), 'MPEG-4');
end

%FIXME: M1/M0 and M2/M0 are subject to aliases at 67 kHz

blur_time_sys = ceil(numFramesInterp / PW_params.pulseAnal_blurScaleFactor);
blur_time_dia = ceil(numFramesInterp / PW_params.pulseAnal_blurScaleFactor);

average_cycle_length = 0;
nb_of_averaged_cycles = 0;

if size(sysIdxList, 2) == 1
    nb_of_averaged_cycles = 1;
else

    for ii = 2:size(sysIdxList, 2)
        average_cycle_length = average_cycle_length + (sysIdxList(ii) - sysIdxList(ii - 1));
        nb_of_averaged_cycles = nb_of_averaged_cycles + 1;
    end

end

% save average arterial pulse wave velocity to txt file
tmp = [interp_t(1:length(avgArterialPulseHz))', avgArterialPulseVelocityInPlane];
%size(tmp)
fileID = fopen(fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'avgPulse.txt')), 'w');
fprintf(fileID, '%f %f \r\n', tmp');
fclose(fileID);

figure("Visible", "off")

for ii = 1:size(cycles_signal, 1)

    if ismember(ii, selectedPulseIdx)
        plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 1);
        hold on
    else
        plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
        hold on
    end

end

hold on
plot(interp_t, movavgvar(squeeze(mean(cycles_signal(:, :), 1)), 5), 'k-', 'LineWidth', 2)
fontsize(gca, 12, "points");
xlabel('Average cardiac cycle duration (s)')
ylabel('RMS Doppler Frequency (kHz)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '7_RMS_Doppler_frequency_for_different_cycles.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '7_RMS_Doppler_frequency_for_different_cycles.eps')))

ArterialResistivityIndex(t, v_RMS, maskArtery, 'v', folder)

%% Arterial pulse wave analysis


%% plot pulses in veins and arteries

% find peak systole index
% sys_index_list_one_cycle = find_systole_index(onePulseVideo2)

[~, idx_sys] = max(avgArterialPulseHz);

if ~isnan(onePulseVideoM0)
    %% diastolic Doppler frequency heatmap : 10% of frames before minimum of diastole

    heatmap_dia_raw = squeeze(mean(onePulseVideominusBKG(:, :, floor(0.9 * numFramesInterp):numFramesInterp), 3));
    % onePulseVideo2 : no background correction
    % heatmap_dia = squeeze(mean(onePulseVideo2(:,:,floor(0.9*Ninterp):Ninterp),3));
    % heatmap_dia = flat_field_correction(heatmap_dia, ceil(PW_params.flatField_gwRatio*size(heatmap_dia,1)), PW_params.flatField_borderDMap);
    heatmap_dia = flat_field_correction(heatmap_dia_raw, ceil(PW_params.flatField_gwRatio * size(heatmap_dia_raw, 1)), PW_params.flatField_border);
    % heatmap_dia = imflatfield(heatmap_dia,PW_params.flatField_gwRatio*size(heatmap_dia,1)/2);

    clear heatmap_sys_raw

    %% systolic Doppler frequency heatmap : 10% of frames around peak systole

    a = max(ceil(idx_sys - 0.05 * numFramesInterp), 1);
    b = min(ceil(idx_sys + 0.05 * numFramesInterp), numFramesInterp);
    heatmap_sys_raw = squeeze(mean(onePulseVideominusBKG(:, :, a:b), 3));
    % onePulseVideo2 : no background correction
    % heatmap_sys = squeeze(mean(onePulseVideo2(:,:,a:b),3));
    % heatmap_sys = flat_field_correction(heatmap_sys, ceil(PW_params.flatField_gwRatio*size(heatmap_sys,1)), PW_params.flatField_borderDMap);
    heatmap_sys = flat_field_correction(heatmap_sys_raw, ceil(PW_params.flatField_gwRatio * size(heatmap_sys_raw, 1)), PW_params.flatField_border);
    % heatmap_sys = imflatfield(heatmap_sys,PW_params.flatField_gwRatio*size(heatmap_sys,1)/2);

    clear onePulseVideo
    clear onePulseVideominusBKG

    %% PLOT

    % diastolic Doppler frequency heatmap
    figure("Visible", "off")
    imagesc(heatmap_dia_raw);
    colormap gray
    title('bottom diastole RMS frequency map');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    range(1:2) = clim;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFig.eps')))

    clear heatmap_dia_raw

    % diastolic Doppler frequency heatmap
    f80 = figure("Visible", "off");
    imagesc(heatmap_dia);
    colormap gray
    title('bottom diastole RMS frequency map flatfield');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    f80.Position = [1100 500 380 420];
    range(1:2) = clim;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFlatfieldFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMapFlatfieldFig.eps')))

    % systolic Doppler frequency heatmap
    f89 = figure("Visible", "off");
    f89.Position = [1100 500 380 420];
    imagesc(heatmap_sys_raw);
    colormap gray
    title('peak systole RMS frequency map');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    range(3:4) = clim;
    % same color axis for systolic and diastolic Doppler heatmaps
    clim([min(range), max(range)]);

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFig.eps')))

    % figure("Visible", "off")
    % clim([min(range),max(range)]);

    % systolic Doppler frequency heatmap
    f90 = figure("Visible", "off");
    f90.Position = [1100 500 380 420];
    imagesc(heatmap_sys);
    colormap gray
    title('peak systole RMS frequency map flatfield');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    range(3:4) = clim;
    % same color axis for systolic and diastolic Doppler heatmaps
    clim([min(range), max(range)]);

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFlatfieldFig.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMapFlatfieldFig.eps')))

    %% SAVING IMAGES

    imwrite(rescale(heatmap_sys), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMap.png')), 'png');
    imwrite(rescale(heatmap_dia), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMap.png')), 'png');

    %%
    % FIXME : replace sys + dia blur by homogenous blur ?
    pulse_arteries_blurred_sys = movavgvar(avgArterialPulseHz(1:idx_sys), blur_time_sys);
    diff_pulse_sys = diff(pulse_arteries_blurred_sys);
    pulse_arteries_blurred_dia = movavgvar(avgArterialPulseHz(idx_sys:end), blur_time_dia);
    % diff_pulse_dia = diff(pulse_arteries_blurred_dia);
    diff_avgPulse = diff(movavgvar(avgArterialPulseHz, blur_time_sys)) * ToolBox.ScalingFactorVelocityInPlane * numFramesInterp / interp_t(end);
    delta = max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:));
    thr = max(pulse_arteries_blurred_dia(:)) - delta ./ exp(1);
    idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < thr);
    [max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux

    if idx_T_diff_max == 1
        acc = abs(max_diff_pulse / (interp_t(idx_T_diff_max) - interp_t(idx_T_diff_max + 1)));
    else
        acc = max_diff_pulse / (interp_t(idx_T_diff_max) - interp_t(idx_T_diff_max - 1));
    end

    % computation of average arterial pulse wave parameters
    T_syst = interp_t(idx_sys);
    systole_area = sum(avgArterialPulseHz(1:idx_sys));
    diastole_area = sum(avgArterialPulseHz(idx_sys:end));
    tmp = systole_area;
    systole_area = systole_area / (diastole_area + systole_area);
    diastole_area = diastole_area / (diastole_area + tmp);
    nb_of_detected_systoles = size(sysIdxList, 2);
    [max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
    %[min_diff_pulse, idx_T_diff_min] = min(diff_pulse_dia);%faire ca mieux

    % txt file output with measured pulse wave parameters
    fileID = fopen(fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, '_pulseWaveParameters.txt')), 'w');
    fprintf(fileID, [ ...
        'Value of pulse derivative at the maximum systolic increase :\n%d\n' ...
        'Maximal acceleration (m/s^(-2)) :\n%d\n' ...
        'Time of maximum systolic increase (s) :\n%d\n' ...
        'Time of systolic peak (s) :\n%d\n' ...
        'Time of the intersection between the diastolic descent and the threshold 1/e (s) :\n%d\n' ...
        'Number of detected systoles :\n%d\nNumber of averaged cycles :\n%d\n' ...
        'Area under the systolic rise curve :\n%d\n' ...
        'Area under the diastolic descent curve  :\n%d\n'], ...
        max_diff_pulse, ...
        acc, ...
        interp_t(idx_T_diff_max + 1), ...
        T_syst, ...
        interp_t(idx_sys + idx_list_threshold_dia(1)), ...
        nb_of_detected_systoles, ...
        nb_of_averaged_cycles, ...
        systole_area, ...
        diastole_area);
    fclose(fileID);

    %
    figure("Visible", "off")
    plot(interp_t(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'k-', LineWidth = 2);
    xline(interp_t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(interp_t(idx_T_diff_max + 1), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (2)', 'FontSize', 14);
    xline(interp_t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(interp_t(idx_sys + idx_list_threshold_dia(1)), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (3)', 'FontSize', 14);
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel('blood flow velocity (mm/s)', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.eps')))

    plot2txt(interp_t(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'AvgArterialPulseVelocityInPlane')

    figure("Visible", "off")

    for ii = 1:size(cycles_signal, 1)

        if ismember(ii, selectedPulseIdx)
            plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 2);
            hold on
        else
            plot(interp_t, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
            hold on
        end

    end

    title('arterial Doppler signal ');
    legend('arterial signal ');
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel('Doppler signal (kHz)');
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    h = findobj(gca, 'Type', 'line');

    l = h.XData;
    m = h.YData;
    plot2txt(l, m, 'ArterialDopplerSignal')

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'all_cycles.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'all_cycles.eps')))

    clear cycles_signal

    figure("Visible", "off")
    plot(interp_t, avgArterialPulseHz, 'k.', ...
        interp_t(1:idx_sys), pulse_arteries_blurred_sys(1:idx_sys), 'k-', ...
        interp_t(idx_sys:numFramesInterp), pulse_arteries_blurred_dia(1:(numFramesInterp - idx_sys + 1)), 'k-', ...
        LineWidth = 2);
    xline(interp_t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(interp_t(idx_T_diff_max + 1), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (2)', 'FontSize', 14);
    xline(interp_t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(interp_t(idx_sys + idx_list_threshold_dia(1)), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (3)', 'FontSize', 14);
    %                 yline(1/exp(1),':',LineWidth=2);
    legend('arterial pulse', 'smoothed line');
    fontsize(gca, 12, "points");
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('frequency (kHz)', 'FontSize', 14);
    axis tight;
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    title('average background-corrected RMS frequency in retinal arteries');

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWaveLabeled.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWaveLabeled.eps')))

    plot2txt(interp_t, avgArterialPulseHz, 'AvgArterialPulseHz')

    figure("Visible", "off")
    plot(interp_t(1:end - 1), diff_avgPulse, 'k-', LineWidth = 2);
    x = 0;
    yline(x, ':', LineWidth = 2);
    xline(interp_t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(interp_t(idx_T_diff_max + 1), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (2)', 'FontSize', 14);
    xline(interp_t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(interp_t(idx_sys + idx_list_threshold_dia(1)), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (3)', 'FontSize', 14);
    legend(' arterial pulse');
    fontsize(gca, 12, "points");
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('time derivative (mm.s^{-2})', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    title('Derivative of average arterial pulse wave');
    axis tight

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWaveDerivative.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWaveDerivative.eps')))

    plot2txt(interp_t(1:end - 1), diff_avgPulse, 'AverageArterialPulseWaveDerivative')
end

clear idx_sys
clear max_diff_pulse
clear acc

clear pulse_arteries_blurred_sys
clear diff_pulse_sys
clear pulse_arteries_blurred_dia
clear diff_pulse_dia
clear diff_avgPulse

close all

end
