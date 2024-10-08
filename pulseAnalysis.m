function [v_OneCycle, v_RMS_video, exec_times, total_time] = pulseAnalysis(numFramesInterp, f_RMS_video, f_AVG_video, M2_data_video, M0_data_video, M0_disp_video, sysIdxList, maskArtery, maskVein, maskBackground, ToolBox, path)

% Variable : LocalBKG_artery, Taille : 10631287200 bytes
% Variable : f_AVG_video, Taille : 10631287200 bytes (DEBUT)
% Variable : f_RMS_video, Taille : 10631287200 bytes (DEBUT)
% Variable : maskArtery, Taille : 18849800 bytes (DEBUT)
% Variable : meanIm, Taille : 18849800 bytes  (DEBUT)
% Variable : maskBackground, Taille : 2356225 bytes (DEBUT)
% Variable : maskVein, Taille : 2356225 bytes (DEBUT)
% Variable : variableInfo, Taille : 12898 bytes

exec_times_id = [];
exec_times_time = [];
total_time = 0;

PW_params = Parameters_json(path);
veinsAnalysis = PW_params.veins_analysis;
[numX, numY, numFrames] = size(f_RMS_video);

mkdir(ToolBox.PW_path_png, 'pulseAnalysis')
mkdir(ToolBox.PW_path_eps, 'pulseAnalysis')

%% 1) Display and save raw heatmaps

strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

%% 1) 1) Doppler AVG frequency heatmap

tic

%  Doppler AVG frequency heatmap
figure(10)
imagesc(f_AVG_video);
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
imwrite(rescale(f_AVG_video), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'frequency_AVG.png')), 'png');

clear f_AVG_video

% Colorbar for AVG image
colorfig = figure(11);
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

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'frequency_AVG_colorbar.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'frequency_AVG_colorbar.eps')))

exec_times_id = [exec_times_id, "Doppler AVG frequency heatmap"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% 1) 2) Doppler RMS frequency heatmap

tic

f_RMS = squeeze(mean(f_RMS_video, 3));

%  Doppler AVG frequency heatmap
figure(20)
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
imwrite(rescale(f_RMS), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'frequency_RMS.png')), 'png');

clear f_RMS

% Colorbar for AVG image
colorfig = figure(21);
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

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarRMSFrequency.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'colorbarRMSFrequency.eps')))

exec_times_id = [exec_times_id, "Doppler RMS frequency heatmap"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% 2 ) Calculate raw signals of arteries, background and veins

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
    figure(20)
    plot(t, arterial_signal, '-k', t, background_signal, ':k', t, venous_signal, '-.k', 'LineWidth', 2);
    title('Arterial Pulse Waveform and Background Signal'); % averaged outside of segmented vessels
    legend('Arterial Pulse', 'Background', 'Venous Signal');
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsRaw.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsRaw.eps')))

    plot2txt(t, arterial_signal, 'arterial_signal', ToolBox)
    plot2txt(t, background_signal, 'background_signal', ToolBox)
    plot2txt(t, venous_signal, 'venous_signal', ToolBox)

else
    figure(20)
    plot(t, arterial_signal, '-k', t, background_signal, ':k', 'LineWidth', 2);
    title('Arterial Pulse Waveform and Background Signal'); % averaged outside of segmented vessels
    legend('Arterial Pulse', 'Background');
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsRaw.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsRaw.eps')))

    plot2txt(t, arterial_signal, 'FullArterialSignal', ToolBox)
    plot2txt(t, background_signal, 'FullBackgroundSignal', ToolBox)

end

exec_times_id = [exec_times_id, "Calculate raw signals"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% 3) Smoothing signals

tic

delta_f_arterial_signal = arterial_signal - background_signal;

if veinsAnalysis
    delta_f_venous_signal = venous_signal - background_signal;
end

% smooth trendline data by iterative local linear regression.

delta_f_arterial_smooth = smoothdata(delta_f_arterial_signal, 'lowess');

if veinsAnalysis
    delta_f_venous_smooth = smoothdata(delta_f_venous_signal, 'lowess');
end

figure(30)

if veinsAnalysis
    plot(t, delta_f_arterial_smooth, '-k', t, delta_f_venous_smooth, '-.k', 'LineWidth', 2);
    title('Smoothed Arterial Signal and Venous Signal'); % averaged outside of segmented vessels
    legend('Arterial Signal', 'Venous Signal');
    plot2txt(t, delta_f_arterial_smooth, 'FullArterialSignalSmoothed', ToolBox)
    plot2txt(t, delta_f_venous_smooth, 'FullVenousSignalSmoothed', ToolBox)
else
    plot(t, delta_f_arterial_smooth, '-k', 'LineWidth', 2);
    title('Smoothed Arterial Signal'); % averaged outside of segmented vessels
    legend('Arterial Signal');
    plot2txt(t, delta_f_arterial_smooth, 'FullArterialSignalSmoothed', ToolBox)
end

fontsize(gca, 14, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel(strYlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsSmoothed.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'signalsSmoothed.eps')))

figure(31)
plot(t, delta_f_arterial_signal, ':k', t, delta_f_arterial_smooth, '-k', 'LineWidth', 2);
title('Smoothed arterial signal');
% legend('Noisy signal', 'Robust linear regression');
fontsize(gca, 12, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel(strYlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseArteryFiltered.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseArteryFiltered.eps')))

if veinsAnalysis
    figure(32)
    hold on
    plot(t, delta_f_venous_signal, ':k', t, delta_f_venous_smooth, '-k', 'LineWidth', 2);
    title('Smoothed venous signal');
    % legend('Noisy signal', 'Robust linear regression');
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    box on
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVeinFiltered.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVeinFiltered.eps')))
end


noise = sqrt(abs(delta_f_arterial_signal .^ 2 - delta_f_arterial_smooth .^ 2));
idxOutNoise = find(noise > PW_params.pulseAnal_outNoiseThreshold * std(noise));

dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(delta_f_arterial_smooth))));
disp(['data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);

if veinsAnalysis
    noise = sqrt(abs(delta_f_venous_signal .^ 2 - delta_f_venous_smooth .^ 2));
    idxOutNoise = find(noise > PW_params.pulseAnal_outNoiseThreshold * std(noise));

    dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(delta_f_arterial_smooth))));
    disp(['data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);
end

exec_times_id = [exec_times_id, "Smoothing signals"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% 4) Calculation of pulse signal derivative and finding/smoothing pulses

tic

fullArterialPulseDerivative = diff(arterial_signal);
fullArterialPulseSmoothDerivative = diff(delta_f_arterial_smooth);

if veinsAnalysis
    fullVenousPulseDerivative = diff(venous_signal);
    fullVenousPulseSmoothDerivative = diff(delta_f_venous_smooth);
end

exec_times_id = [exec_times_id, "Calculate pulse derivative"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

% now cleanup dataCube to create_one_cycle()
% strategy : use average pulse profiles to detect and
% find  noisy frames @ >3 std from zero-mean
% replace them with cleaned data

figure(40)

if veinsAnalysis
    hold on
    plot(t(1:length(fullArterialPulseDerivative)), fullArterialPulseDerivative, ':k', t(1:length(fullArterialPulseSmoothDerivative)), fullArterialPulseSmoothDerivative, '-k', 'LineWidth', 2);
    text(t(sysIdxList) - 0.3, fullArterialPulseSmoothDerivative(sysIdxList) + 0.03, num2str((1:numel(sysIdxList))'))
    plot(t(1:length(fullVenousPulseDerivative)), fullVenousPulseSmoothDerivative, ':k', t(1:length(fullVenousPulseSmoothDerivative)), fullVenousPulseSmoothDerivative, '-k', 'LineWidth', 2);
    text(t(sysIdxList) - 0.3, fullVenousPulseSmoothDerivative(sysIdxList) + 0.03, num2str((1:numel(sysIdxList))'))
    title('derivative of the arterial and venous pulse waveform');
    legend('\delta <p(t)> - <b(t)>', 'from smoothed data');
else
    plot(t(1:length(fullArterialPulseDerivative)), fullArterialPulseDerivative, ':k', t(1:length(fullArterialPulseSmoothDerivative)), fullArterialPulseSmoothDerivative, '-k', 'LineWidth', 2);
    text(t(sysIdxList) - 0.3, fullArterialPulseSmoothDerivative(sysIdxList) + 0.03, num2str((1:numel(sysIdxList))'))
    title('derivative of the arterial pulse waveform');
    legend('\delta <p(t)> - <b(t)>', 'from smoothed data');
end

fontsize(gca, 12, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel('A.U.', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'derivative.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'derivative.eps')))

plot2txt(t(1:length(fullArterialPulseDerivative)), fullArterialPulseDerivative, 'FullArterialPulseDerivative', ToolBox)
plot2txt(t(1:length(fullArterialPulseSmoothDerivative)), fullArterialPulseSmoothDerivative, 'FullArterialPulseSmoothDerivative', ToolBox)

clear fullArterialPulseDerivative fullArterialPulseSmoothDerivative

if veinsAnalysis
    plot2txt(t(1:length(fullVenousPulseDerivative)), fullVenousPulseDerivative, 'FullVenousPulseDerivative', ToolBox)
    plot2txt(t(1:length(fullVenousPulseSmoothDerivative)), fullVenousPulseSmoothDerivative, 'FullVenousPulseSmoothDerivative', ToolBox)
    clear fullVenousPulseDerivative fullVenousPulseSmoothDerivative
end

figure(41)
plot(t(1:length(delta_f_arterial_smooth)), delta_f_arterial_smooth, ':k', t(1:length(noise)), noise, '-k', ...
    'LineWidth', 2);
yline(0, ':', {''}, LineWidth = 2);
yline(std(noise), ':', {'1 std'}, LineWidth = 2);
yline(2 * std(noise), ':', {'2 std'}, LineWidth = 2);
yline(3 * std(noise), ':', {'3 std'}, LineWidth = 2);
title('signal vs. noise');
legend('filtered arterial pulse', 'residual');
fontsize(gca, 12, "points");
xlabel(strXlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'filteredPulseVsResidual.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'filteredPulseVsResidual.eps')))

plot2txt(t(1:length(delta_f_arterial_smooth)), delta_f_arterial_smooth, 'FilteredArterialPulse', ToolBox)
plot2txt(t(1:length(noise)), noise, 'ResidualArterialPulse', ToolBox)

%% Local BKG Artery and Veins %~1min

tic

if veinsAnalysis
    local_mask_vessel = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width));
else
    local_mask_vessel = imdilate(maskArtery, strel('disk', PW_params.local_background_width));
end

f_RMS_background = zeros(numX, numY, numFrames);

parfor frameIdx = 1:numFrames
    f_RMS_background(:, :, frameIdx) = single(regionfill(f_RMS_video(:, :, frameIdx), local_mask_vessel));
end

if PW_params.DiffFirstCalculationsFlag == 0 %SIGNED DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    delta_f_RMS = sign(tmp) .* sqrt(abs(tmp));
    clear tmp

elseif PW_params.DiffFirstCalculationsFlag == 1 % DIFFERENCE FIRST

    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    tmp = tmp .* (tmp > 0);
    delta_f_RMS = sqrt(tmp);
    clear tmp

elseif PW_params.DiffFirstCalculationsFlag == 2 % TO BE TESTED

    LocalBKG_vesselM2 = zeros(numX, numY, numFrames);
    LocalBKG_vesselM0 = zeros(numX, numY, numFrames);

    parfor frameIdx = 1:numFrames
        LocalBKG_vesselM2(:, :, frameIdx) = single(regionfill(M2_data_video(:, :, frameIdx), local_mask_vessel));
        LocalBKG_vesselM0(:, :, frameIdx) = single(regionfill(M0_data_video(:, :, frameIdx), local_mask_vessel));
    end

    tmpM2 = M2_data_video - LocalBKG_vesselM2;
    tmpM0 = M0_data_video - LocalBKG_vesselM0;
    tmpM2M0 = tmpM2 ./ tmpM0;
    delta_f_RMS = sign(tmpM2M0) .* sqrt(abs(tmpM2M0));
    clear tmpM2 tmpM0 tmpM2M0

else % DIFFERENCE LAST

    delta_f_RMS = f_RMS_video - f_RMS_background;

end

exec_times_id = [exec_times_id, "Local Backgrounds"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

f18 = figure(18);
f18.Position = [1100 485 350 420];
LocalBackground_in_vessels = mean(f_RMS_background, 3) .* local_mask_vessel + ones(numX, numY) * mean(f_RMS_background, 'all') .* ~local_mask_vessel;
imagesc(LocalBackground_in_vessels);
colormap gray
title('Local Background in vessels');
fontsize(gca, 14, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
imwrite(rescale(LocalBackground_in_vessels), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_in_vessels.png')))

range(1:2) = clim;

parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.avi')));
parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.mp4')), 'MPEG-4');

clear f_RMS_background

%% Creation of the avg pulse for In-plane arteries ~5min

tic

[onePulseVideo, ~, ~, onePulseVideoM0] = createOneCycle(f_RMS_video, M0_data_video, maskArtery, sysIdxList, numFramesInterp, path, ToolBox);

clear f_RMS_video f_RMS_video

[onePulseVideominusBKG, selectedPulseIdx, cycles_signal, ~] = createOneCycle(delta_f_RMS, M0_data_video, maskArtery, sysIdxList, numFramesInterp, path, ToolBox);

%% Velocity Calculation

v_RMS_video = ToolBox.ScalingFactorVelocityInPlane * delta_f_RMS * ToolBox.NormalizationFactor;

clear fullVideoM2M0minusBKGClean delta_f_RMS

%% 

avgArterialPulseHz = squeeze(sum(onePulseVideominusBKG .* maskArtery, [1 2])) / nnz(maskArtery);
avgArterialPulseVelocityInPlane = avgArterialPulseHz * ToolBox.ScalingFactorVelocityInPlane;

v_OneCycle = (onePulseVideominusBKG .* maskArtery + onePulseVideo .* ~maskArtery) * ToolBox.ScalingFactorVelocityInPlane;

% avi
parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideo), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.avi')));
parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideoM0), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.avi')));

% mp4

parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideo), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.mp4')), 'MPEG-4');
parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(onePulseVideoM0), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.mp4')), 'MPEG-4');

%FIXME: M1/M0 and M2/M0 are subject to aliases at 67 kHz

blur_time_sys = ceil(numFramesInterp / PW_params.pulseAnal_blurScaleFactor);
blur_time_dia = ceil(numFramesInterp / PW_params.pulseAnal_blurScaleFactor);

average_cycle_length = 0;
nb_of_averaged_cycles = 0;

if size(sysIdxList, 2) == 1
    average_cycle_length = numFramesInterp;
    nb_of_averaged_cycles = 1;
else

    for ii = 2:size(sysIdxList, 2)
        average_cycle_length = average_cycle_length + (sysIdxList(ii) - sysIdxList(ii - 1));
        nb_of_averaged_cycles = nb_of_averaged_cycles + 1;
    end

    average_cycle_length = average_cycle_length / (length(sysIdxList) - 1);
end

t = linspace(0, (ToolBox.stride / ToolBox.fs * average_cycle_length) / 1000, numFramesInterp); % /1000 to get time in s

% save average arterial pulse wave velocity to txt file
tmp = [t(1:length(avgArterialPulseHz))', avgArterialPulseVelocityInPlane];
%size(tmp)
fileID = fopen(fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'avgPulse.txt')), 'w');
fprintf(fileID, '%f %f \r\n', tmp');
fclose(fileID);

figure(33)

for ii = 1:size(cycles_signal, 1)

    if ismember(ii, selectedPulseIdx)
        plot(t, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 1);
        hold on
    else
        plot(t, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
        hold on
    end

end

hold on
plot(t, movavgvar(squeeze(mean(cycles_signal(:, :), 1)), 5), 'k-', 'LineWidth', 2)
fontsize(gca, 12, "points");
xlabel('Average cardiac cycle duration (s)')
ylabel('RMS Doppler Frequency (kHz)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'RMS_Doppler_frequency_for_different_cycles.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'RMS_Doppler_frequency_for_different_cycles.eps')))

exec_times_id = [exec_times_id, "Average pulse for In-plane arteries"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% Construct Velocity video ~3min

tic

flowVideoRGB = zeros(numX, numY, 3, numFrames);
referenceMean = mean(M0_disp_video, 3);
v_RMS_mean = mat2gray(squeeze(mean(v_RMS_video(:, :, :), 3)));

if veinsAnalysis
    [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(v_RMS_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
    [hue_vein_mean, sat_vein_mean, val_vein_mean, ~] = createHSVmap(v_RMS_mean, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
    val_mean = v_RMS_mean .* (~(maskArtery + maskVein)) + val_artery_mean .* maskArtery + val_vein_mean .* maskVein  - (val_artery_mean + val_vein_mean)./2.*(maskArtery&maskVein);
    hue_mean = (hue_artery_mean + hue_vein_mean)  - (hue_artery_mean + hue_vein_mean)./2.*(maskArtery&maskVein);
    sat_mean = (sat_artery_mean + sat_vein_mean)  - (sat_artery_mean + sat_vein_mean)./2.*(maskArtery&maskVein);
    flowVideoRGB_mean = hsv2rgb(hue_mean, sat_mean, val_mean);
    flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery + maskVein) + rescale(referenceMean) .* ~(maskArtery + maskVein);
    imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png')))

    parfor ii = 1:numFrames
        v = mat2gray(squeeze(v_RMS_video(:, :, ii)));
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        [hue_vein, sat_vein, val_vein, ~] = createHSVmap(v, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
        val = v .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein - (val_artery + val_vein)./2.*(maskArtery&maskVein);
        hue = (hue_artery + hue_vein)  - (hue_artery + hue_vein)./2.*(maskArtery&maskVein);
        sat = (sat_artery + sat_vein)  - (sat_artery + sat_vein)./2.*(maskArtery&maskVein);
        flowVideoRGB(:, :, :, ii) = hsv2rgb(hue, sat, val);
        flowVideoRGB(:, :, :, ii) = flowVideoRGB(:, :, :, ii) .* (maskArtery + maskVein - maskArtery&maskVein) + rescale(M0_disp_video(:, :, ii)) .* ~(maskArtery + maskVein);
    end

else
    [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(f_RMS_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
    val_mean = f_RMS_mean .* ~maskArtery + val_artery_mean .* maskArtery;
    flowVideoRGB_mean = hsv2rgb(hue_artery_mean, sat_artery_mean, val_mean);
    flowVideoRGB_mean = flowVideoRGB_mean .* maskArtery + rescale(referenceMean) .* ~maskArtery;
    imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png')))

    parfor ii = 1:numFrames
        v = mat2gray(squeeze(v_RMS_video(:, :, ii)));
        [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val = v .* (~(maskArtery)) + val_artery .* maskArtery;
        flowVideoRGB(:, :, :, ii) = hsv2rgb(hue_artery, sat_artery, val);
        flowVideoRGB(:, :, :, ii) = flowVideoRGB(:, :, :, ii) .* (maskArtery) + rescale(M0_disp_video(:, :, ii)) .* ~(maskArtery);
    end

end

clear hue_artery sat_artery val_artery hue_vein sat_vein val_vein

%% Saving video
% avi

parfeval(backgroundPool, @writeVideoOnDisc, 0, flowVideoRGB, fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')));

timePeriod = ToolBox.stride / ToolBox.fs / 1000;

writeGifOnDisc(flowVideoRGB, fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "flowMap")), timePeriod);

% f_RMS_artery = sum(flowVideoRGB .* maskArtery, [1 2]) / nnz(maskArtery);
% f_RMS_vein = sum(flowVideoRGB .* maskVein, [1 2]) / nnz(maskVein);
% f_RMS_max_Arteries = max(f_RMS_artery(:));
% f_RMS_max_Veins = max(f_RMS_vein(:));
% f_RMS_min_Arteries = min(f_RMS_artery(:));
% f_RMS_min_Veins = min(f_RMS_vein(:));

% mp4
parfeval(backgroundPool, @writeVideoOnDisc, 0, flowVideoRGB, fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')), 'MPEG-4');

exec_times_id = [exec_times_id, "Construct velocity video"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% Arterial pulse wave analysis

disp('arterial pulse wave analysis...');

tic

%%
%zero-mean local-to-average arterial pulse cross-correlation
% reference video : average arterial pulse everywhere
% avgArterialPulse_3d = zeros(size(one_pulse_video));
% for mm = 1:size(one_pulse_video, 1)
%     for pp = 1:size(one_pulse_video, 2)
%         avgArterialPulse_3d(mm,pp,:) = avgArterialPulse;
%     end
% end
% [max_C_3, id_max] = zmXCorrVideos(one_pulse_video,avgArterialPulse_3d);

%%

% %     Time lags from cross-correlation between the average arterial pulse and local pulse waveforms
% figure(99);
% imagesc(id_max);
% title('local-to-average arterial pulse lag');
% colormap default
% colorbar;
% axis image
% axis off

% %     cross-correlation value between the average arterial pulse and local pulse waveforms
% figure(77)
% imagesc(log(abs(max_C_3)))
% title('local-to-average arterial pulse maximum cross-correlation value');
% colormap gray;
% colorbar ;
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image

%% artery mask by lag-correlation with average pulse
% arteries = zeros(size(id_max));
% % arteries (or(id_max < 0.1*Ninterp, id_max > 0.9 * Ninterp)) = 1 ;
% arteries (or(id_max <= ceil(0.01*Ninterp), id_max >= floor(0.95 * Ninterp))) = 1 ;
% %     arteries = createArteryMask(one_pulse_video) ;
% figure(10)
% %denoising
% %arteries = medfilt2(arteries,[3 3]);
% imagesc(arteries) ;
% colormap gray
% title('segmented arteries by lag-correlation with average pulse');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image

%% vein mask by lag-correlation with average pulse
% veins = zeros(size(id_max)) ;
% veins(and(id_max>=ceil(0.1*Ninterp),id_max<=floor(0.3*Ninterp))) = 1 ;
% %denoising
% %veins = medfilt2(veins,[3 3]);
% figure(11)
% imagesc(veins) ;
% colormap gray
% title('segmented veins by cross-correlation with average arterial pulse');
% fontsize(gca,12,"points") ;
% set(gca, 'LineWidth', 2);
% axis off
% axis image

%% calculation of average pulse in arteries
% I_arteries = one_pulse_video .* maskArtery;
% avgArterialPulse = squeeze(sum(I_arteries, [1 2]))/nnz(maskArtery);
% avgArterialPulse = avgArterialPulse - min(avgArterialPulse) ;
% [~,idx] = min(avgArterialPulse,[],1) ;
% avgArterialPulse = circshift(avgArterialPulse,-idx);

% %% calculation of pulse wave in veins
% I_veins = one_pulse_video .* veins ;
% pulse_veins = squeeze(sum(I_veins, [1 2]))/nnz(veins);
% pulse_veins = pulse_veins - min(pulse_veins) ;
% pulse_veins = circshift(pulse_veins,-idx) ;
%
% max_plot = max(max(avgArterialPulse(:)),max(pulse_veins(:))) ;
% avgArterialPulse = avgArterialPulse ./ max_plot ;
% pulse_veins = pulse_veins ./ max_plot ;

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

exec_times_id = [exec_times_id, "Arterial Pulsewave analysis"];
exec_times_time = [exec_times_time, toc];
total_time = total_time + toc;

%% PLOT

% diastolic Doppler frequency heatmap
figure(79)
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
f80 = figure(80);
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
f89 = figure(89);
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

% figure(85)
% clim([min(range),max(range)]);

% systolic Doppler frequency heatmap
f90 = figure(90);
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

imwrite(heatmap_sys, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'heatmap_RMS_systol.png')), 'png');
imwrite(heatmap_dia, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'heatmap_RMS_diastole.png')), 'png');
imwrite(rescale(heatmap_sys), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'systoleHeatMap.png')), 'png');
imwrite(rescale(heatmap_dia), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'diastoleHeatMap.png')), 'png');

%%
% FIXME : replace sys + dia blur by homogenous blur ?
pulse_arteries_blurred_sys = movavgvar(avgArterialPulseHz(1:idx_sys), blur_time_sys);
diff_pulse_sys = diff(pulse_arteries_blurred_sys);
pulse_arteries_blurred_dia = movavgvar(avgArterialPulseHz(idx_sys:end), blur_time_dia);
% diff_pulse_dia = diff(pulse_arteries_blurred_dia);
diff_avgPulse = diff(movavgvar(avgArterialPulseHz, blur_time_sys)) * ToolBox.ScalingFactorVelocityInPlane * numFramesInterp / t(end);
delta = max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:));
thr = max(pulse_arteries_blurred_dia(:)) - delta ./ exp(1);
idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < thr);
[max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux

if idx_T_diff_max == 1
    acc = abs(max_diff_pulse / (t(idx_T_diff_max) - t(idx_T_diff_max + 1)));
else
    acc = max_diff_pulse / (t(idx_T_diff_max) - t(idx_T_diff_max - 1));
end

% computation of average arterial pulse wave parameters
T_syst = t(idx_sys);
systole_area = sum(avgArterialPulseHz(1:idx_sys));
diastole_area = sum(avgArterialPulseHz(idx_sys:end));
tmp = systole_area;
systole_area = systole_area / (diastole_area + systole_area);
diastole_area = diastole_area / (diastole_area + tmp);
nb_of_detected_systoles = size(sysIdxList, 2);
[max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux
%[min_diff_pulse, idx_T_diff_min] = min(diff_pulse_dia);%faire ca mieux
if (~exist("nb_of_averaged_cycles"))
    nb_of_averaged_cycles = 0;
end

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
    t(idx_T_diff_max + 1), ...
    T_syst, ...
    t(idx_sys + idx_list_threshold_dia(1)), ...
    nb_of_detected_systoles, ...
    nb_of_averaged_cycles, ...
    systole_area, ...
    diastole_area);
fclose(fileID);

%
figure(100)
plot( ...
    t(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'k-', ...
    LineWidth = 2);
xline(t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
text(t(idx_T_diff_max + 1), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
xline(T_syst, ':', {}, LineWidth = 2);
text(T_syst, min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (2)', 'FontSize', 14);
xline(t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
text(t(idx_sys + idx_list_threshold_dia(1)), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (3)', 'FontSize', 14);
fontsize(gca, 12, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel('blood flow velocity (mm/s)', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.eps')))

plot2txt(t(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'AvgArterialPulseVelocityInPlane', ToolBox)

figure(101)

for ii = 1:size(cycles_signal, 1)

    if ismember(ii, selectedPulseIdx)
        plot( ...
            t, movavgvar(cycles_signal(ii, :), 5), 'k-', ...
            'LineWidth', 2);
        hold on
    else
        plot( ...
            t, movavgvar(cycles_signal(ii, :), 5), 'k--', ...
            'LineWidth', 1);
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
plot2txt(l, m, 'ArterialDopplerSignal', ToolBox)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'all_cycles.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'all_cycles.eps')))

clear cycles_signal

figure(102)
plot(t, avgArterialPulseHz, 'k.', ...
    t(1:idx_sys), pulse_arteries_blurred_sys(1:idx_sys), 'k-', ...
    t(idx_sys:numFramesInterp), pulse_arteries_blurred_dia(1:(numFramesInterp - idx_sys + 1)), 'k-', ...
    LineWidth = 2);
xline(t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
text(t(idx_T_diff_max + 1), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
xline(T_syst, ':', {}, LineWidth = 2);
text(T_syst, min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (2)', 'FontSize', 14);
xline(t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
text(t(idx_sys + idx_list_threshold_dia(1)), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (3)', 'FontSize', 14);
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

plot2txt(t, avgArterialPulseHz, 'AvgArterialPulseHz', ToolBox)

figure(103)
plot(t(1:end - 1), diff_avgPulse, 'k-', LineWidth = 2);
x = 0;
yline(x, ':', LineWidth = 2);
xline(t(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
text(t(idx_T_diff_max + 1), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
xline(T_syst, ':', {}, LineWidth = 2);
text(T_syst, min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (2)', 'FontSize', 14);
xline(t(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
text(t(idx_sys + idx_list_threshold_dia(1)), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (3)', 'FontSize', 14);
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

plot2txt(t(1:end - 1), diff_avgPulse, 'AverageArterialPulseWaveDerivative', ToolBox)
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

exec_times = [exec_times_id; exec_times_time];

return;

%% Analysis CRA & CRV Kept in case we want to study CRA/CRV again
%
% fullCRAPulse = fullVideoM1M0 .* maskCRA;
% fullCRAPulse = squeeze(sum(fullCRAPulse, [1 2]))/nnz(maskCRA);
%
% fullCRVPulse = fullVideoM1M0 .* maskCRV;
% fullCRVPulse = squeeze(sum(fullCRVPulse, [1 2]))/nnz(maskCRV);
%
% %
% fullBackgroundM1M0Signal = fullVideoM1M0 .* maskBackgroundM1M0;
% fullBackgroundM1M0Signal = squeeze(sum(fullBackgroundM1M0Signal, [1 2]))/nnz(maskBackgroundM1M0);
%
%
% figure(102)
% plot(fullTime,fullCRAPulse,'-k', fullTime,fullBackgroundM1M0Signal,':k','LineWidth',2) ;
% title('CRA pulse waveform and backgroundAVG signal'); % averaged outside of segmented vessels
% legend('CRA pulse','backgroundAVG') ;
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
%
% figure(103)
% plot(fullTime,fullCRVPulse,'-k', fullTime,fullBackgroundM1M0Signal,':k','LineWidth',2) ;
% title('CRV pulse waveform and backgroundAVG signal'); % averaged outside of segmented vessels
% legend('CRV pulse','backgroundAVG') ;
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
%
% % Renormalization to avoid bias
% fullCRAPulseMinusBackground = fullCRAPulse - fullBackgroundM1M0Signal;
% fullCRAPulse = fullCRAPulseMinusBackground;
% fullCRVPulseMinusBackground = fullCRVPulse - fullBackgroundM1M0Signal;
% fullCRVPulse = fullCRVPulseMinusBackground;
%
% % le plus simple fonctionne bien : soustraire le bkg.
% A = ones(size(fullVideoM1M0));
% % B = ones(size(dataCube));
% for pp = 1:size(fullVideoM1M0,3)
% %     A(:,:,pp) = A(:,:,pp) * fullBackgroundSignal(pp) * avgFullArterialPulse / avgFullBackgroundSignal;
%       A(:,:,pp) = A(:,:,pp) * fullBackgroundM1M0Signal(pp);
% %     B(:,:,pp) = B(:,:,pp) * fullBackgroundSignal(pp);
% end
% fullVideoM1M0 = fullVideoM1M0 - A ;
% % dataCube = dataCube - A .* maskVessel;
% % dataCube = dataCube - B .* maskBackground;
%
%
% % remove outliers
% % 1st pass
% disp('remove outliers... 1st pass.');
% [idxOutPw,fullCRAPulseRmOut] = discardPulseWaveOutliers(fullCRAPulse,3);
% [idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundM1M0Signal,3);
% dataReliabilityIndex1 = ceil(100*(1-0.5*(length(idxOutPw)/length(fullCRAPulse) + length(idxOutBkg)/length(fullBackgroundM1M0Signal))));
% disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);
%
% [idxOutPw,fullCRVPulseRmOut] = discardPulseWaveOutliers(fullCRVPulse,3);
% [idxOutBkg,fullBackgroundSignalRmOut] = discardPulseWaveOutliers(fullBackgroundM1M0Signal,3);
% dataReliabilityIndex1 = ceil(100*(1-0.5*(length(idxOutPw)/length(fullCRVPulse) + length(idxOutBkg)/length(fullBackgroundM1M0Signal))));
% disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);
%
% % smooth trendline data by iterative local linear regression.
% % fullArterialPulseClean = smoothdata(fullArterialPulse,'rloess');
% fullCRAPulseClean = smoothdata(fullCRAPulseRmOut,'lowess');
% fullCRVPulseClean = smoothdata(fullCRVPulseRmOut,'lowess');
% fullBackgroundSignalClean = smoothdata(fullBackgroundSignalRmOut,'lowess');
%
%
% figure(108)
% plot(fullTime,fullCRAPulse,':k', ...
%     fullTime,fullCRAPulseClean,'-k', ...
%     'LineWidth',2) ;
% title('CRA pulse minus backgroundAVG vs. filtered pulse');
% legend('<p(t)> - <b(t)>','local linear regression');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;
%
% figure(110)
% plot(fullTime,fullCRVPulse,':k', ...
%     fullTime,fullCRVPulseClean,'-k', ...
%     'LineWidth',2) ;
% title('CRV pulse minus backgroundAVG vs. filtered pulse');
% legend('<p(t)> - <b(t)>','local linear regression');
% fontsize(gca,12,"points") ;
% xlabel(strXlabel,'FontSize',14) ;
% pbaspect([1.618 1 1]) ;
% set(gca, 'LineWidth', 2);
% axis tight;

%
% % png
% print('-f108','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRAfilteredPulse.png'))) ;
% print('-f110','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVfilteredPulse.png'))) ;
% print('-f24','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'flattenedDopplerHeatMap.png'))) ;
% print('-f77','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'zeroLagXcorr.png'))) ;
% print('-f99','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'timeLags.png'))) ;
%
% % eps
% print('-f108','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRAfilteredPulse.eps'))) ;
% print('-f110','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVfilteredPulse.eps'))) ;
% print('-f77','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'zeroLagXcorr.eps'))) ;
% print('-f99','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'timeLags.eps'))) ;
%

end
