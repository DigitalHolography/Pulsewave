function [v_RMS_video, exec_times] = pulseAnalysis(numFramesInterp, f_RMS_video, f_AVG_video, M2_data_video, M0_data_video, M0_ff_video, sysIdxList, maskArtery, maskVein, maskBackground, flag_ExtendedPulseWave_analysis, ToolBox, path)

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

PW_params = Parameters_json(path);
veinsAnalysis = PW_params.veins_analysis;
entirePulseAnalysis = flag_ExtendedPulseWave_analysis;
exportVideos = PW_params.exportVideos;
f_AVG_mean = mean(f_AVG_video, 3);

mkdir(ToolBox.PW_path_png, 'pulseAnalysis')
mkdir(ToolBox.PW_path_eps, 'pulseAnalysis')
folder = 'pulseAnalysis';

%% 1) Display and save raw heatmaps

[numX, numY, numFrames] = size(f_RMS_video);
strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

if entirePulseAnalysis

    %% 1) 1) Doppler AVG frequency heatmap

    tic

    %  Doppler AVG frequency heatmap
    figure(10)
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
    imwrite(rescale(f_AVG_mean), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_frequency_AVG.png')), 'png');

    clear f_AVG_mean

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

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_frequency_AVG_colorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_frequency_AVG_colorbar.eps')))

    exec_times_id = [exec_times_id, "Doppler AVG frequency heatmap"];
    exec_times_time = [exec_times_time, toc];

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
    imwrite(rescale(f_RMS), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_frequency_RMS.png')), 'png');

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

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_colorbarRMSFrequency.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, '1_colorbarRMSFrequency.eps')))

    exec_times_id = [exec_times_id, "Doppler RMS frequency heatmap"];
    exec_times_time = [exec_times_time, toc];

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
        graphSignal(ToolBox, '2_signalsRaw', folder, ...
            t, arterial_signal, '-', cArtery, ...
            t, background_signal, ':', cBlack, ...
            t, venous_signal, '-', cVein, ...
            Title = 'Arterial Pulse Waveform and Background Signal', xlabel = strXlabel, ylabel = strYlabel, ...
            Legends = {'Artery', 'Background', 'Vein'}, TxtName = {'FullArterialSignal', 'FullBackgroundSignal', 'FullVenousSignal'})
    else
        graphSignal(ToolBox, '2_signalsRaw', folder, ...
            t, arterial_signal, '-', cArtery, ...
            t, background_signal, ':', cBlack, ...
            Title = 'Arterial Pulse Waveform and Background Signal', xlabel = strXlabel, ylabel = strYlabel, ...
            Legends = {'Artery', 'Background'}, TxtName = {'FullArterialSignal', 'FullBackgroundSignal'})
    end

    exec_times_id = [exec_times_id, "Calculate raw signals"];
    exec_times_time = [exec_times_time, toc];

    %% 3) Smoothing signals

    tic

    delta_f_arterial_signal = arterial_signal - background_signal;
    delta_f_arterial_smooth = smoothdata(delta_f_arterial_signal, 'lowess');

    graphSignal(ToolBox, '3_arterialSignalSmoothed', folder, ...
        t, delta_f_arterial_signal, ':', cArtery, ...
        t, delta_f_arterial_smooth, '-', cArtery, ...
        Title = 'Arterial Signal Smoothing', xlabel = strXlabel, ylabel = strYlabel, ...
        Legends = {'Noisy', 'Robust linear regression'})

    if veinsAnalysis
        delta_f_venous_signal = venous_signal - background_signal;
        delta_f_venous_smooth = smoothdata(delta_f_venous_signal, 'lowess');

        graphSignal(ToolBox, '3_signalsSmoothed', folder, ...
            t, delta_f_arterial_smooth, '-', cArtery, ...
            t, delta_f_venous_smooth, '-.', cVein, ...
            Title = 'Smoothed Arterial Signal and Venous Signal', xlabel = strXlabel, ylabel = strYlabel, ...
            Legends = {'Artery', 'Vein'}, TxtName = {'FullArterialSignalSmoothed', 'FullVenousSignalSmoothed'})

        graphSignal(ToolBox, '3_venousSignalSmoothed', folder, ...
            t, delta_f_venous_signal, ':', cVein, ...
            t, delta_f_venous_smooth, '-', cVein, ...
            Title = 'Smoothed venous signal', xlabel = strXlabel, ylabel = strYlabel, ...
            Legends = {'Noisy', 'Robust linear regression'})
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

    graphSignal(ToolBox, '3_filteredPulseVsResidual', folder, ...
        t, delta_f_arterial_smooth, ':', cArtery, ...
        t, noise, '-', cBlack, ...
        Title = 'signal vs. noise', xlabel = strXlabel, ...
        Legend = {'filtered arterial pulse', 'residual'}, ...
        yLines = [0 std(noise) 2 * std(noise) 3 * std(noise)], yLineLabels = {'', '1 std', '2 std', '3 std'}, ...
        TxtName = {'FilteredArterialPulse', 'ResidualArterialPulse'})

    exec_times_id = [exec_times_id, "Smoothing signals"];
    exec_times_time = [exec_times_time, toc];

    %% 4) Calculation of pulse signal derivative and finding/smoothing pulses

    tic

    fullArterialPulseDerivative = gradient(arterial_signal);
    fullArterialPulseSmoothDerivative = gradient(delta_f_arterial_smooth);

    if veinsAnalysis
        fullVenousPulseDerivative = gradient(venous_signal);
        fullVenousPulseSmoothDerivative = gradient(delta_f_venous_smooth);
    end

    exec_times_id = [exec_times_id, "Calculate pulse derivative"];
    exec_times_time = [exec_times_time, toc];

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

        graphSignal(ToolBox, '4_derivative', folder, ...
            t, fullArterialPulseDerivative, ':', cArtery, ...
            t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
            t, fullVenousPulseDerivative, ':', cVein, ...
            t, fullVenousPulseSmoothDerivative, '-', cVein, ...
            Title = 'derivative of the arterial and venous pulse waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
            Legend = {'\delta <p(t)> - <b(t)>', 'from smoothed data'}, ...
            TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative', 'fullVenousPulseDerivative', 'FullVenousPulseSmoothDerivative'}, ...
            TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts)

    else
        textsX = t(sysIdxList) - 0.3;
        textsY = fullArterialPulseSmoothDerivative(sysIdxList)' + 0.03;
        texts = cell(1, size(sysIdxList));

        for n = 1:numel(sysIdxList)
            texts{n} = num2str(n);
        end

        graphSignal(ToolBox, '4_derivative', folder, ...
            t, fullArterialPulseDerivative, ':', cArtery, ...
            t, fullArterialPulseSmoothDerivative, '-', cArtery, ...
            Title = 'derivative of the arterial waveform', xlabel = strXlabel, ylabel = 'A.U.', ...
            Legend = {'\delta <p(t)> - <b(t)>', 'from smoothed data'}, ...
            TxtName = {'FullArterialPulseDerivative', 'FullArterialPulseSmoothDerivative'}, ...
            TxtFigX = textsX, TxtFigY = textsY, TxtFigString = texts)
    end

end

%% Local BKG Artery and Veins %~1min

tic
exec_times_id = [exec_times_id, "Local BKG Artery and Veins"];

if veinsAnalysis
    maskVesselDilated = imdilate(maskArtery | maskVein, strel('disk', PW_params.local_background_width));
    imwrite(maskVesselDilated - (maskArtery | maskVein), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');

else
    maskVesselDilated = imdilate(maskArtery, strel('disk', PW_params.local_background_width));
    imwrite(maskVesselDilated - (maskArtery), fullfile(ToolBox.PW_path_png, 'mask', sprintf("%s_%s", ToolBox.main_foldername, 'maskVesselDilated.png')), 'png');
end

f_RMS_background = zeros(numX, numY, numFrames, 'single');

parfor frameIdx = 1:numFrames
    f_RMS_background(:, :, frameIdx) = single(regionfill(f_RMS_video(:, :, frameIdx), maskVesselDilated));
end

graphSignal(ToolBox, '5_Arteries_fRMS', folder, ...
    t, squeeze(sum(f_RMS_video .* maskArtery, [1, 2]) / nnz(maskArtery)), '-', cArtery, ...
    t, squeeze(sum(f_RMS_background .* maskArtery, [1, 2]) / nnz(maskArtery)), '--', cBlack, ...
    Title = 'Average f_{RMS} in Arteries', xlabel = strXlabel, ylabel = strYlabel, ...
    Legend = {'Arteries', 'Background'})

if veinsAnalysis
    graphSignal(ToolBox, '5_Veins_fRMS', folder, ...
        t, squeeze(sum(f_RMS_video .* maskVein, [1, 2]) / nnz(maskVein)), '-', cVein, ...
        t, squeeze(sum(f_RMS_background .* maskVein, [1, 2]) / nnz(maskVein)), '--', cBlack, ...
        Title = 'Average f_{RMS} in Veins', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Veins', 'Background'})
end

exec_times_time = [exec_times_time, toc];

%% Difference calculation

tic
exec_times_id = [exec_times_id, "Difference calculation"];

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
        LocalBKG_vesselM2(:, :, frameIdx) = single(regionfill(M2_data_video(:, :, frameIdx), maskVesselDilated));
        LocalBKG_vesselM0(:, :, frameIdx) = single(regionfill(M0_data_video(:, :, frameIdx), maskVesselDilated));
    end

    tmpM2 = M2_data_video - LocalBKG_vesselM2;
    tmpM0 = M0_data_video - LocalBKG_vesselM0;
    tmpM2M0 = tmpM2 ./ mean(tmpM0, [1 2]);
    delta_f_RMS = sign(tmpM2M0) .* sqrt(abs(tmpM2M0));
    clear tmpM2 tmpM0 tmpM2M0 M2_data_video

else % DIFFERENCE LAST

    delta_f_RMS = f_RMS_video - f_RMS_background;

end

v_RMS_video = ToolBox.ScalingFactorVelocityInPlane * delta_f_RMS;

if veinsAnalysis
    graphSignal(ToolBox, '6_Vessels_velocity', folder, ...
        t, squeeze(sum(v_RMS_video .* maskArtery, [1, 2]) / nnz(maskArtery)), '-', cArtery, ...
        t, squeeze(sum(v_RMS_video .* maskVein, [1, 2]) / nnz(maskVein)), '-', cVein, ...
        Title = 'Average estimated velocity in Arteries and Veins', xlabel = strXlabel, ylabel = 'mm/s')
else
    graphSignal(ToolBox, '6_Arteries_velocity', folder, ...
        t, squeeze(sum(v_RMS_video .* maskArtery, [1, 2]) / nnz(maskArtery)), '-', cArtery, ...
        Title = 'Average estimated velocity in Arteries', xlabel = strXlabel, ylabel = 'mm/s')
end

exec_times_time = [exec_times_time, toc];

f18 = figure(18);
f18.Position = [1100 485 350 420];
LocalBackground_in_vessels = mean(f_RMS_background, 3) .* maskVesselDilated + ones(numX, numY) * mean(f_RMS_background, 'all') .* ~maskVesselDilated;
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

if exportVideos
    timePeriod = ToolBox.stride / ToolBox.fs / 1000;

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.avi')));
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_vessels.mp4')), 'MPEG-4');

    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_video), fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'f_AVG_vessels.avi')));
    parfeval(backgroundPool, @writeVideoOnDisc, 0, mat2gray(f_RMS_video), fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'f_AVG_vessels.mp4')), 'MPEG-4');

    parfeval(backgroundPool, @writeGifOnDisc, 0, mat2gray(f_RMS_background), fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.main_foldername, 'LocalBackground_vessels')), timePeriod);
    parfeval(backgroundPool, @writeGifOnDisc, 0, mat2gray(f_RMS_video), fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.main_foldername, 'f_AVG_vessels')), timePeriod);
end

clear LocalBackground_in_vessels f_RMS_background

if entirePulseAnalysis
    %% Creation of the avg pulse for In-plane arteries ~5min

    tic

    fprintf("Average Pulse\n")
    [onePulseVideo, ~, ~, onePulseVideoM0] = createOneCycle(f_RMS_video, M0_data_video, maskArtery, sysIdxList, numFramesInterp, path, ToolBox);

    clear f_RMS_video f_RMS_video

    fprintf("Average Pulse minus Background\n")
    [onePulseVideominusBKG, selectedPulseIdx, cycles_signal, ~] = createOneCycle(delta_f_RMS, M0_data_video, maskArtery, sysIdxList, numFramesInterp, path, ToolBox);

    clear delta_f_RMS

    avgArterialPulseHz = squeeze(sum(onePulseVideominusBKG .* maskArtery, [1 2])) / nnz(maskArtery);
    avgArterialPulseVelocityInPlane = avgArterialPulseHz * ToolBox.ScalingFactorVelocityInPlane;

    v_OneCycle = (onePulseVideominusBKG .* maskArtery + onePulseVideo .* ~maskArtery) * ToolBox.ScalingFactorVelocityInPlane;
    ArterialResistivityIndex(v_OneCycle, M0_ff_video, maskArtery, ToolBox, path)

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

    %% Arterial pulse wave analysis

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
        plot(t(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'k-', LineWidth = 2);
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
                plot(t, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 2);
                hold on
            else
                plot(t, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
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

end

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
