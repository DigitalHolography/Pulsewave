function [v_RMS_one_cycle, v_RMS_all, flowVideoRGB, exec_times, total_time] = pulseAnalysis(Ninterp, fullVideoM2M0, fullVideoM1M0, fullVideoM2, ~, fullVideoM0, sys_index_list, meanIm, maskArtery, maskVein, maskBackground, ToolBox, path)

    % Variable : LocalBKG_artery, Taille : 10631287200 bytes
    % Variable : fullVideoM1M0, Taille : 10631287200 bytes (DEBUT)
    % Variable : fullVideoM2M0, Taille : 10631287200 bytes (DEBUT)
    % Variable : maskArtery, Taille : 18849800 bytes (DEBUT)
    % Variable : meanIm, Taille : 18849800 bytes  (DEBUT)
    % Variable : maskBackground, Taille : 2356225 bytes (DEBUT)
    % Variable : maskVein, Taille : 2356225 bytes (DEBUT)
    % Variable : variableInfo, Taille : 12898 bytes

    exec_times_id = [];
    exec_times_time = [];
    total_time = 0;

    PW_params = Parameters_json(path);

    veins_analysis = PW_params.veins_analysis;

    N_frame = size(fullVideoM2M0, 3);

    meanIm = mat2gray(meanIm);
    tolVal = [0.02, 0.98];
    meanIm = mat2gray(imadjust(meanIm, stretchlim(meanIm, tolVal)));

    %% Display images and figures
    strXlabel = 'Time(s)'; %createXlabelTime(1);

    strYlabel = 'frequency (kHz)';
    range0(1:2) = clim;
    fullTime = linspace(0, N_frame * ToolBox.stride / ToolBox.fs / 1000, N_frame);

    %% Doppler AVG frequency heatmap

    tic

    heatmap_AVG_raw = squeeze(mean(fullVideoM1M0, 3));

    clear fullVideoM1M0
    mkdir(ToolBox.PW_path_png, 'pulseAnalysis')
    mkdir(ToolBox.PW_path_eps, 'pulseAnalysis')
    imwrite(heatmap_AVG_raw, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'heatmap_AVG_raw.png')), 'png');

    %  Doppler AVG frequency heatmap
    figure(78)
    imagesc(heatmap_AVG_raw);
    colormap gray
    title('AVG frequency map RAW');
    fontsize(gca, 12, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    range(1:2) = clim;

    clear heatmap_AVG_raw

    exec_times_id = [exec_times_id, "Doppler AVG frequency heatmap"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    %% calculate raw signals of arteries, background and veins

    tic

    fullBackgroundSignal = fullVideoM2M0 .* maskBackground;
    fullBackgroundSignal = squeeze(sum(fullBackgroundSignal, [1 2])) / nnz(maskBackground);

    fullArterialPulse = fullVideoM2M0 .* maskArtery;
    fullArterialPulse = squeeze(sum(fullArterialPulse, [1 2])) / nnz(maskArtery);
    fullArterialPulseMinusBackground = fullArterialPulse - fullBackgroundSignal;

    % fullBackgroundM2 = fullVideoM2 .* maskBackground;
    % fullBackgroundM2 = squeeze(sum(fullBackgroundM2, [1 2])) / nnz(maskBackground);

    % fullArterialPulseM2 = fullVideoM2 .* maskArtery;
    % fullArterialPulseM2 = squeeze(sum(fullArterialPulseM2, [1 2])) / nnz(maskArtery);
    % fullArterialPulseMinusBackgroundM2 = fullArterialPulseM2 - fullBackgroundM2;

    % fullBackgroundM0 = fullVideoM0 .* maskBackground;
    % fullBackgroundM0 = squeeze(sum(fullBackgroundM0, [1 2])) / nnz(maskBackground);

    % fullArterialPulseM0 = fullVideoM0 .* maskArtery;
    % fullArterialPulseM0 = squeeze(sum(fullArterialPulseM0, [1 2])) / nnz(maskArtery);
    % fullArterialPulseMinusBackgroundM0 = fullArterialPulseM0 - fullBackgroundM0;

    if veins_analysis
        fullVenousSignal = fullVideoM2M0 .* maskVein;
        fullVenousSignal = squeeze(sum(fullVenousSignal, [1 2])) / nnz(maskVein);
        fullVenousSignalMinusBackground = fullVenousSignal - fullBackgroundSignal;

        % fullVenousSignalM2 = fullVideoM2 .* maskVein;
        % fullVenousSignalM2 = squeeze(sum(fullVenousSignalM2, [1 2])) / nnz(maskVein);
        % fullVenousSignalMinusBackgroundM2 = fullVenousSignalM2 - fullBackgroundM2;

        % fullVenousSignalM0 = fullVideoM0 .* maskVein;
        % fullVenousSignalM0 = squeeze(sum(fullVenousSignalM0, [1 2])) / nnz(maskVein);
        % fullVenousSignalMinusBackgroundM0 = fullVenousSignalM0 - fullBackgroundM0;
    end

    exec_times_id = [exec_times_id, "Calculate raw signals"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    %% cleaning signals

    tic

    % remove outliers
    % 1st pass
    disp('remove outliers... 1st pass.');
    [idxOutPw, fullArterialPulseRmOut] = discardPulseWaveOutliers(fullArterialPulseMinusBackground, 3);
    [idxOutBkg, ~] = discardPulseWaveOutliers(fullBackgroundSignal, 3);

    if veins_analysis
        [~, fullVenousSignalRmOut] = discardPulseWaveOutliers(fullVenousSignalMinusBackground, 3);
    end

    % FIXME

    dataReliabilityIndex1 = ceil(100 * (1 - PW_params.pulseAnal_dataReliabilityFactor * (length(idxOutPw) / length(fullArterialPulse) + length(idxOutBkg) / length(fullBackgroundSignal))));
    disp(['data reliability index 1 : ' num2str(dataReliabilityIndex1) ' %']);

    % smooth trendline data by iterative local linear regression.

    fullArterialPulseClean = smoothdata(fullArterialPulseRmOut, 'lowess');

    if veins_analysis
        fullVenousSignalClean = smoothdata(fullVenousSignalRmOut, 'lowess');
    end

    % FIXME: usecase?
    % fullBackgroundSignalClean = smoothdata(fullBackgroundSignalRmOut,'lowess');

    exec_times_id = [exec_times_id, "Cleaning signals"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    % figure(2)
    % plot(fullBackgroundSignal)
    % hold on
    % plot(fullArterialPulse)
    % hold on
    % plot(fullBackgroundSignalClean)
    % hold off
    % legend

    %% PLOTS

    figure(30)
    plot(fullTime, fullArterialPulseMinusBackground, '-k', fullTime, fullArterialPulseClean, '-k', 'LineWidth', 2);
    %title('arterial pulse minus background vs. filtered pulse');
    %legend('<p(t)> - <b(t)>', 'local linear regression');
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'filteredPulse.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'filteredPulse.eps')))

    plot2txt(fullTime, fullArterialPulseMinusBackground, 'FullArterialPulseMinusBackground', ToolBox)
    plot2txt(fullTime, fullArterialPulseClean, 'FullArterialPulseClean', ToolBox)

    if veins_analysis
        figure(31)
        plot(fullTime, fullVenousSignalMinusBackground, ':k', ...
            fullTime, fullVenousSignalClean, '-k', ...
            'LineWidth', 2);
        title('venous signal minus background vs. filtered pulse');
        legend('<p(t)> - <b(t)>', 'local linear regression');
        fontsize(gca, 12, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;
    end

    %% calculate the pulse derivative and finding/cleaning pulses

    tic

    fullArterialPulseDerivative = diff(fullArterialPulse);
    fullArterialPulseCleanDerivative = diff(fullArterialPulseClean);

    % now cleanup dataCube to create_one_cycle()
    % strategy : use average pulse profiles to detect and
    % find  noisy frames @ >3 std from zero-mean
    % replace them with cleaned data
    noise = sqrt(abs(abs(fullArterialPulseMinusBackground) .^ 2 - abs(fullArterialPulseClean) .^ 2));
    idxOutNoise = find(noise > PW_params.pulseAnal_outNoiseThreshold * std(noise));

    dataReliabilityIndex2 = ceil(100 * (1 - (length(idxOutNoise) / length(fullArterialPulseClean))));
    disp(['data reliability index 2 : ' num2str(dataReliabilityIndex2) ' %']);

    exec_times_id = [exec_times_id, "Calculate pulse derivative"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    %% PLOTS

    figure(40)
    plot( ...
        fullTime(1:length(fullArterialPulseDerivative)), fullArterialPulseDerivative, ':k', ...
        fullTime(1:length(fullArterialPulseCleanDerivative)), fullArterialPulseCleanDerivative, '-k', ...
        'LineWidth', 2);
    text(fullTime(sys_index_list) - 0.3, fullArterialPulseCleanDerivative(sys_index_list) + 0.03, num2str((1:numel(sys_index_list))'))
    title('derivative of the arterial pulse waveform');
    legend('\delta <p(t)> - <b(t)>', 'from smoothed data');
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel('A.U.', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'derivative.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'derivative.eps')))

    plot2txt(fullTime(1:length(fullArterialPulseDerivative)), fullArterialPulseDerivative, 'FullArterialPulseDerivative', ToolBox)
    plot2txt(fullTime(1:length(fullArterialPulseCleanDerivative)), fullArterialPulseCleanDerivative, 'FullArterialPulseCleanDerivative', ToolBox)

    clear fullArterialPulseDerivative fullArterialPulseCleanDerivative

    % c = colorbar('southoutside');
    %
    % Colorbar for raw/flattened image
    colorfig = figure(2410);
    colorfig.Units = 'normalized';
    colormap(c);
    colormap gray
    clim(range0);
    hCB = colorbar('north');
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 15, "points");
    colorTitleHandle = get(hCB, 'Title');
    titleString = 'RMS Doppler frequency (kHz)';
    set(colorTitleHandle, 'String', titleString);

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'RMS_frequency_colorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'RMS_frequency_colorbar.eps')))

    figure(44)
    plot( ...
        fullTime(1:length(fullArterialPulseClean)), fullArterialPulseClean, ':k', ...
        fullTime(1:length(noise)), noise, '-k', ...
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

    plot2txt(fullTime(1:length(fullArterialPulseClean)), fullArterialPulseClean, 'FilteredArterialPulse', ToolBox)
    plot2txt(fullTime(1:length(noise)), noise, 'ResidualArterialPulse', ToolBox)

    % % C'EST QUOI?
    % figure(85)
    % clim([min(range),max(range)]);

    %% Local BKG Artery and Veins %~1min

    tic

    local_mask_artery = imdilate(maskArtery, strel('disk', PW_params.local_background_width));
    % local_mask_artery = and(local_mask_artery, not(maskVein));
    LocalBKG_artery = zeros(size(fullVideoM2M0));
    LocalBKG_arteryM2 = zeros(size(fullVideoM2));
    LocalBKG_arteryM0 = zeros(size(fullVideoM0));

    if veins_analysis
        local_mask_vein = imdilate(maskVein, strel('disk', PW_params.local_background_width));
        % local_mask_vein = and(local_mask_vein, not(maskArtery));
        LocalBKG_vein = zeros(size(fullVideoM2M0));
        LocalBKG_veinM2 = zeros(size(fullVideoM2));
        LocalBKG_veinM0 = zeros(size(fullVideoM0));
    end

    if veins_analysis

        parfor nn = 1:N_frame
            LocalBKG_artery(:, :, nn) = single(regionfill(fullVideoM2M0(:, :, nn), local_mask_artery));
            LocalBKG_vein(:, :, nn) = single(regionfill(fullVideoM2M0(:, :, nn), local_mask_vein));
            LocalBKG_arteryM2(:, :, nn) = single(regionfill(fullVideoM2(:, :, nn), local_mask_artery));
            LocalBKG_veinM2(:, :, nn) = single(regionfill(fullVideoM2(:, :, nn), local_mask_vein));
            LocalBKG_arteryM0(:, :, nn) = single(regionfill(fullVideoM0(:, :, nn), local_mask_artery));
            LocalBKG_veinM0(:, :, nn) = single(regionfill(fullVideoM0(:, :, nn), local_mask_vein));
        end

    else

        parfor nn = 1:N_frame
            LocalBKG_artery(:, :, nn) = single(regionfill(fullVideoM2M0(:, :, nn), local_mask_artery));
            LocalBKG_arteryM2(:, :, nn) = single(regionfill(fullVideoM2(:, :, nn), local_mask_artery));
            LocalBKG_arteryM0(:, :, nn) = single(regionfill(fullVideoM0(:, :, nn), local_mask_artery));
        end

    end

    if veins_analysis
        % Remove arteries (always calculated) from local BKG
        if PW_params.DiffFirstCalculationsFlag == 0 %SIGNED DIFFERENCE FIRST

            tmp = fullVideoM2M0 .^ 2 - (LocalBKG_vein .* ~local_mask_artery + LocalBKG_artery .* local_mask_artery) .^ 2;
            fullVideoM2M0minusBKG = sign(tmp) .* sqrt(abs(tmp));
            clear tmp

        elseif PW_params.DiffFirstCalculationsFlag == 1 % DIFFERENCE FIRST

            tmp = fullVideoM2M0 .^ 2 - (LocalBKG_vein .* ~local_mask_artery + LocalBKG_artery .* local_mask_artery) .^ 2;
            tmp = tmp .* (tmp > 0);
            fullVideoM2M0minusBKG = sqrt(abs(tmp));
            clear tmp

        elseif PW_params.DiffFirstCalculationsFlag == 2 % TO BE TESTED

            fullVideoM2eff = fullVideoM2 - (LocalBKG_veinM2 .* ~local_mask_artery + LocalBKG_arteryM2 .* local_mask_artery);
            fullVideoM0eff = fullVideoM0 - (LocalBKG_veinM0 .* ~local_mask_artery + LocalBKG_arteryM0 .* local_mask_artery);
            fullVideoM2M0eff = fullVideoM2eff ./ fullVideoM0eff;
            fullVideoM2M0minusBKG = sqrt(fullVideoM2M0eff .* (fullVideoM2M0eff > 0));

        else % DIFFERENCE LAST

            fullVideoM2M0minusBKG = fullVideoM2M0 - (LocalBKG_vein .* ~local_mask_artery + LocalBKG_artery .* local_mask_artery);

        end

    else

        if PW_params.DiffFirstCalculationsFlag == 0 %SIGNED DIFFERENCE FIRST

            tmp = fullVideoM2M0 .^ 2 - LocalBKG_artery .^ 2;
            fullVideoM2M0minusBKG = sign(tmp) .* sqrt(abs(tmp));
            clear tmp

        elseif PW_params.DiffFirstCalculationsFlag == 1 % DIFFERENCE FIRST

            tmp = fullVideoM2M0 .^ 2 - LocalBKG_artery .^ 2;
            tmp = tmp .* (tmp > 0);
            fullVideoM2M0minusBKG = sqrt(tmp);
            clear tmp

        elseif PW_params.DiffFirstCalculationsFlag == 2 % TO BE TESTED

            fullVideoM2eff = fullVideoM2 - LocalBKG_arteryM2;
            fullVideoM0eff = fullVideoM0 - LocalBKG_arteryM0;
            fullVideoM2M0eff = fullVideoM2eff ./ fullVideoM0eff;
            fullVideoM2M0minusBKG = sqrt(fullVideoM2M0eff .* (fullVideoM2M0eff > 0));

        else % DIFFERENCE LAST

            fullVideoM2M0minusBKG = fullVideoM2M0 - LocalBKG_artery;

        end

    end

    exec_times_id = [exec_times_id, "Local Backgrounds"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    %% FIXME : compute true regularied cube by replacing bad frames

    fullArterialPulseRegularized = squeeze(sum(fullVideoM2M0minusBKG .* maskArtery, [1 2])) / nnz(maskArtery);

    %% PLOT

    figure(43)
    plot( ...
        fullTime(1:length(fullArterialPulseRegularized)), fullArterialPulseRegularized, '-k', ...
        fullTime(1:length(fullArterialPulse)), fullArterialPulseMinusBackground, ':k', ...
        'LineWidth', 2);
    yline(0, ':', {''}, LineWidth = 2);
    title('Regularized pulse wave signal averaged in retinal arteries');
    legend('Regularized arterial pulse', '<p(t)> - <b(t)>');
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'regularizedPulse.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'regularizedPulse.eps')))

    plot2txt(fullTime(1:length(fullArterialPulseRegularized)), fullArterialPulseRegularized, 'FullArterialPulseRegularized', ToolBox)
    plot2txt(fullTime(1:length(fullArterialPulse)), fullArterialPulseMinusBackground, 'FullArterialPulseMinusBackground', ToolBox)

    f18 = figure(18);
    f18.Position = [1100 485 350 420];
    LocalBackground_in_arteries = mean(LocalBKG_artery, 3) .* local_mask_artery + ones(size(LocalBKG_artery, 1)) * mean(LocalBKG_artery, 'all') .* ~local_mask_artery;
    imagesc(LocalBackground_in_arteries);
    colormap gray
    title('Local Background in arteries');
    fontsize(gca, 14, "points");
    set(gca, 'LineWidth', 2);
    c = colorbar('southoutside');
    c.Label.String = 'RMS Doppler frequency (kHz)';
    c.Label.FontSize = 12;
    axis off
    axis image
    imwrite(rescale(LocalBackground_in_arteries), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_in_arteries.png')))

    range(1:2) = clim;

    clear local_mask_artery

    if veins_analysis
        f19 = figure(19);
        f19.Position = [1100 485 350 420];
        LocalBackground_in_veins = mean(LocalBKG_vein, 3) .* local_mask_vein + ones(size(LocalBKG_vein, 1)) * mean(LocalBKG_vein, 'all') .* ~local_mask_vein;
        imagesc(LocalBackground_in_veins);
        colormap gray
        title('Local Background in vein');
        fontsize(gca, 14, "points");
        set(gca, 'LineWidth', 2);
        c = colorbar('southoutside');
        c.Label.String = 'RMS Doppler frequency (kHz)';
        c.Label.FontSize = 12;
        axis off
        axis image
        imwrite(rescale(LocalBackground_in_veins), fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'LocalBackground_in_veins.png')));

        range(1:2) = clim;
    end

    if veins_analysis
        figure(20)
        plot(fullTime, fullArterialPulse, '-k', fullTime, fullBackgroundSignal, ':k', fullTime, fullVenousSignal, '-.k', 'LineWidth', 2);
        title('Arterial Pulse Waveform and Background Signal'); % averaged outside of segmented vessels
        legend('Arterial Pulse', 'Background', 'Venous Signal');
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVsBackground.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVsBackground.eps')))

        plot2txt(fullTime, fullArterialPulse, 'FullArterialPulse', ToolBox)
        plot2txt(fullTime, fullBackgroundSignal, 'FullBackgroundSignal', ToolBox)
        plot2txt(fullTime, fullVenousSignal, 'FullVenousSignal', ToolBox)

    else
        figure(20)
        plot(fullTime, fullArterialPulse, '-k', fullTime, fullBackgroundSignal, ':k', 'LineWidth', 2);
        title('Arterial Pulse Waveform and Background Signal'); % averaged outside of segmented vessels
        legend('Arterial Pulse', 'Background');
        fontsize(gca, 14, "points");
        xlabel(strXlabel, 'FontSize', 14);
        ylabel(strYlabel, 'FontSize', 14);
        pbaspect([1.618 1 1]);
        set(gca, 'LineWidth', 2);
        axis tight;
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVsBackground.png')))
        exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'pulseVsBackground.eps')))

        plot2txt(fullTime, fullArterialPulse, 'FullArterialPulse', ToolBox)
        plot2txt(fullTime, fullBackgroundSignal, 'FullBackgroundSignal', ToolBox)

    end

    clear LocalBKG_artery LocalBKG_vein

    % v_RMS_all = ToolBox.ScalingFactorVelocityInPlane*fullVideoM2M0minusBKG;

    %% Global BKG for beautiful images

    A = ones(size(fullVideoM2M0));

    for pp = 1:N_frame
        A(:, :, pp) = A(:, :, pp) * fullBackgroundSignal(pp);
    end

    fullVideoM2M0 = fullVideoM2M0 - A;

    clear fullBackgroundSignal fullVenousSignal A

    %% Construct Velocity video ~3min

    tic

    flowVideoRGB = zeros(size(fullVideoM2M0, 1), size(fullVideoM2M0, 2), 3, size(fullVideoM2M0, 3));
    fullVideoM0_norm = rescale(fullVideoM0);

    if veins_analysis

        v_mean = mat2gray(squeeze(mean(fullVideoM2M0(:, :, :), 3)));
        fullVideoM0_normMean = mean(fullVideoM0_norm, 3);
        [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(v_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        [hue_vein_mean, sat_vein_mean, val_vein_mean, ~] = createHSVmap(v_mean, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
        val_mean = v_mean .* (~(maskArtery + maskVein)) + val_artery_mean .* maskArtery + val_vein_mean .* maskVein;
        flowVideoRGB_mean = hsv2rgb(hue_artery_mean + hue_vein_mean, sat_artery_mean + sat_vein_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* (maskArtery + maskVein) + ones(size(flowVideoRGB_mean)) .* fullVideoM0_normMean .* ~(maskArtery + maskVein);
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png')))

        for ii = 1:size(fullVideoM2M0, 3)
            v = mat2gray(squeeze(fullVideoM2M0(:, :, ii)));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            [hue_vein, sat_vein, val_vein, ~] = createHSVmap(v, maskVein, 0.68, 0.5); %0.5/0.68 for cyan-dark blue range
            val = v .* (~(maskArtery + maskVein)) + val_artery .* maskArtery + val_vein .* maskVein;
            flowVideoRGB(:, :, :, ii) = hsv2rgb(hue_artery + hue_vein, sat_artery + sat_vein, val);
            flowVideoRGB(:, :, :, ii) = flowVideoRGB(:, :, :, ii) .* (maskArtery + maskVein) + ones(size(flowVideoRGB(:, :, :, ii))) .* fullVideoM0_norm(:, :, ii) .* ~(maskArtery + maskVein);
        end

    else

        v_mean = mat2gray(squeeze(mean(fullVideoM2M0(:, :, :), 3)));
        fullVideoM0_normMean = mean(fullVideoM0_norm, 3);
        [hue_artery_mean, sat_artery_mean, val_artery_mean, ~] = createHSVmap(v_mean, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
        val_mean = v_mean .* ~maskArtery + val_artery_mean .* maskArtery;
        flowVideoRGB_mean = hsv2rgb(hue_artery_mean, sat_artery_mean, val_mean);
        flowVideoRGB_mean = flowVideoRGB_mean .* maskArtery + ones(size(flowVideoRGB_mean)) .* fullVideoM0_normMean .* ~maskArtery;
        imwrite(flowVideoRGB_mean, fullfile(ToolBox.PW_path_png, sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideo.png')))

        for ii = 1:size(fullVideoM2M0, 3)
            v = mat2gray(squeeze(fullVideoM2M0(:, :, ii)));
            [hue_artery, sat_artery, val_artery, ~] = createHSVmap(v, maskArtery, 0, 0.18); % 0 / 0.18 for orange-yellow range
            val = v .* ~maskArtery + val_artery .* maskArtery;
            flowVideoRGB(:, :, :, ii) = hsv2rgb(hue_artery, sat_artery, val);
            flowVideoRGB(:, :, :, ii) = flowVideoRGB(:, :, :, ii) .* (maskArtery) + ones(size(flowVideoRGB(:, :, :, ii))) .* fullVideoM0_norm(:, :, ii) .* ~(maskArtery);
        end

    end

    clear hue_artery sat_artery val_artery hue_vein sat_vein val_vein

    % save video
    % avi

    parfeval(backgroundPool,@writeVideoOnDisc,0,flowVideoRGB,fullfile(ToolBox.PW_path_avi, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')));

    timePeriod = ToolBox.stride / ToolBox.fs / 1000;

    writeGifOnDisc(flowVideoRGB,fullfile(ToolBox.PW_path_gif, sprintf("%s_%s.gif", ToolBox.PW_folder_name, "flowMap")),timePeriod);


    % f_RMS_artery = sum(flowVideoRGB .* maskArtery, [1 2]) / nnz(maskArtery);
    % f_RMS_vein = sum(flowVideoRGB .* maskVein, [1 2]) / nnz(maskVein);
    % f_RMS_max_Arteries = max(f_RMS_artery(:));
    % f_RMS_max_Veins = max(f_RMS_vein(:));
    % f_RMS_min_Arteries = min(f_RMS_artery(:));
    % f_RMS_min_Veins = min(f_RMS_vein(:));

    % mp4
    parfeval(backgroundPool,@writeVideoOnDisc,0,flowVideoRGB,fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'flowVideo')),'MPEG-4');


    exec_times_id = [exec_times_id, "Construct velocity video"];
    exec_times_time = [exec_times_time, toc];
    total_time = total_time + toc;

    %% Creation of the avg pulse for In-plane arteries ~5min

    tic

    [onePulseVideo, ~, ~, onePulseVideoM0] = create_one_cycle(fullVideoM2M0, fullVideoM0, maskArtery, sys_index_list, Ninterp, path, ToolBox);

    clear fullVideoM2M0

    [onePulseVideominusBKG, selectedPulseIdx, cycles_signal, ~] = create_one_cycle(fullVideoM2M0minusBKG, fullVideoM0, maskArtery, sys_index_list, Ninterp, path, ToolBox);

    %% Velocity Calculation

    v_RMS_all = ToolBox.ScalingFactorVelocityInPlane * fullVideoM2M0minusBKG * ToolBox.NormalizationFactor;
    % fprintf('Factor Normalization was performed: %f\n',ToolBox.NormalizationFactor)
    % v_RMS_all = ToolBox.ScalingFactorVelocityInPlane *
    % fullVideoM2M0minusBKG + ToolBox.NormalizationOffset;
    % fprintf('Offset Normalization was performed: +%d µL/min',ToolBox.NormalizationOffset)

    clear fullVideoM2M0minusBKG

    avgArterialPulseHz = squeeze(sum(onePulseVideominusBKG .* maskArtery, [1 2])) / nnz(maskArtery);
    avgArterialPulseVelocityInPlane = avgArterialPulseHz * ToolBox.ScalingFactorVelocityInPlane;

    v_RMS_one_cycle = (onePulseVideominusBKG .* maskArtery + onePulseVideo .* ~maskArtery) * ToolBox.ScalingFactorVelocityInPlane;

    % avi
    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(onePulseVideo),fullfile(ToolBox.PW_path_avi,sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.avi')));
    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(onePulseVideoM0),fullfile(ToolBox.PW_path_avi,sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.avi')));


    % mp4
    
    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(onePulseVideo),fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycle.mp4')), 'MPEG-4');
    parfeval(backgroundPool,@writeVideoOnDisc,0,mat2gray(onePulseVideoM0),fullfile(ToolBox.PW_path_mp4, sprintf("%s_%s", ToolBox.main_foldername, 'one_cycleM0.mp4')), 'MPEG-4');

    

    %FIXME: M1/M0 and M2/M0 are subject to aliases at 67 kHz

    blur_time_sys = ceil(Ninterp / PW_params.pulseAnal_blurScaleFactor);
    blur_time_dia = ceil(Ninterp / PW_params.pulseAnal_blurScaleFactor);

    average_cycle_length = 0;
    nb_of_averaged_cycles = 0;

    if size(sys_index_list, 2) == 1
        average_cycle_length = Ninterp;
        nb_of_averaged_cycles = 1;
    else

        for ii = 2:size(sys_index_list, 2)
            average_cycle_length = average_cycle_length + (sys_index_list(ii) - sys_index_list(ii - 1));
            nb_of_averaged_cycles = nb_of_averaged_cycles + 1;
        end

        average_cycle_length = average_cycle_length / (length(sys_index_list) - 1);
    end

    T = linspace(0, (ToolBox.stride / ToolBox.fs * average_cycle_length) / 1000, Ninterp); % /1000 to get time in s

    % save average arterial pulse wave velocity to txt file
    tmp = [T(1:length(avgArterialPulseHz))', avgArterialPulseVelocityInPlane];
    %size(tmp)
    fileID = fopen(fullfile(ToolBox.PW_path_txt, sprintf("%s_%s", ToolBox.main_foldername, 'avgPulse.txt')), 'w');
    fprintf(fileID, '%f %f \r\n', tmp');
    fclose(fileID);

    figure(33)

    for ii = 1:size(cycles_signal, 1)

        if ismember(ii, selectedPulseIdx)
            plot(T, movavgvar(cycles_signal(ii, :), 5), 'k-', 'LineWidth', 1);
            hold on
        else
            plot(T, movavgvar(cycles_signal(ii, :), 5), 'k--', 'LineWidth', 1);
            hold on
        end

    end

    hold on
    plot(T, movavgvar(squeeze(mean(cycles_signal(:, :), 1)), 5), 'k-', 'LineWidth', 2)
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

    %% diastolic Doppler frequency heatmap : 10% of frames before minimum of diastole

    heatmap_dia_raw = squeeze(mean(onePulseVideominusBKG(:, :, floor(0.9 * Ninterp):Ninterp), 3));
    % onePulseVideo2 : no background correction
    % heatmap_dia = squeeze(mean(onePulseVideo2(:,:,floor(0.9*Ninterp):Ninterp),3));
    % heatmap_dia = flat_field_correction(heatmap_dia, ceil(PW_params.flatField_gwRatio*size(heatmap_dia,1)), PW_params.flatField_borderDMap);
    heatmap_dia = flat_field_correction(heatmap_dia_raw, ceil(PW_params.flatField_gwRatio * size(heatmap_dia_raw, 1)), PW_params.flatField_border);
    % heatmap_dia = imflatfield(heatmap_dia,PW_params.flatField_gwRatio*size(heatmap_dia,1)/2);

    clear heatmap_sys_raw

    %% systolic Doppler frequency heatmap : 10% of frames around peak systole

    a = max(ceil(idx_sys - 0.05 * Ninterp), 1);
    b = min(ceil(idx_sys + 0.05 * Ninterp), Ninterp);
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

    % Colorbar for AVG image
    colorfig = figure(1230);
    clim(range0);
    colorfig.Units = 'normalized';
    colormap(c);
    colormap gray
    hCB = colorbar('north');
    set(gca, 'Visible', false)
    set(gca, 'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca, 15, "points");
    colorTitleHandle = get(hCB, 'Title');
    titleString = 'AVG Doppler frequency (kHz)';
    set(colorTitleHandle, 'String', titleString);

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVG_frequency_colorbar.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'AVG_frequency_colorbar.eps')))

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
    diff_avgPulse = diff(movavgvar(avgArterialPulseHz, blur_time_sys)) * ToolBox.ScalingFactorVelocityInPlane * Ninterp / T(end);
    delta = max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:));
    thr = max(pulse_arteries_blurred_dia(:)) - delta ./ exp(1);
    idx_list_threshold_dia = find(pulse_arteries_blurred_dia(:) < thr);
    [max_diff_pulse, idx_T_diff_max] = max(diff_pulse_sys); %faire ca mieux

    if idx_T_diff_max == 1
        acc = abs(max_diff_pulse / (T(idx_T_diff_max) - T(idx_T_diff_max + 1)));
    else
        acc = max_diff_pulse / (T(idx_T_diff_max) - T(idx_T_diff_max - 1));
    end

    % computation of average arterial pulse wave parameters
    T_syst = T(idx_sys);
    systole_area = sum(avgArterialPulseHz(1:idx_sys));
    diastole_area = sum(avgArterialPulseHz(idx_sys:end));
    tmp = systole_area;
    systole_area = systole_area / (diastole_area + systole_area);
    diastole_area = diastole_area / (diastole_area + tmp);
    nb_of_detected_systoles = size(sys_index_list, 2);
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
        T(idx_T_diff_max + 1), ...
        T_syst, ...
        T(idx_sys + idx_list_threshold_dia(1)), ...
        nb_of_detected_systoles, ...
        nb_of_averaged_cycles, ...
        systole_area, ...
        diastole_area);
    fclose(fileID);

    %
    figure(100)
    plot( ...
        T(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'k-', ...
        LineWidth = 2);
    xline(T(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(T(idx_T_diff_max + 1), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (2)', 'FontSize', 14);
    xline(T(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(T(idx_sys + idx_list_threshold_dia(1)), min(avgArterialPulseVelocityInPlane(:)) + 0.1 * (max(avgArterialPulseVelocityInPlane(:)) - min(avgArterialPulseVelocityInPlane(:))), '  (3)', 'FontSize', 14);
    fontsize(gca, 12, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel('blood flow velocity (mm/s)', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;

    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'pulseAnalysis', sprintf("%s_%s", ToolBox.main_foldername, 'avgPulseWave.eps')))

    plot2txt(T(1:length(avgArterialPulseHz)), avgArterialPulseVelocityInPlane, 'AvgArterialPulseVelocityInPlane', ToolBox)

    figure(101)

    for ii = 1:size(cycles_signal, 1)

        if ismember(ii, selectedPulseIdx)
            plot( ...
                T, movavgvar(cycles_signal(ii, :), 5), 'k-', ...
                'LineWidth', 2);
            hold on
        else
            plot( ...
                T, movavgvar(cycles_signal(ii, :), 5), 'k--', ...
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
    plot(T, avgArterialPulseHz, 'k.', ...
        T(1:idx_sys), pulse_arteries_blurred_sys(1:idx_sys), 'k-', ...
        T(idx_sys:Ninterp), pulse_arteries_blurred_dia(1:(Ninterp - idx_sys + 1)), 'k-', ...
        LineWidth = 2);
    xline(T(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(T(idx_T_diff_max + 1), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (2)', 'FontSize', 14);
    xline(T(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(T(idx_sys + idx_list_threshold_dia(1)), min(pulse_arteries_blurred_dia(:)) + 0.1 * (max(pulse_arteries_blurred_dia(:)) - min(pulse_arteries_blurred_dia(:))), '  (3)', 'FontSize', 14);
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

    plot2txt(T, avgArterialPulseHz, 'AvgArterialPulseHz', ToolBox)

    figure(103)
    plot(T(1:end - 1), diff_avgPulse, 'k-', LineWidth = 2);
    x = 0;
    yline(x, ':', LineWidth = 2);
    xline(T(idx_T_diff_max + 1), ':', {}, LineWidth = 2)
    text(T(idx_T_diff_max + 1), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), ' (1)', 'FontSize', 14); %Display at minimum+x %
    xline(T_syst, ':', {}, LineWidth = 2);
    text(T_syst, min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (2)', 'FontSize', 14);
    xline(T(idx_sys + idx_list_threshold_dia(1)), ':', {}, LineWidth = 2);
    text(T(idx_sys + idx_list_threshold_dia(1)), min(diff_avgPulse(:)) + 0.1 * (max(diff_avgPulse(:)) - min(diff_avgPulse(:))), '  (3)', 'FontSize', 14);
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

    plot2txt(T(1:end - 1), diff_avgPulse, 'AverageArterialPulseWaveDerivative', ToolBox)

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
    % print('-f102','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRApulseVsBackground.png'))) ;
    % print('-f103','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVpulseVsBackground.png'))) ;
    % print('-f108','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRAfilteredPulse.png'))) ;
    % print('-f110','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVfilteredPulse.png'))) ;
    % print('-f24','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'flattenedDopplerHeatMap.png'))) ;
    % print('-f77','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'zeroLagXcorr.png'))) ;
    % print('-f99','-dpng',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'timeLags.png'))) ;
    %
    % % eps
    % print('-f102','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRApulseVsBackground.eps'))) ;
    % print('-f103','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVpulseVsBackground.eps'))) ;
    % print('-f108','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRAfilteredPulse.eps'))) ;
    % print('-f110','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'CRVfilteredPulse.eps'))) ;
    % print('-f77','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'zeroLagXcorr.eps'))) ;
    % print('-f99','-depsc',fullfile(one_cycle_dir, 'pulseAnalysis', sprintf("%s_%s",ToolBox.main_foldername, 'timeLags.eps'))) ;
    %

end
