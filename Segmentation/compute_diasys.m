function [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_ff_video, firstMaskArteryClean)

    TB = getGlobalToolBox;
    [~, ~, numFrames] = size(M0_ff_video);
    fullTime = linspace(0, numFrames * TB.stride / TB.fs / 1000, numFrames);

    [~, fullPulse, sys_max_list, sys_min_list] = find_systole_index(M0_ff_video, firstMaskArteryClean);

    figure('Visible', 'off')
    hold on
    plot(fullTime, fullPulse, 'k--', 'LineWidth', 2)

    numSysMax = numel(sys_max_list); % number of systoles
    numSysMin = numel(sys_min_list); % number of systoles
    fpCycleMax = round(numFrames / numSysMax); % Frames per cycle
    fpCycleMin = round(numFrames / numSysMin); % Frames per cycle

    sysindexes = [];
    diasindexes = [];

    for idx = 1:numSysMax
        try
            % Calculate sysindexes and ensure the values stay within the valid range
            start_idx = sys_max_list(idx)- round(fpCycleMax * 0.05);
            end_idx = sys_max_list(idx) + round(fpCycleMax * 0.1);
            sys_range = start_idx:min(end_idx, numFrames);
            sysindexes = [sysindexes, sys_range];
            plot(fullTime(sys_range), fullPulse(sys_range), 'r-', 'LineWidth', 2)
        catch
        end
    end

    for idx = 1:numSysMin
        try
            % Calculate diaindexes and ensure the values stay within the valid range
            start_idx = sys_min_list(idx) - round(fpCycleMin * 0.15);
            end_idx = sys_min_list(idx);
            dias_range = max(start_idx, 1):min(end_idx, numFrames);
            diasindexes = [diasindexes, dias_range];
            plot(fullTime(dias_range), fullPulse(dias_range), 'b-', 'LineWidth', 2)
        catch
        end
    end

    % Ensure sysindexes and diaindexes are within the bounds of the video size
    sysindexes = sysindexes(sysindexes >= 1 & sysindexes <= numFrames);
    diasindexes = diasindexes(diasindexes >= 1 & diasindexes <= numFrames);

    % Compute the mean images
    M0_Systole_img = mean(M0_ff_video(:, :, sysindexes), 3);
    M0_Diastole_img = mean(M0_ff_video(:, :, diasindexes), 3);

    % Adjust axes
    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), axP(3), axP(4) * 1.07])
    box on
    set(gca, 'LineWidth', 2)
    set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
    title('diastole and systole')

    exportgraphics(gca, fullfile(TB.path_png, 'mask', 'steps', ...
        sprintf('%s_vessel_20_plot_diasys.png', TB.main_foldername)))

end