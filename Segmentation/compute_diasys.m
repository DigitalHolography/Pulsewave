function [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_ff_video, firstMaskArteryClean)

    [~, ~, numFrames] = size(M0_ff_video);
    [sys_index_list, ~] = find_systole_index(M0_ff_video, firstMaskArteryClean);

    numSys = numel(sys_index_list); % number of systoles
    fpCycle = round(numFrames / numSys); % Frames per cycle

    sysindexes = [];
    diaindexes = [];

    for idx = 1:numSys
        try
            % Calculate sysindexes and ensure the values stay within the valid range
            start_idx = sys_index_list(idx);
            end_idx = sys_index_list(idx) + round(fpCycle * 0.1);
            sysindexes = [sysindexes, start_idx:min(end_idx, size(M0_ff_video, 3))];
        catch
        end

        try
            % Calculate diaindexes and ensure the values stay within the valid range
            start_idx = sys_index_list(idx) - round(fpCycle * 0.15);
            end_idx = sys_index_list(idx) - round(fpCycle * 0.05);
            diaindexes = [diaindexes, max(start_idx, 1):min(end_idx, size(M0_ff_video, 3))];
        catch
        end
    end

    % Ensure sysindexes and diaindexes are within the bounds of the video size
    sysindexes = sysindexes(sysindexes >= 1 & sysindexes <= size(M0_ff_video, 3));
    diaindexes = diaindexes(diaindexes >= 1 & diaindexes <= size(M0_ff_video, 3));

    % Compute the mean images
    M0_Systole_img = mean(M0_ff_video(:, :, sysindexes), 3);
    M0_Diastole_img = mean(M0_ff_video(:, :, diaindexes), 3);
end