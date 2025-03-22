function [oneCycleVideo, selectedPulseIdx, cyclesSignal, oneCycleVideoM0] = createOneCycle(video, videoM0, mask, sysIdxList, numInterp)
% createOneCycle - Identifies pulse cycles and averages them into one video.
% Inputs:
%   video: Input video data (3D or 4D array).
%   videoM0: M0 video data (3D or 4D array).
%   mask: Mask to isolate the region of interest.
%   sysIdxList: List of systolic indices in the video.
%   numInterp: Number of interpolation points for each cycle.
% Outputs:
%   oneCycleVideo: Averaged video of one cardiac cycle.
%   selectedPulseIdx: Indices of selected pulses.
%   cyclesSignal: Signal for each cycle.
%   oneCycleVideoM0: Averaged M0 video of one cardiac cycle.

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

tic;

% Initialize variables
numSys = length(sysIdxList) - 1; % Number of detected systoles
[numX, numY, ~] = size(video);
cyclesSignal = zeros(numSys, numInterp); % Preallocate for cycle signals

if numSys > 0 % At least two systoles detected
    % Calculate the signal for each cycle
    for ii = 1:numSys
        interp_range = linspace(sysIdxList(ii), sysIdxList(ii + 1) - 1, numInterp);
        frame_range = sysIdxList(ii):sysIdxList(ii + 1) - 1;
        tmp = squeeze(sum(video(:, :, frame_range) .* mask, [1 2])) / nnz(mask);
        cyclesSignal(ii, :) = interp1(frame_range, tmp, interp_range);
    end

    % Check pulse integrity using cross-correlation
    cleanSignal = squeeze(mean(cyclesSignal, 1)); % Average signal across cycles
    cleanSignal = cleanSignal - mean(cleanSignal); % Zero-mean normalization
    dataReliabilityIndex = zeros(numSys, 1); % Preallocate for reliability index

    for ii = 1:numSys
        signal = cyclesSignal(ii, :);
        signal = signal - mean(signal); % Zero-mean normalization
        C = sum(signal .* cleanSignal) / numInterp; % Cross-correlation
        dataReliabilityIndex(ii) = C; % Store reliability index
    end

    % Normalize reliability index to percentage
    dataReliabilityIndex = dataReliabilityIndex ./ max(dataReliabilityIndex) * 100;

    % Display reliability index for each pulse
    for ii = 1:numSys
        fprintf("        Data reliability index for pulse %d: %0.2f %%\n", ii, dataReliabilityIndex(ii));
    end

    % Select pulses based on reliability threshold
    cycles_accepted = (dataReliabilityIndex > params.oneCycle_outNoiseThreshold);
    selectedPulseIdx = find(cycles_accepted);

    if ~isempty(selectedPulseIdx) % At least one cycle accepted
        % Preallocate for single cycles
        single_cycles = zeros(numX, numY, numInterp, numSys);
        single_cyclesM0 = zeros(numX, numY, numInterp, numSys);

        % Interpolate and average selected cycles
        for jj = 1:length(selectedPulseIdx)
            ii = selectedPulseIdx(jj);
            interp_range = linspace(sysIdxList(ii), sysIdxList(ii + 1) - 1, numInterp);
            frame_range = sysIdxList(ii):sysIdxList(ii + 1) - 1;

            parfor id_x = 1:numX

                for id_y = 1:numY
                    single_cycles(id_x, id_y, :, ii) = interp1(frame_range, squeeze(video(id_x, id_y, frame_range)), interp_range);
                    single_cyclesM0(id_x, id_y, :, ii) = interp1(frame_range, squeeze(videoM0(id_x, id_y, frame_range)), interp_range);
                end

            end

        end

        % Average selected cycles
        oneCycleVideo = squeeze(sum(single_cycles, 4)) / length(selectedPulseIdx);
        oneCycleVideoM0 = squeeze(sum(single_cyclesM0, 4)) / length(selectedPulseIdx);
    else
        warning("No cycles met the reliability threshold. Using all cycles.");
        selectedPulseIdx = 1:numSys; % Use all cycles if none meet the threshold
        oneCycleVideo = squeeze(mean(single_cycles, 4));
        oneCycleVideoM0 = squeeze(mean(single_cyclesM0, 4));
    end

else % Less than two systoles detected
    warning("Less than two systoles detected. Using entire video as one cycle.");
    cyclesSignal = squeeze(sum(video(:, :, :) .* mask, [1 2])) / nnz(mask);
    cyclesSignal = interp1(1:size(video, 3), cyclesSignal, linspace(1, size(video, 3), numInterp));
    oneCycleVideo = interp3(video, size(video, 1), size(video, 2), numInterp);
    oneCycleVideoM0 = interp3(videoM0, size(video, 1), size(video, 2), numInterp);
    selectedPulseIdx = 0; % Indicate no valid pulses
end

% Shift the cycle to start and end with bottom diastole
oneP = squeeze(sum(oneCycleVideo .* mask, [1 2])) / nnz(mask);
[~, shift] = min(oneP); % Find bottom systole
oneCycleVideo = circshift(oneCycleVideo, -shift, 3);
oneCycleVideoM0 = circshift(oneCycleVideoM0, -shift, 3);
cyclesSignal = circshift(cyclesSignal, -shift, 2); % Shift pulse to start & end with bottom diastole

fprintf('    - CreateOneCycle took %ds\n', round(toc));
end
