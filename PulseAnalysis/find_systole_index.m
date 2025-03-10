function [sys_index_list, fullPulse, sys_max_list, sys_min_list] = find_systole_index(video, maskArtery)
% This function extracts systole and diastole indices from a noisy video signal.
% Inputs:
%   - video: 3D video data (height x width x time)
%   - maskArtery: Binary mask indicating the artery region
% Outputs:
%   - sys_index_list: Indices of systole peaks
%   - fullPulse: The extracted pulse signal
%   - sys_max_list: Indices of local maxima within each cycle
%   - sys_min_list: Indices of local minima within each cycle

% Get global toolbox and parameters
TB = getGlobalToolBox;
params = TB.getParams;

% Step 1: Extract the pulse signal from the video using the artery mask
fullPulse = squeeze(sum(video .* maskArtery, [1 2]) / nnz(maskArtery));

% Step 2: Preprocess the pulse signal
pulse_init = detrend(fullPulse); % Remove linear trend
pulse_init = pulse_init - mean(pulse_init, "all"); % Remove mean to center the signal

% Step 1: Detect and replace outliers
mad_threshold = 3;
median_pulse = median(pulse_init);
mad_pulse = mad(pulse_init, 1);
outlier_indices = abs(pulse_init - median_pulse) > mad_threshold * mad_pulse;
pulse_init(outlier_indices) = interp1(find(~outlier_indices), pulse_init(~outlier_indices), find(outlier_indices), 'linear');

% Step 2: Smooth the signal
window_size = 15;
polynomial_order = 3;
pulse_init = sgolayfilt(pulse_init, polynomial_order, window_size);

% Step 3: Handle outliers using Hampel filter
% The Hampel filter replaces outliers with the median of a sliding window
window_size = 3; % Half-window size for Hampel filter
fullArterialPulseRmOut = hampel(pulse_init, window_size);

% Step 4: Smooth the signal to reduce noise
fullArterialPulseClean = smoothdata(fullArterialPulseRmOut, 'lowess');

% Step 5: Compute the derivative of the smoothed signal
% The derivative helps identify sharp changes corresponding to systole peaks
diff_signal = diff(fullArterialPulseClean);

% Step 6: Detect systole peaks using findpeaks
% Peaks in the derivative correspond to systole events
min_peak_height = max(diff_signal) * params.systoleThreshold; % Adaptive threshold
min_peak_distance = 10; % Minimum distance between peaks (to avoid duplicates)
[~, sys_index_list] = findpeaks(diff_signal, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

% Step 7: Validate and clean up the detected peaks
% Remove cycles that are too short (less than 10 points apart)
i = 1;
while i < numel(sys_index_list)
    if sys_index_list(i+1) - sys_index_list(i) < 10
        sys_index_list(i+1) = []; % Remove the next peak if too close
    else
        i = i + 1; % Move to the next peak
    end
end

% Step 8: Find local maxima and minima within each cycle
sys_max_list = []; % List of local maxima (systole peaks)
sys_min_list = []; % List of local minima (diastole troughs)

% Use parallel processing for efficiency
parfor i = 1:(numel(sys_index_list) - 1)
    % Find the maximum within the current cycle
    [~, amax] = max(pulse_init(sys_index_list(i):sys_index_list(i+1)));
    sys_max_list(i) = sys_index_list(i) + amax - 1;

    % Find the minimum within the current cycle
    [~, amin] = min(pulse_init(sys_index_list(i):sys_index_list(i+1)));
    sys_min_list(i) = sys_index_list(i) + amin - 1;
end

% Step 10: Error handling
if isempty(sys_index_list)
    error('No systole peaks detected. Check signal quality or adjust parameters.');
end

%     figure(3)
%     plot(diff_signal);
%     title('pulse init derivate');
%     1;

%     spectrum_signal = fft(diff_signal);
%     w = 1:length(spectrum_signal);
%     figure(41);
%     plot(w, abs(spectrum_signal));
%     signal_center = floor(length(spectrum_signal)/2);
%     [~, omega1] = max(abs(spectrum_signal(1:signal_center)));
%     [~, omega2] = max(abs(spectrum_signal(signal_center:end)));
%     omega2 = omega2 + signal_center - 1;
%     tmp = zeros(1, length(spectrum_signal));
%     tmp(omega1) = spectrum_signal(omega1);
%     tmp(omega2) = spectrum_signal(omega2);
%     tmp = ifft(tmp);
%     figure(42)
%     plot(w, real(tmp), w, diff_signal);
%     pseudo_period = length(w) / omega1;
%     distance_between_peaks = 0.7 * pseudo_period;
%
%     [~,index_list] = findpeaks(diff_signal,1:length(diff_signal),'MinPeakDistance',distance_between_peaks,'SortStr','descend');
%     [~,index_list_tmp] = findpeaks(tmp,1:length(tmp),'MinPeakDistance',distance_between_peaks,'SortStr','descend');
%     index_list = sort(index_list(1:length(index_list)));
%     fprintf("num Cycles: %d\n", length(index_list)-1);
%     sys_index_list = zeros(1, length(index_list_tmp));
%     for ii = 1:length(sys_index_list)
%         tab = abs(index_list-index_list_tmp(ii));
%         [~, new_index] = min(tab);
%         sys_index_list(ii) = index_list(new_index);
%     end
%     sys_index_list = sort(sys_index_list, 'ascend');
end
