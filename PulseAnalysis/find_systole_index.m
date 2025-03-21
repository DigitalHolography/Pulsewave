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

% Step 1: Extract the pulse signal from the video using the artery mask
fullPulse = squeeze(sum(video .* maskArtery, [1 2]) / nnz(maskArtery));

% Step 2: Preprocess the pulse signal
pulse_init = double(detrend(fullPulse)); % Remove linear trend
pulse_init = pulse_init - mean(pulse_init); % Remove mean to center the signal

% Step 3: Detect and replace outliers using a manual MAD calculation
mad_threshold = 3;
median_pulse = median(pulse_init); % Compute median
mad_pulse = median(abs(pulse_init - median_pulse)); % Manual MAD calculation

outlier_indices = abs(pulse_init - median_pulse) > mad_threshold * mad_pulse;
pulse_init(outlier_indices) = interp1(find(~outlier_indices), pulse_init(~outlier_indices), find(outlier_indices), 'linear', 'extrap');

% Step 4: Smooth the signal using Savitzky-Golay filter (Requires Signal Processing Toolbox)
window_size = 15;
polynomial_order = 3;
pulse_init = sgolayfilt(pulse_init, polynomial_order, window_size);

% Step 5: Replace Hampel filter with a manual implementation
window_size = 3;
fullArterialPulseRmOut = hampel_manual(pulse_init, window_size);

% Step 6: Smooth the signal to reduce noise (Replace Lowess smoothing with Moving Average)
fullArterialPulseClean = movmean(fullArterialPulseRmOut, window_size);

% Step 7: Compute the derivative of the smoothed signal
diff_signal = diff(fullArterialPulseClean);

% Step 8: Detect systole peaks using findpeaks
min_peak_height = max(diff_signal) * 0.3; % Adaptive threshold
min_peak_distance = 10; % Minimum distance between peaks
[~, sys_index_list] = findpeaks(diff_signal, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

% Step 9: Validate and clean up the detected peaks
sys_index_list = validate_peaks(sys_index_list, 10);

% Step 10: Find local maxima and minima within each cycle
sys_max_list = zeros(1, numel(sys_index_list) - 1);
sys_min_list = zeros(1, numel(sys_index_list) - 1);

for i = 1:(numel(sys_index_list) - 1)
    % Find the maximum within the current cycle
    [~, amax] = max(pulse_init(sys_index_list(i):sys_index_list(i + 1)));
    sys_max_list(i) = sys_index_list(i) + amax - 1;

    % Find the minimum within the current cycle
    [~, amin] = min(pulse_init(sys_index_list(i):sys_index_list(i + 1)));
    sys_min_list(i) = sys_index_list(i) + amin - 1;
end

% Step 11: Error handling
if isempty(sys_index_list)
    error('No systole peaks detected. Check signal quality or adjust parameters.');
end

end

%% **Manual Hampel Filter (Replacing hampel())**
function y = hampel_manual(x, k)
n = length(x);
y = x; % Initialize output

for i = (k + 1):(n - k)
    window = x(i - k:i + k);
    med = median(window);
    sigma = median(abs(window - med));
    threshold = 3 * sigma;

    if abs(x(i) - med) > threshold
        y(i) = med;
    end

end

end

%% **Validate Peaks (Removes peaks that are too close)**
function sys_index_list = validate_peaks(sys_index_list, min_distance)
i = 1;

while i < numel(sys_index_list)

    if sys_index_list(i + 1) - sys_index_list(i) < min_distance
        sys_index_list(i + 1) = [];
    else
        i = i + 1;
    end

end

end
