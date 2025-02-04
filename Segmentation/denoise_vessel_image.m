function [M0_nlm_img, M0_binary_img] = denoise_vessel_image(M0_ff_img, diaphragm, name, ToolBox)
    % Input: M0_ff_img - Noisy vessel image (numX x numY)
    % Output: M0_denoised_img - Denoised vessel image

    PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);

    % Step 1: Preprocessing
    % Convert to grayscale if the image is RGB
    if size(M0_ff_img, 3) == 3
        M0_gray_img = rgb2gray(M0_ff_img);
    else
        M0_gray_img = M0_ff_img;
    end

    % Normalize the image to [0, 1]
    M0_norm_img = rescale(M0_gray_img);

    % Step 2: Denoising
    % Apply Gaussian smoothing
    sigma = PW_params.gauss_filt_size_for_barycenter; % Standard deviation for Gaussian filter
    M0_gaussian_img = imgaussfilt(M0_norm_img, sigma);

    % Apply Non-Local Means (NLM) denoising
    patch_size = 5; % Size of patches for NLM
    search_window = 11; % Size of search window for NLM
    M0_nlm_img = imnlmfilt(M0_gaussian_img, 'DegreeOfSmoothing', 0.05, ...
                           'SearchWindowSize', search_window, ...
                           'ComparisonWindowSize', patch_size);

    % Step 3: Vessel Enhancement
    % Apply Frangi Vesselness Filter
    M0_vesselness_img = vesselness_filter(M0_nlm_img, PW_params.masks_vesselness_sigma, PW_params.masks_vesselness_beta);

    % Thresholding to segment vessels (optional)

    M0_vesselness_img = M0_vesselness_img .* diaphragm;
    threshold = graythresh(M0_vesselness_img); % Otsu's method
    M0_binary_img = imbinarize(M0_vesselness_img, threshold);

    saveImage(rescale(M0_nlm_img), ToolBox, sprintf('%s_M0_desnoised_img.png', name), isStep = true)
    saveImage(M0_binary_img, ToolBox, sprintf('%s_maskVessel.png', name), isStep = true)
end