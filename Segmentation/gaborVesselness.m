function [vessel_mask, combined_response] = gaborVesselness(input_image, name, ToolBox)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Step 1: Preprocessing
% Convert to grayscale if necessary
if size(input_image, 3) == 3
    input_image = rgb2gray(input_image);
end

% Normalize the image to the range [0, 1]
input_image = im2double(input_image);

% Apply Gaussian smoothing to reduce noise
smoothed_image = imgaussfilt(input_image, 2); % Adjust sigma as needed

% Step 2: Gabor Filtering
% Define Gabor filter parameters
params = ToolBox.getParams;
range = params.json.Mask.VesselnessGaborRange;
step = params.json.Mask.VesselnessGaborStep;
wavelength = (range(1):step:range(2)); % Adjust based on vessel thickness
orientation = 0:30:150; % Multiple orientations to capture all vessel directions
gabor_bank = gabor(wavelength, orientation);

% Apply Gabor filters
gabor_magnitude = imgaborfilt(smoothed_image, gabor_bank);

% Combine responses from all orientations
combined_response = mean(gabor_magnitude, 3);

% Step 3: Postprocessing
% Normalize the combined response
combined_response = mat2gray(combined_response);

% Thresholding to create a binary mask
vessel_mask = imbinarize(combined_response, 'adaptive');

% Step 4: Morphological Operations
% Clean up the mask using morphological operations
vessel_mask = bwareaopen(vessel_mask, 50); % Remove small objects
vessel_mask = imclose(vessel_mask, strel('disk', 2)); % Close small gaps

if nargin == 3
    saveImage(combined_response, ToolBox, sprintf('%s_gabor_img.png', name), isStep = true)
    saveImage(vessel_mask, ToolBox, sprintf('%s_gabor_mask.png', name), isStep = true)
end

end
