function If = vesselness_filter(I, sigma, beta)
% VESSELNESS_FILTER filters a 2D image resulting in the enhancement of
% line-like features that correspond to blood vessels.
%
% Performs a 'vesselness filtering' operation.
%
% INPUTS:
%       I (double matrix) - image to be filtered.
%
%       sigma (double scalar, > 0)  - Standard deviations in each direction
%       of the Gaussian Kernel (determines the width or scale of the
%       Gaussian). This in turn controls the scale of the vessels
%       highlighted in the output image. A higher sigma value will
%       highlight vessels of a wider diameter, and vice versa.
%
%       beta (double scalar) - Controls the sensitivity of the line filter
%       to the blobness measure and is set to 0.5 in the paper by Frangi et
%       al.  A higher beta value will give more weighting in the output
%       image to pixels containing blob-like features, and vice versa.
%
% OUTPUTS:
%       If (double matrix) - Filtered image.
%
% FUNCTION DEPENDENCIES:
%	    1) gaussian_imfilter_sep - filters an image using an nth order
%       Gaussian derivative.
%
% AUTHOR:
%       Edward James, Medical Imaging CDT MRes , UCL, 2017

%Compute dimensions of input image
dims = size(I);
cols = dims(2);
rows = dims(1);

clear dims % clear redundant variables from the memory to speed up function

gamma = 1; % as no particular scaling parameter is preferred in this case

% Calculate the 4 second order Gaussian derivatives
A = sigma ^ gamma;

% For Ixx
Ixx = gaussian_imfilter_sep(I, 2, 0, sigma);
Ixx = A * Ixx;

% For Ixy
Ixy = gaussian_imfilter_sep(I, 1, 1, sigma);
Ixy = A * Ixy;

% For Iyx
Iyx = Ixy;

% For Iyy
Iyy = gaussian_imfilter_sep(I, 0, 2, sigma);
Iyy = A * Iyy;

clear I sigma gamma A

% Insert these values into a 4D matrix containing (rows x cols) by (2x2)
% Hessian matrices
H(1, 1, 1:rows, 1:cols) = Ixx(1:rows, 1:cols);
H(1, 2, 1:rows, 1:cols) = Ixy(1:rows, 1:cols);
H(2, 1, 1:rows, 1:cols) = Iyx(1:rows, 1:cols);
H(2, 2, 1:rows, 1:cols) = Iyy(1:rows, 1:cols);

clear Ixx Ixy Iyx  Iyy

% Find the eigenvalues of each (2x2) Hessian matrix

% Pre-allocate matrices of zeros to lambda_1 and lambda_2 for each pixel
lambda_1 = zeros(rows, cols);
lambda_2 = zeros(rows, cols);

for y = 1:rows

    for x = 1:cols
        % Compute the two eigen values at each pixel
        eig_val = eig(H(:, :, y, x));
        % Assign these eigen values to lambda_1 or lambda_2 so that the
        % magnitude of lambda_1 is less than or equal to the magnitude of
        % lambda_2
        if abs(eig_val) == sort(abs(eig_val))
            lambda_1(y, x) = eig_val(1);
            lambda_2(y, x) = eig_val(2);
        else
            lambda_1(y, x) = eig_val(2);
            lambda_2(y, x) = eig_val(1);
        end

    end

end

clear H eig_val x y

% Assign an R_blob value for each pixel (i.e. a measure of blobness in 2D)
R_blob = abs(lambda_1) ./ abs(lambda_2);

% Assign a temporary c value for each pixel
c_temp = sqrt((lambda_1 .^ 2) + (lambda_2 .^ 2));
% Assign a final value of c (which is the maximum value of c_temp computed
% over the whole image)
c = max(max(c_temp));

% Calculate the value of each pixel of If according c_temp, R_blob, beta
% and c values

A = R_blob .^ 2;
B = 2 * beta ^ 2;
C = exp(- (A ./ B));
D = 2 * c ^ 2;
E = c_temp ./ D; %using original expression in exercise sheet #5
F = 1 - exp(-E);
If = C .* F;

clear A B C D E F R_blob beta c_temp lambda_1

% Convert the pixels of If to zero where lambda_2 is greater than zero

% Create a matrix to contain ones and zeros depending on values in lambda_2
lambda_2_logical = zeros(rows, cols);
% Convert all corresponding values of lambda_2_logical to 1, where lambda_2
% is equal to or less than zero
lambda_2_logical(lambda_2 <= -0.005) = 1;
%lambda_2_logical(lambda_2<=0)=1;
clear lambda_2 rows cols
%Multiply the resulting matrix of ones and zeros element wise by If
If = If .* lambda_2_logical;
end
