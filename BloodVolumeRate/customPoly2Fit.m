function [p1, p2, p3, rsquared, p1_err, p2_err, p3_err] = customPoly2Fit(x, y)
% Fits data to a polynomial of order 2 without using polyfit or fit
% Improves conditioning by centering and scaling the input data
% Input:
%   x - vector of n data points
%   y - vector of n data points
% Output:
%   p1, p2, p3 - coefficients of the polynomial y = p1*x^2 + p2*x + p3
%   rsquared - coefficient of determination (R^2) value
%   p1_err, p2_err, p3_err - standard errors of the coefficients

% Ensure x and y are column vectors
x = x(:);
y = y(:);

% Center and scale x
x_mean = mean(x);
x_std = std(x);
x_scaled = (x - x_mean) / x_std;

% Set up the design matrix for a 2nd order polynomial fit
X = [x_scaled.^2, x_scaled, ones(size(x_scaled))];

% Solve for the coefficients using the normal equation
p = (X' * X) \ (X' * y);

% Transform coefficients back to the original scale
p1 = p(1) / (x_std^2); % Coefficient for x^2
p2 = p(2) / x_std - 2 * p1 * x_mean; % Coefficient for x
p3 = p(3) - p2 * x_mean + p1 * x_mean^2; % Constant term

% Calculate fitted y values
y_fit = p1 * x.^2 + p2 * x + p3;

% Calculate R-squared value
ss_res = sum((y - y_fit).^2); % sum of squares of residuals
ss_tot = sum((y - mean(y)).^2); % total sum of squares
rsquared = 1 - (ss_res / ss_tot); % coefficient of determination

n = length(y); % number of data points
p_len = length(p); % number of parameters (3 for quadratic fit)

if n == p_len % The polynomial fit is perfect
    p1_err = 0;
    p2_err = 0;
    p3_err = 0;
else
    % Estimate the variance of residuals
    variance = ss_res / (n - p_len); % estimated variance of residuals

    % Calculate covariance matrix
    cov_matrix = variance .* inv(X' * X);

    % Extract standard errors (square root of diagonal elements of covariance matrix)
    p1_err = sqrt(cov_matrix(1,1)) / (x_std^2);
    p2_err = sqrt(cov_matrix(2,2)) / x_std;
    p3_err = sqrt(cov_matrix(3,3));
end
