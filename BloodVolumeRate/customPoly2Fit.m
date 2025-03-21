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

% Center and scale x to improve conditioning
x_mean = mean(x);
x_std = std(x);
x_scaled = (x - x_mean) / x_std;

% Construct the design matrix
X = [x_scaled .^ 2, x_scaled, ones(length(x_scaled), 1)];

% Solve for the coefficients using normal equations
p = (X' * X) \ (X' * y);

% Transform coefficients back to original scale
p1 = p(1) / (x_std ^ 2); % Rescale p1
p2 = p(2) / x_std; % Rescale p2
p3 = p(3) - (p1 * x_mean ^ 2 + p2 * x_mean); % Adjust p3

% Compute fitted values
y_fit = p1 * x .^ 2 + p2 * x + p3;

% Compute R-squared value
ss_res = sum((y - y_fit) .^ 2); % Sum of squares of residuals
ss_tot = sum((y - mean(y)) .^ 2); % Total sum of squares
rsquared = 1 - (ss_res / ss_tot); % Coefficient of determination

% Compute standard errors of the coefficients
n = length(y);
p_len = length(p);

if n > p_len % Avoid division by zero in perfect fits
    variance = ss_res / (n - p_len); % Estimated variance of residuals
    cov_matrix = variance .* inv(X' * X);

    % Extract standard errors
    p1_err = sqrt(cov_matrix(1, 1)) / (x_std ^ 2);
    p2_err = sqrt(cov_matrix(2, 2)) / x_std;
    p3_err = sqrt(cov_matrix(3, 3));
else
    p1_err = 0;
    p2_err = 0;
    p3_err = 0;
end

end
