function [p1, p2, p3, rsquared, p1_err, p2_err, p3_err] = customPoly2Fit(x, y)
    % Fits data to a polynomial of order 2 without using polyfit or fit
    % Calculates coefficients and their standard errors
    % Input:
    %   x - vector of x data points
    %   y - vector of y data points
    % Output:
    %   p1, p2, p3 - coefficients of the polynomial y = p1*x^2 + p2*x + p3
    %   rsquared - coefficient of determination (R^2) value
    %   p1_err, p2_err, p3_err - standard errors of the coefficients

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Set up the design matrix for a 2nd order polynomial fit
    X = [x.^2, x, ones(size(x))];
    
    % Solve for the coefficients using the normal equation
    p = (X' * X) \ (X' * y);
    
    % Extract coefficients
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    
    % Calculate fitted y values
    y_fit = X * p;
    
    % Calculate R-squared value
    ss_res = sum((y - y_fit).^2); % sum of squares of residuals
    ss_tot = sum((y - mean(y)).^2); % total sum of squares
    rsquared = 1 - (ss_res / ss_tot); % coefficient of determination
    
    % Estimate the variance of residuals
    n = length(y); % number of data points
    p_len = length(p); % number of parameters (3 for quadratic fit)
    variance = ss_res / (n - p_len); % estimated variance of residuals
    
    % Calculate covariance matrix
    cov_matrix = variance * inv(X' * X);
    
    % Extract standard errors (square root of diagonal elements of covariance matrix)
    p1_err = sqrt(cov_matrix(1,1));
    p2_err = sqrt(cov_matrix(2,2));
    p3_err = sqrt(cov_matrix(3,3));
end
