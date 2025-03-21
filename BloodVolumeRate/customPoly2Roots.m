function [r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err)
% Calculates roots of a quadratic polynomial and propagates errors
% Input:
%   p1, p2, p3 - coefficients of the polynomial y = p1*x^2 + p2*x + p3
%   p1_err, p2_err, p3_err - standard errors of the coefficients
% Output:
%   r1, r2 - roots of the polynomial
%   r1_err, r2_err - propagated errors of the roots

% Calculate the discriminant and its square root
discriminant = p2 ^ 2 - 4 * p1 * p3;

if discriminant < 0
    %         disp('The polynomial has complex roots.');
    r1 = 0;
    r2 = 0;
    r1_err = 0;
    r2_err = 0;
    return
end

sqrt_discriminant = sqrt(discriminant);

% Compute the roots
r1 = (-p2 + sqrt_discriminant) / (2 * p1);
r2 = (-p2 - sqrt_discriminant) / (2 * p1);

% Partial derivatives for error propagation
dr1_dp1 = (2 * p3 / sqrt_discriminant - p2 / p1) / (2 * p1 ^ 2);
dr1_dp2 = (1 / (2 * p1)) * (1 + p2 / sqrt_discriminant);
dr1_dp3 = -1 / (sqrt_discriminant * p1);

dr2_dp1 =- (2 * p3 / sqrt_discriminant + p2 / p1) / (2 * p1 ^ 2);
dr2_dp2 = (1 / (2 * p1)) * (1 - p2 / sqrt_discriminant);
dr2_dp3 = -1 / (sqrt_discriminant * p1);

% Propagate errors for r1 and r2 using the partial derivatives
r1_err = sqrt((dr1_dp1 * p1_err) ^ 2 + (dr1_dp2 * p2_err) ^ 2 + (dr1_dp3 * p3_err) ^ 2);
r2_err = sqrt((dr2_dp1 * p1_err) ^ 2 + (dr2_dp2 * p2_err) ^ 2 + (dr2_dp3 * p3_err) ^ 2);
end
