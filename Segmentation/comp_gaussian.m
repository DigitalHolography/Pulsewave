function [G] = comp_gaussian(x, sigma)
% Function for exercise sheet 1 - computes the value of a Gaussian
% function.

%DESCRIPTION: G = comp_gaussian(x, sigma)
%       Computes the value of a zero-mean Gaussian function
%       (i.e. normal  distribution)

% INPUTS:
%       x (double vector) - Function input values (x-axis)
%       sigma (double scalar)  - Standard deviation (determines the 'width' of the Gaussian)

%OUTPUTS:
%       G (double vector) - Evaluated Gaussian function

%Function dependencies: None

%Author: Edward James, Medical Imaging CDT MRes, 2017

%Compute the Gaussian function
G = (1 / (sigma * (sqrt(2 * pi))) * exp((-1 * (x .^ 2)) / (2 * sigma ^ 2)));
end
