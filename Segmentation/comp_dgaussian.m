function [Gn] = comp_dgaussian(x, sigma, n)
% Function for exercise sheet 3 - computes the nth order derivative of a
% Gaussian function

% INPUTS:
%       x (double vector) - Function input values (x-axis)
%       sigma (double scalar)  - Standard deviation
%       n (integer scalar) - determines the order of the derivative
%       (n >= 0)

%OUTPUTS:
%       Gn (double vector) - Evaluated nth order derivative of a Gaussian function

%Function dependencies:
%       comp_gaussian(x,sigma) - computes the value of a Gaussian function
%       - G(x,sigma)
%       comp_hermite_rec(n,x) - omputes the physicist's Hermite polynomial
%       of order n, using a recursive method - Hn (n,x)
%       (therefore note when n = 0, Gn = G(x,sigma)

%Author: Edward James, Medical Imaging CDT MRes, 2017

%First form expressions for the three terms in the expression for an nth order
%derivative of a Gaussian Function:

A = (-1 / (sigma * sqrt(2))) ^ n;

B = comp_hermite_rec(n, x / (sigma * sqrt(2)));

C = comp_gaussian(x, sigma);

%Multiply all three terms together to obtain Gn:
Gn = A * B .* C;
end
