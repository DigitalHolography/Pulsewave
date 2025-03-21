function [mask, R] = diskMask(numX, numY, R1, R2, opt)

% numX, numY : size of the mask
% If nargin = 3:
% D1: upper bound of the diameter
% If nargin = 4:
% D1: lower bound of the diameter
% D2: upper bound of the diameter

arguments
    numX
    numY
    R1
    R2 = 2
    opt.center = [1/2, 1/2]
end

if R1 > R2
    error("R1 must be lower than R2")
end

x_c = opt.center(1);
y_c = opt.center(2);

x = linspace(0, 1, numX);
y = linspace(0, 1, numY);

[X, Y] = meshgrid(y, x);
R = (((X - x_c)) .^ 2 + ((Y - y_c)) .^ 2);

if nargin == 3
    % Sup & Equal so that r2 = 0.5 includes the radius
    mask = R <= R1 ^ 2;
elseif nargin == 4
    % Inf & Equal so that r1 = 0 excludes the whole image
    % Sup & Equal so that r2 = 0.5 includes the radius
    mask = (R > R1 ^ 2) & (R <= R2 ^ 2);
end

end
