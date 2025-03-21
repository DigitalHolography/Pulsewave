function [fitParams, fittedCurve] = fitGaussian(binCenters, radialAverage)
% FITGAUSSIAN Fits a Gaussian function to radial average data.
%
% Inputs:
%   binCenters: Centers of the radial bins (1D array).
%   radialAverage: Radial average intensities corresponding to binCenters (1D array).
%
% Outputs:
%   fitParams: Fitted Gaussian parameters [A, mu, sigma, C].
%   fittedCurve: Fitted Gaussian curve evaluated at binCenters.

% Define the Gaussian function
gaussianFunc = @(p, x) p(1) * exp(- ((x - p(2)) / p(3)) .^ 2) + p(4);

% Initial guess for the Gaussian parameters [A, mu, sigma, C]
[~, maxIdx] = max(radialAverage);
initialGuess = [max(radialAverage), binCenters(maxIdx), (max(binCenters) - min(binCenters)) * 0.1, min(radialAverage)];

% Define bounds for [A, mu, sigma, C]
lb = [0, min(binCenters), 0, -Inf]; % Lower bounds
ub = [Inf, max(binCenters), Inf, Inf]; % Upper bounds

% Set optimization options
options = optimoptions('lsqcurvefit', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000, 'Display', 'off');

% Suppress warnings
warning('off', 'optim:lsqcurvefit:JacobianIllConditioned');
warning('off', 'optim:lsqcurvefit:IterationLimitExceeded');

% Perform the fit
fitParams = lsqcurvefit(gaussianFunc, initialGuess, binCenters, radialAverage, lb, ub, options);

% Turn warnings back on (optional)
warning('on', 'optim:lsqcurvefit:JacobianIllConditioned');
warning('on', 'optim:lsqcurvefit:IterationLimitExceeded');

% Compute the fitted Gaussian curve
fittedCurve = gaussianFunc(fitParams, binCenters);
end
