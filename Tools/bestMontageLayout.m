function [rows, cols] = bestMontageLayout(n)
% bestMontageLayout: Computes the optimal layout for displaying 'n' figures
% on a 16:9 screen while keeping the aspect ratio balanced.
%
% Inputs:
%   n - Number of figures
% Outputs:
%   rows - Optimal number of rows
%   cols - Optimal number of columns

% Start with a square root approach
bestRatio = inf; % Initialize with a high value
bestRows = 1;
bestCols = n;

for rows = 1:ceil(sqrt(n))
    cols = ceil(n / rows);
    aspectRatio = (cols / rows) / (16/9); % Normalize to 16:9

    % Closer to 1 means better fit to 16:9
    if abs(aspectRatio - 1) < bestRatio
        bestRatio = abs(aspectRatio - 1);
        bestRows = rows;
        bestCols = cols;
    end

end

rows = bestRows;
cols = bestCols;
end
