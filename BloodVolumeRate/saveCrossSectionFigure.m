function saveCrossSectionFigure(subImg, crossSectionWidth, ToolBox, name_section)
% Save a figure showing the cross-section of the blood vessel.
%
% Inputs:
%   subImg              - 2D array, the sub-image of the blood vessel.
%   tmp_section         - 1D array, the cross-section profile.
%   crossSectionWidth   - Scalar, the width of the cross-section.
%   sectionIdx          - Scalar, index of the current section.
%   ToolBox             - Struct, contains parameters and paths.
%   insert              - String, additional identifier for the filename.
%   name_section        - String, name of the section.

% Create figure
f = figure('Visible', 'off');

% Compute Section Cut
section_cut = squeeze(sum(subImg, 1));
section_cut(section_cut < 0) = 0;
cross_section_profile = (section_cut ./ max(section_cut)) * size(section_cut, 2);

% Define x-axis values
xAx = linspace(0, size(subImg, 1), size(subImg, 1));

% Display the sub-image
imagesc(xAx, xAx, subImg);
colormap("gray");
axis image;
hold on;

% Set aspect ratio
set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);

% Plot the cross-section profile
p = plot(xAx, cross_section_profile);
p.LineWidth = 2;
p.Color = 'red';
p.LineStyle = ':';

% Plot the cross-section width line
x = [round(size(subImg, 1) / 2) - round(crossSectionWidth / 2), ...
         round(size(subImg, 1) / 2) + round(crossSectionWidth / 2)];
y = [round(size(subImg, 1) / 2), round(size(subImg, 1) / 2)];
line(x, y, 'Color', 'red', 'LineWidth', 3);

% Turn off axes
axis off;

% Capture the figure
frame = getframe(gca);

% Save the figure
outputPath = fullfile(ToolBox.path_png, 'volumeRate', 'crossSection', ...
    sprintf('%s_cross_section_%s.png', ToolBox.main_foldername, name_section));
imwrite(frame.cdata, outputPath);

% Close the figure
close(f);
end
