function generateHealthReport(ToolBox, vRMS, maskArtery, ~)
% Define file paths
dataFilePath = fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_EF_main_outputs.txt'));
pdfPath = fullfile(ToolBox.path_pdf, sprintf("%s_EyeFlowReport.pdf", ToolBox.main_foldername));
directoryName = ToolBox.main_foldername;
%     params = ToolBox.getParams;
%     veinsAnalysis = params.json.VeinsAnalysis;
[~, ~, numFrames] = size(vRMS);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

% Set A4 dimensions in centimeters
a4Width = 21.0; % A4 width in cm
a4Height = 29.7; % A4 height in cm

% Define margins (2 cm on all sides)
margin = 2; % in cm

% Create a new figure for the report
fig = figure('Units', 'centimeters', 'Position', [0 0 a4Width a4Height], ...
    'Color', 'white', 'Visible', 'on', ...
    'PaperSize', [a4Width a4Height], 'PaperPosition', [0 0 a4Width a4Height], ...
    'PaperPositionMode', 'manual'); % Set paper size and position

% Read and parse the data file
data = parseDataFile(dataFilePath);

% Add logo to the top-right corner
addLogo(fig, 'eyeflow_logo.png', a4Width, a4Height);

% Add title with directory name
reportTitle = sprintf('EyeFlow Report - %s', directoryName);
addTitle(fig, reportTitle, margin, a4Width, a4Height);

% Add subtitle with current date/time
currentDateTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
subtitleText = sprintf('Generated on: %s', currentDateTime);
addSubtitle(fig, subtitleText, margin, a4Width, a4Height);

% Add data fields
yPos = 0.85 - (margin / a4Height); % Starting vertical position for annotations (normalized units)
yPos = addField(fig, sprintf('Heart Beat: %.1f bpm', data.heartBeat), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Systoles: %s', mat2str(round(data.systoleIndices * ToolBox.stride / ToolBox.fs / 1000, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Number of Cycles: %d', data.numCycles), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Max Systole Indices: %s', mat2str(round(data.maxSystoleIndices * ToolBox.stride / ToolBox.fs / 1000, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Min Systole Indices: %s', mat2str(round(data.minSystoleIndices * ToolBox.stride / ToolBox.fs / 1000, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Mean Blood Volume Rate (Artery): %.1f ± %.1f µL/min', data.meanBloodVolumeRateArtery, data.stdBloodVolumeRateArtery / 2), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Mean Blood Volume Rate (Vein): %.1f ± %.1f µL/min', data.meanBloodVolumeRateVein, data.stdBloodVolumeRateVein / 2), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Max Systole Blood Volume Rate (Artery): %.1f µL/min', data.maxSystoleBloodVolumeRateArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Min Diastole Blood Volume Rate (Artery): %.1f µL/min', data.minDiastoleBloodVolumeRateArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Stroke Volume (Artery): %.1f nL', data.strokeVolumeArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Total Volume (Artery): %.1f nL', data.totalVolumeArtery), yPos, margin, a4Width);

% Add signal plot
arterySignal = sum(vRMS .* maskArtery, [1 2]) / nnz(maskArtery);
addSignalPlot(fig, squeeze(arterySignal), t, yPos, margin, a4Width, a4Height);

% Save the figure as a PDF using print
print(fig, pdfPath, '-dpdf', '-fillpage'); % Export to PDF with A4 size and margins

% Close the figure
close(fig);
end

% Helper function to parse the data file
function data = parseDataFile(dataFilePath)
% Read the data from the text file
fileContent = fileread(dataFilePath);

% Extract values using regular expressions
data.heartBeat = str2double(regexp(fileContent, 'Heart beat: ([\d.]+)', 'tokens', 'once'));

% Extract Systole Indices and remove trailing comma
systoleIndicesStr = regexp(fileContent, 'Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.systoleIndices = str2num(systoleIndicesStr{1}); %#ok<ST2NM>

data.numCycles = str2double(regexp(fileContent, 'Number of Cycles: ([\d]+)', 'tokens', 'once'));

% Extract Max Systole Indices and remove trailing comma
maxSystoleIndicesStr = regexp(fileContent, 'Max Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.maxSystoleIndices = str2num(maxSystoleIndicesStr{1}); %#ok<ST2NM>

% Extract Min Systole Indices and remove trailing comma
minSystoleIndicesStr = regexp(fileContent, 'Min Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.minSystoleIndices = str2num(minSystoleIndicesStr{1}); %#ok<ST2NM>

data.meanBloodVolumeRateArtery = str2double(regexp(fileContent, 'Flow Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.stdBloodVolumeRateArtery = str2double(regexp(fileContent, 'Flow Rate Standard Deviation Artery : ([\d.]+)', 'tokens', 'once'));
data.meanBloodVolumeRateVein = str2double(regexp(fileContent, 'Flow Rate Vein : ([\d.]+)', 'tokens', 'once'));
data.stdBloodVolumeRateVein = str2double(regexp(fileContent, 'Flow Rate Standard Deviation Vein : ([\d.]+)', 'tokens', 'once'));
data.maxSystoleBloodVolumeRateArtery = str2double(regexp(fileContent, 'MaxSystole Blood Volume Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.minDiastoleBloodVolumeRateArtery = str2double(regexp(fileContent, 'MinDiastole Blood Volume Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.strokeVolumeArtery = str2double(regexp(fileContent, 'Stroke Volume Artery : ([\d.]+)', 'tokens', 'once'));
data.totalVolumeArtery = str2double(regexp(fileContent, 'Total Volume Artery : ([\d.]+)', 'tokens', 'once'));
end

% Helper function to add a title
function addTitle(fig, titleText, margin, a4Width, a4Height)
% Calculate normalized position for the title
titleX = margin / a4Width; % Normalized x position (2 cm margin)
titleY = 1 - (margin / a4Height); % Normalized y position (2 cm margin)

annotation(fig, 'textbox', [titleX titleY - 0.1 1 - 2 * titleX 0.1], 'String', titleText, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Interpreter', 'none');
end

% Helper function to add a data field
function yPos = addField(fig, fieldText, yPos, margin, a4Width)
% Calculate normalized position for the field
fieldX = margin / a4Width; % Normalized x position (2 cm margin)
fieldY = yPos; % Normalized y position (2 cm margin)

annotation(fig, 'textbox', [fieldX fieldY 1 - 2 * fieldX 0.05], 'String', fieldText, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'EdgeColor', 'none');
yPos = yPos - 0.02; % Move down for the next field
end

% Helper function to add the signal plot
function addSignalPlot(fig, signal, t, yPos, margin, a4Width, a4Height)
% Calculate normalized position for the plot
plotX = margin / a4Width; % Normalized x position (2 cm margin)
plotY = yPos - 0.3 - (margin / a4Height); % Normalized y position (2 cm margin)
plotWidth = 1 - 2 * plotX; % Normalized width
plotHeight = 0.25; % Normalized height

axes('Parent', fig, 'Position', [plotX plotY plotWidth plotHeight]); % Position the plot
plot(t, signal, 'r', 'LineWidth', 2);
title('Arterial Signal');
xlabel('Time (s)');
ylabel('Velocity (mm/s)');

pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4) * 1.07])
box on
end

% Helper function to add the logo
function addLogo(fig, logoPath, a4Width, a4Height)
% Load the logo image
logoImg = imread(logoPath);

% Calculate normalized position for the logo (top-right corner)
logoWidth = 3; % Width of the logo in cm
logoHeight = size(logoImg, 1) / size(logoImg, 2) * logoWidth; % Maintain aspect ratio
logoX = 1 - (logoWidth / a4Width); % Normalized x position
logoY = 1 - (logoHeight / a4Height); % Normalized y position

% Create axes for the logo
axes('Parent', fig, 'Units', 'normalized', 'Position', [logoX logoY logoWidth / a4Width logoHeight / a4Height]);
h = imshow(logoImg);
axis off;

% Set the transparency (alpha) of the logo to 50%
alpha(h, 0.5);
end

function addSubtitle(fig, subtitleText, margin, a4Width, a4Height)
% Calculate normalized position for the subtitle
subtitleX = margin / a4Width; % Normalized x position (2 cm margin)
subtitleY = 1 - (margin / a4Height) - 0.1; % Normalized y position (below the title)

% Add subtitle with italic font and dark gray color
annotation(fig, 'textbox', [subtitleX subtitleY 1 - 2 * subtitleX 0.05], 'String', subtitleText, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'FontWeight', 'normal', 'FontAngle', 'italic', ...
    'Color', [0.3 0.3 0.3], 'EdgeColor', 'none', 'Interpreter', 'none'); % Dark gray color
end
