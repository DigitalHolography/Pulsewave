function generateHealthReport(pdfPath, dataFilePath, signal)
    % Set A4 dimensions in centimeters
    a4Width = 21.0; % A4 width in cm
    a4Height = 29.7; % A4 height in cm

    % Define margins (2 cm on all sides)
    margin = 2; % in cm

    % Create a new figure for the report
    fig = figure('Units', 'centimeters', 'Position', [0 0 a4Width a4Height], 'Visible', 'on');

    % Read and parse the data file
    data = parseDataFile(dataFilePath);

    % Add title
    addTitle(fig, 'Health Report for Practitioners', margin, a4Width, a4Height);

    % Add data fields
    yPos = 0.9 - (margin / a4Height); % Starting vertical position for annotations (normalized units)
    yPos = addField(fig, sprintf('Heart Beat: %.2f bpm', data.heartBeat), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Systole Indices: %s', mat2str(data.systoleIndices)), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Number of Cycles: %d', data.numCycles), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Max Systole Indices: %s', mat2str(data.maxSystoleIndices)), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Min Systole Indices: %s', mat2str(data.minSystoleIndices)), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Mean Blood Volume Rate (Artery): %.2f µL/min', data.meanBloodVolumeRateArtery), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Std Blood Volume Rate (Artery): %.2f µL/min', data.stdBloodVolumeRateArtery), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Mean Blood Volume Rate (Vein): %.2f µL/min', data.meanBloodVolumeRateVein), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Std Blood Volume Rate (Vein): %.2f µL/min', data.stdBloodVolumeRateVein), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Max Systole Blood Volume Rate (Artery): %.2f µL/min', data.maxSystoleBloodVolumeRateArtery), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Min Diastole Blood Volume Rate (Artery): %.2f µL/min', data.minDiastoleBloodVolumeRateArtery), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Stroke Volume (Artery): %.2f nL', data.strokeVolumeArtery), yPos, margin, a4Width);
    yPos = addField(fig, sprintf('Total Volume (Artery): %.2f nL', data.totalVolumeArtery), yPos, margin, a4Width);

    % Add signal plot
    addSignalPlot(fig, signal, yPos, margin, a4Width, a4Height);

    % Save the figure as a PDF
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');

    % Close the figure
    % close(fig);
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
    
    data.meanBloodVolumeRateArtery = str2double(regexp(fileContent, 'Mean Blood Volume Rate Artery : ([\d.]+)', 'tokens', 'once'));
    data.stdBloodVolumeRateArtery = str2double(regexp(fileContent, 'Std Blood Volume Rate Artery : ([\d.]+)', 'tokens', 'once'));
    data.meanBloodVolumeRateVein = str2double(regexp(fileContent, 'Mean Blood Volume Rate Vein : ([\d.]+)', 'tokens', 'once'));
    data.stdBloodVolumeRateVein = str2double(regexp(fileContent, 'Std Blood Volume Rate Vein : ([\d.]+)', 'tokens', 'once'));
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
    
    annotation(fig, 'textbox', [titleX titleY-0.1 1-2*titleX 0.1], 'String', titleText, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none');
end

% Helper function to add a data field
function yPos = addField(fig, fieldText, yPos, margin, a4Width)
    % Calculate normalized position for the field
    fieldX = margin / a4Width; % Normalized x position (2 cm margin)
    fieldY = yPos ; % Normalized y position (2 cm margin)
    
    annotation(fig, 'textbox', [fieldX fieldY 1-2*fieldX 0.05], 'String', fieldText, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 12, 'EdgeColor', 'none');
    yPos = yPos - 0.02; % Move down for the next field
end

% Helper function to add the signal plot
function addSignalPlot(fig, signal, yPos, margin, a4Width, a4Height)
    % Calculate normalized position for the plot
    plotX = margin / a4Width; % Normalized x position (2 cm margin)
    plotY = yPos - 0.3 - (margin / a4Height); % Normalized y position (2 cm margin)
    plotWidth = 1 - 2 * plotX; % Normalized width
    plotHeight = 0.25; % Normalized height
    
    axes('Parent', fig, 'Position', [plotX plotY plotWidth plotHeight]); % Position the plot
    plot(signal, 'r', 'LineWidth', 2);
    title('Arterial Signal');
    xlabel('Time');
    ylabel('Amplitude');

    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);

    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), 0, axP(4)*1.07])
    box on
end