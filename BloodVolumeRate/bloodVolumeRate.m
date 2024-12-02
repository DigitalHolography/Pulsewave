function [] = bloodVolumeRate(ToolBox, maskArtery, maskVein, v_RMS, M0_disp_video, xy_barycenter, flagBloodVelocityProfile)

tic

PW_params = Parameters_json(ToolBox.PW_path);

veins_analysis = PW_params.veins_analysis;
exportVideos = PW_params.exportVideos;
force_width = [];

if ~isempty(PW_params.forcewidth)
    force_width = PW_params.forcewidth;
end

[x_barycenter, y_barycenter] = xy_barycenter{:};

mkdir(ToolBox.PW_path_png, 'volumeRate')
mkdir(ToolBox.PW_path_eps, 'volumeRate')

[numX, numY, numFrames] = size(v_RMS);
[X, Y] = meshgrid(1:numY, 1:numX);
Color_std = [0.7 0.7 0.7];

fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));

v_RMS_AVG = mean(v_RMS, 3);
L = (numY + numX) / 2;

%% 0) Change mask section

radius1 = (PW_params.radius_ratio - PW_params.radius_gap) * (numY + numX) / 2;
radius2 = (PW_params.radius_ratio + PW_params.radius_gap) * (numY + numX) / 2;
cercle_mask1 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= radius1;
cercle_mask2 = sqrt((X - x_barycenter) .^ 2 + (Y - y_barycenter) .^ 2) <= radius2;

maskSection = xor(cercle_mask1, cercle_mask2);

figure(1)
imagesc(maskSection), axis image

%% 1) Find the locations of the sections
%% 1) 1) Arteries

maskSectionArtery = maskSection .* maskArtery;

figure(110)
imagesc(maskSectionArtery .* v_RMS_AVG), axis image

maskSectionArtery = bwlabel(maskSectionArtery);

figure(111)
imagesc(maskSectionArtery), axis image

numSectionsArtery = max(maskSectionArtery, [], 'all');
masksSectionsArtery = zeros(numX, numY, numSectionsArtery);

parfor sectionIdx = 1:numSectionsArtery
    
    masksSectionsArtery(:, :, sectionIdx) = (maskSectionArtery == sectionIdx);
    
end

SubImg_locs_artery = zeros(numSectionsArtery, 2);
SubImg_width_artery = zeros(numSectionsArtery, 1);

for sectionIdx = 1:numSectionsArtery
    [row, col] = find(masksSectionsArtery(:, :, sectionIdx));
    SubImg_locs_artery(sectionIdx, 1) = round(mean(row));
    SubImg_locs_artery(sectionIdx, 2) = round(mean(col));
    SubImg_width_artery(sectionIdx) = 0.01 * numX;
end

%% 1) 2) Veins

if veins_analysis
    
    maskSectionVein = maskSection .* maskVein;
    
    figure(120)
    imagesc(maskSectionVein .* v_RMS_AVG), axis image
    
    maskSectionVein = bwlabel(maskSectionVein);
    
    figure(121), axis image
    imagesc(maskSectionVein)
    
    numSectionsVein = max(maskSectionVein, [], 'all');
    masksSectionsVein = zeros(numX, numY, numSectionsVein);
    
    parfor sectionIdx = 1:numSectionsVein
        masksSectionsVein(:, :, sectionIdx) = (maskSectionVein == sectionIdx);
    end
    
    SubImg_locs_vein = zeros(numSectionsVein, 2);
    SubImg_width_vein = zeros(numSectionsVein, 1);
    
    for sectionIdx = 1:numSectionsVein
        [row, col] = find(masksSectionsVein(:, :, sectionIdx));
        SubImg_locs_vein(sectionIdx, 1) = round(mean(row));
        SubImg_locs_vein(sectionIdx, 2) = round(mean(col));
        SubImg_width_vein(sectionIdx) = 0.01 * numX;
    end
    
end

%% 2) Compute blood volume rate: Cross_section_analysis

strXlabel = 'Time(s)';
strYlabel = 'Velocity (mm.s-1)';

%% Arteries

[avgVolumeRateArtery, stdVolumeRateArtery, crossSectionAreaArtery, avgVelocityArtery, stdVelocityArtery, crossSectionMaskArtery, ~, ~, ~, crossSectionWidthArtery] = crossSectionAnalysis2(ToolBox, SubImg_locs_artery, SubImg_width_artery, maskArtery, v_RMS, PW_params.flowRate_sliceHalfThickness, 'artery', flagBloodVelocityProfile, [], force_width, 1);

labelsArteries = cell(numSectionsArtery, 1);
avgVolumeRateArtery_total = sum(avgVolumeRateArtery, 1);
stdVolumeRateArtery_total = sqrt(sum(stdVolumeRateArtery .^ 2, 1));
avgVelocityArtery_total = sum(avgVelocityArtery, 1);
stdVelocityArtery_total = sqrt(sum(stdVelocityArtery .^ 2, 1));
velocityArtery_mean = mean(avgVelocityArtery_total);

for sectionIdx = 1:numSectionsArtery
    labelsArteries{sectionIdx} = sprintf("A%d", sectionIdx);
end

% Velocity Over Sections
figure(220)
plot(fullTime, avgVelocityArtery, 'LineWidth', 2)
title('Blood velocity in artery sections');
legend(labelsArteries);
fontsize(gca, 14, "points");
xlabel(strXlabel, 'FontSize', 14);
ylabel(strYlabel, 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);
axis tight;

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArteryEverySection.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArteryEverySection.eps')))

% Mean Velocity
velocityArtery_plot = figure(221);
velocityArtery_plot.Position = [200 475 600 300];

graphSignalStd(ToolBox, velocityArtery_plot, avgVelocityArtery_total, stdVelocityArtery_total, numFrames, strYlabel, strXlabel, "Total velocity in arteries", 'mm/s')

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArterySection.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityArterySection.eps')))

plot2txt(ToolBox, fullTime, avgVelocityArtery_total, 'AVGVelocityArteriesTotal')
plot2txt(ToolBox, fullTime, stdVelocityArtery_total, 'STDVelocityArteriesTotal')
plot2txt(ToolBox, fullTime, avgVolumeRateArtery_total, 'AVGVolumeRateArteriesTotal')
plot2txt(ToolBox, fullTime, stdVolumeRateArtery_total, 'STDVolumeRateArteriesTotal')

for sectionIdx = 1:numSectionsArtery
    plot2txt(ToolBox, fullTime, avgVolumeRateArtery(sectionIdx, :), strcat('volumeRate_artery_A', num2str(sectionIdx)))
    plot2txt(ToolBox, fullTime, stdVolumeRateArtery(sectionIdx, :), strcat('volumeRate_artery_std_A', num2str(sectionIdx)))
    plot2txt(ToolBox, fullTime, avgVelocityArtery(sectionIdx, :), strcat('avg_velocity_artery_A', num2str(sectionIdx)))
    plot2txt(ToolBox, fullTime, stdVelocityArtery(sectionIdx, :), strcat('std_velocity_artery_A', num2str(sectionIdx)))
end

%% 2) 2) Veins
if veins_analysis
    [avgVolumeRateVein, stdVolumeRateVein, crossSectionAreaVein, avgVelocityVein, stdVelocityVein, crossSectionMaskVein, ~, ~, ~, crossSectionWidthVein] = crossSectionAnalysis2(ToolBox, SubImg_locs_vein, SubImg_width_vein, maskVein, v_RMS, PW_params.flowRate_sliceHalfThickness, 'vein', flagBloodVelocityProfile, [], force_width, 1);
    
    labelsVeins = cell(numSectionsVein, 1);
    avgVolumeRateVein_total = sum(avgVolumeRateVein, 1);
    stdVolumeRateVein_total = sqrt(sum(stdVolumeRateVein .^ 2, 1));
    avgVelocityVein_total = sum(avgVelocityVein, 1);
    stdVelocityVein_total = sqrt(sum(stdVelocityVein .^ 2, 1));
    velocityVein_mean = mean(avgVelocityVein_total);
    
    for sectionIdx = 1:numSectionsVein
        labelsVeins{sectionIdx} = sprintf("V%d", sectionIdx);
    end
    
    % Velocity Over Sections
    figure(220)
    plot(fullTime, avgVelocityVein, 'LineWidth', 2)
    title('Blood velocity in vein sections');
    legend(labelsVeins);
    fontsize(gca, 14, "points");
    xlabel(strXlabel, 'FontSize', 14);
    ylabel(strYlabel, 'FontSize', 14);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2);
    axis tight;
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinEverySection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinEverySection.eps')))
    
    % Mean Velocity
    velocityVein_plot = figure(221);
    velocityVein_plot.Position = [200 475 600 300];
    
    graphSignalStd(ToolBox, velocityVein_plot, avgVelocityVein_total, stdVelocityVein_total, numFrames, strYlabel, strXlabel, "Total velocity in arteries", 'mm/s')
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinSection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'velocityVeinSection.eps')))
    
    plot2txt(ToolBox, fullTime, avgVelocityVein_total, 'AVGVelocityVeinsTotal')
    plot2txt(ToolBox, fullTime, stdVelocityVein_total, 'STDVelocityVeinsTotal')
    plot2txt(ToolBox, fullTime, avgVolumeRateVein_total, 'AVGVolumeRateVeinsTotal')
    plot2txt(ToolBox, fullTime, stdVolumeRateVein_total, 'STDVolumeRateVeinsTotal')
    
    for sectionIdx = 1:numSectionsVein
        plot2txt(ToolBox, fullTime, avgVolumeRateVein(sectionIdx, :), strcat('volumeRate_vein_V', num2str(sectionIdx)))
        plot2txt(ToolBox, fullTime, stdVolumeRateVein(sectionIdx, :), strcat('volumeRate_vein_std_V', num2str(sectionIdx)))
        plot2txt(ToolBox, fullTime, avgVelocityVein(sectionIdx, :), strcat('avg_velocity_vein_V', num2str(sectionIdx)))
        plot2txt(ToolBox, fullTime, stdVelocityVein(sectionIdx, :), strcat('std_velocity_vein_V', num2str(sectionIdx)))
    end
    
end

%% 3) Vein and artery numerotation
M0_disp_image = rescale(mean(M0_disp_video, 3));
numerotation_plot = figure(300);
numerotation_plot.Position = [100 100 600 600];
x_center = x_barycenter;
y_center = y_barycenter;
graphMaskTags(300, M0_disp_image, crossSectionMaskArtery, SubImg_locs_artery, labelsArteries, x_center, y_center)
exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'arteries_numerotation.png')))

if veins_analysis
    numerotationv_plot = figure(301);
    numerotationv_plot.Position = [100 100 600 600];
    graphMaskTags(301, M0_disp_image, crossSectionMaskVein, SubImg_locs_vein, labelsVeins, x_center, y_center, Color = [0 0 1]);
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'veins_numerotation.png')));
    graphMaskTags(301, M0_disp_image, crossSectionMaskArtery, SubImg_locs_artery, labelsArteries, x_center, y_center);
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'arteries_veins_numerotation.png')));
end

%% 4) Volume Rate FigurePosition

%% 4) 1) Video Generation

%% 4) 2) Average BVR figure and Section widths

volume_rate_plot = figure(420);
volume_rate_plot.Position = [200 200 600 600];
etiquettes_frame_values = round(mean(avgVolumeRateArtery(:, :), 2), 1);
graphMaskTags(volume_rate_plot, M0_disp_image, crossSectionMaskArtery, SubImg_locs_artery, etiquettes_frame_values, x_center, y_center);
title(sprintf("%s : %02.0f %s", 'Average total blood volume rate in arteries', round(mean(avgVolumeRateArtery_total), 1), 'µL/min'));
set(gca, 'FontSize', 14)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArteryImage.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArteryImage.eps')))

section_width_plot = figure(430);
section_width_plot.Position = [200 200 600 600];
etiquettes_frame_values = append(string(round(crossSectionWidthArtery * PW_params.cropSection_pixelSize / (2 ^ PW_params.k) * 1000, 1)), "µm");
graphMaskTags(section_width_plot, M0_disp_image, crossSectionMaskArtery, SubImg_locs_artery, etiquettes_frame_values, x_center, y_center, Fontsize = 12);
title(sprintf("%s", 'Cross section width in arteries (µm)'));
set(gca, 'FontSize', 14)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'crossSectionWidthArteryImage.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'crossSectionWidthArteryImage.eps')))

if veins_analysis
    
    volume_rate_plot = figure(520);
    volume_rate_plot.Position = [200 200 600 600];
    etiquettes_frame_values = round(mean(avgVolumeRateVein(:, :), 2), 1);
    graphMaskTags(volume_rate_plot, M0_disp_image, crossSectionMaskVein, SubImg_locs_vein, etiquettes_frame_values, x_center, y_center, Color = [0 0 1]);
    title(sprintf("%s : %02.0f %s", 'Average total blood volume rate in veins', round(mean(avgVolumeRateVein_total), 1), 'µL/min'));
    set(gca, 'FontSize', 14)
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinImage.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinImage.eps')))
    
    section_width_plot = figure(530);
    section_width_plot.Position = [200 200 600 600];
    etiquettes_frame_values = append(string(round(crossSectionWidthVein * PW_params.cropSection_pixelSize / (2 ^ PW_params.k) * 1000, 1)), "µm");
    graphMaskTags(section_width_plot, M0_disp_image, crossSectionMaskVein, SubImg_locs_vein, etiquettes_frame_values, x_center, y_center, Color = [0 0 1], Fontsize = 12);
    title(sprintf("%s", 'Cross section width in veins (µm)'));
    set(gca, 'FontSize', 14)
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'crossSectionWidthVeinImage.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'crossSectionWidthVeinImage.eps')))
    
end

%% 5) Plot Volume Rate
%% 5) 1) Artery plot

volumeRateArtery_plot = figure(510);
volumeRateArtery_plot.Position = [200 200 600 300];

graphSignalStd(ToolBox, volumeRateArtery_plot, avgVolumeRateArtery_total, stdVolumeRateArtery_total, numFrames, 'Blood volume rate (µL/min)', 'Time (s)', "Total blood volume rate in arteries", "µL/min", ylimm = [0 max(avgVolumeRateArtery_total)])

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArterySection.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateArterySection.eps')))

plot2txt(ToolBox, fullTime, avgVolumeRateArtery_total, 'AVGVolumeRateArteryTotal')
plot2txt(ToolBox, fullTime, stdVolumeRateArtery_total, 'STDVolumeRateArteryTotal')

%% 5) 2) Vein plot

if veins_analysis
    
    volumeRateVein_plot = figure(520);
    volumeRateVein_plot.Position = [200 475 600 300];
    
    graphSignalStd(ToolBox, volumeRateVein_plot, avgVolumeRateVein_total, stdVolumeRateVein_total, numFrames, 'Blood volume rate (µL/min)', 'Time (s)', "Total blood volume rate in veins", "µL/min", ylimm = [0 max(avgVolumeRateArtery_total)])
    
    exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinSection.png')))
    exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateVeinSection.eps')))
    
    plot2txt(ToolBox, fullTime, avgVolumeRateArtery_total, 'AVGVolumeRateVeinTotal')
    plot2txt(ToolBox, fullTime, stdVolumeRateArtery_total, 'STDVolumeRateVeinTotal')
    
end

%% 5) 3) Artery Plot & Vein Plot progression

graphCombined(ToolBox, M0_disp_video, crossSectionMaskArtery, SubImg_locs_artery, avgVolumeRateArtery, avgVolumeRateArtery_total, stdVolumeRateArtery_total, xy_barycenter, [], 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in arteries', 'µL/min', skip = ~exportVideos);

if veins_analysis
    
    graphCombined(ToolBox, M0_disp_video, crossSectionMaskVein, SubImg_locs_vein, avgVolumeRateVein, avgVolumeRateVein_total, stdVolumeRateVein_total, xy_barycenter, [], 'Blood Volume Rate (µL/min)', 'Time (s)', 'Total Blood Volume Rate in veins', 'µL/min', skip = ~exportVideos, Color = [0 0 1]);
    
end

%% txt file output with measured pulse wave parameters
volumeRateArtery_mean = mean(avgVolumeRateArtery_total);

if veins_analysis
    volumeRateVein_mean = mean(avgVolumeRateVein_total);
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, [ ...
        'Value of total arterial blood volume rate (µL/min) :\n%d\n' ...
        'Value of total venous blood volume rate (µL/min) :\n%d\n'], ...
        volumeRateArtery_mean, ...
        volumeRateVein_mean);
    fclose(fileID);
else
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, ...
        'Value of total arterial blood volume rate (µL/min) :\n%f\n', ...
        volumeRateArtery_mean);
    fclose(fileID);
end

for sectionIdx = 1:numSectionsArtery
    fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
    fprintf(fileID, [ ...
        'Artery n°A%d : cross_section (mm^2) : \n %d \n' ...
        'Artery n°A%d : vessel diameter (µm) : \n %d \n' ...
        'Artery n°A%d : average velocity (mm/s) : \n %d \n' ...
        'Artery n°A%d : blood volume rate (µL/min) : \n %d \n \n'], ...
        sectionIdx, ...
        crossSectionAreaArtery(sectionIdx), ...
        sectionIdx, ...
        2 * sqrt(crossSectionAreaArtery(sectionIdx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
        sectionIdx, ...
        avgVelocityArtery(sectionIdx), ...
        sectionIdx, ...
        avgVolumeRateArtery(sectionIdx)); % mm^3/s -> µL/min
    fclose(fileID);
end

if veins_analysis
    
    for sectionIdx = 1:numSectionsVein
        fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_pulseWaveOutputParameters.txt')), 'a');
        fprintf(fileID, [ ...
            'Vein n°V%d : cross_section (mm^2) : \n %d \n ' ...
            'Vein n°V%d : vessel diameter (µm) : \n %d \n ' ...
            'Vein n°V%d : average velocity (mm/s) : \n %d \n ' ...
            'Vein n°V%d : blood volume rate (µL/min) : \n %d \n \n'], ...
            sectionIdx, ...
            crossSectionAreaVein(sectionIdx), ...
            sectionIdx, ...
            2 * sqrt(crossSectionAreaVein(sectionIdx) / pi) * 1000, ... % calculation of the diameter knowing the disc area
            sectionIdx, ...
            avgVelocityVein(sectionIdx), ...
            sectionIdx, ...
            avgVolumeRateVein(sectionIdx)); % mm^3/s -> µL/min
        fclose(fileID);
    end
    
end

close all

%% 6) Arterial Resistivity with Volume Rate
maxVolumeRate = max(avgVolumeRateArtery_total(:));
minVolumeRate = min(avgVolumeRateArtery_total(:));
meanVolumeRate = mean(avgVolumeRateArtery_total(:));

%% 6) 1) Arterial Resisitivity Index

ARI = (maxVolumeRate - minVolumeRate) / maxVolumeRate;

figure(610)

curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

hold on

fill(fullTime2, inBetween, Color_std);
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
yline(volumeRateArtery_mean, '--k', 'LineWidth', 2)
yline(maxVolumeRate, '--r', 'Linewidth', 2)
yline(minVolumeRate, '--r', 'Linewidth', 2)
axis tight;

volumeRateArtery_ax = axis;
volumeRateArtery_ax(3) = 0;
hold off

ylabel('Blood volume rate (µL/min)')
xlabel('Time (s)')
title(sprintf("Arterial Resistivity Index %0.2f ", ARI))
axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
fontsize(gca, 14, "points");
box on
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateARI.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateARI.eps')))

%% 6) 2) Arterial Pulsatility Index

API = (maxVolumeRate - minVolumeRate) / meanVolumeRate;

figure(620)

curve1 = avgVolumeRateArtery_total + 0.5 * stdVolumeRateArtery_total;
curve2 = avgVolumeRateArtery_total - 0.5 * stdVolumeRateArtery_total;
fullTime2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

hold on

fill(fullTime2, inBetween, Color_std);
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, avgVolumeRateArtery_total, '-k', 'LineWidth', 2);
yline(maxVolumeRate, '--r', 'Linewidth', 2)
yline(meanVolumeRate, '--r', 'Linewidth', 2)
yline(minVolumeRate, '--r', 'Linewidth', 2)
axis tight;

volumeRateArtery_ax = axis;
volumeRateArtery_ax(3) = 0;
hold off

ylabel('Blood volume rate (µL/min)')
xlabel('Time (s)')
title(sprintf("Arterial Pulsatility Index %0.2f ", API))
axis([volumeRateArtery_ax(1) volumeRateArtery_ax(2) volumeRateArtery_ax(3) volumeRateArtery_ax(4)]);
fontsize(gca, 14, "points");
box on
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateAPI.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'volumeRateAPI.eps')))

close all

end
