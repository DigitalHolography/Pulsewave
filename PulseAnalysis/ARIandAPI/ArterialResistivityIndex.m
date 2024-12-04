function [] = ArterialResistivityIndex(t, arterial_signal, M0_ff_video, maskArtery)

ToolBox = getGlobalToolBox;
PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);
exportVideos = PW_params.exportVideos;
folder = 'volumeRate';
cArtery = [255 22 18] / 255;
[numX, numY, numFrames] = size(M0_ff_video);

M0_ff_video = rescale(M0_ff_video);
M0_ff_image = mean(M0_ff_video, 3);

arterial_signal = filloutliers(arterial_signal, 'center');
arterial_signal_smooth = smoothdata(arterial_signal, 'rlowess');

%% ARI CALC

min_vRMS = min(arterial_signal_smooth);
max_vRMS = max(arterial_signal_smooth);
ARI = (max_vRMS - min_vRMS) / max_vRMS;

%% ARI FIG

% ARI Graph
graphSignal('6_ARI', folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('ARI = %0.2f', ARI), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, min_vRMS, max_vRMS], yLineLabels = {'', '', ''})

% ARI Map
[hue_ARI, sat_ARI, val_ARI, cmapARI] = ARI2HSVmap(ARI, M0_ff_image, maskArtery);
ARImapRGB = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
ARImapRGB = ARImapRGB .* maskArtery + M0_ff_image .* ~maskArtery;

figure("Visible", "off")
imagesc(ARImapRGB);
title(strcat('Arterial Resistivity Index avg. : ', sprintf(" %3.2f", ARI)));
axis image
axis off
set(gca, 'LineWidth', 2);
fontsize(gca, 14, "points");
c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
c.Label.String = 'Arterial resistivity index';
c.Label.FontSize = 14;
colormap(cmapARI);
exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'ARImapFig.eps')))

if exportVideos
    % ARI GIF
    
    ARIvideoRGB = zeros(numX, numY, 3, numFrames);
    gifWriter = GifWriter("ARI", numFrames);
    f = figure("Visible", "off");
    f.Position = [300, 300, 570, 630];
    
    parfor frameIdx = 1:numFrames
        [hue_ARI, sat_ARI, val_ARI] = ARI2HSVmap(ARI, M0_ff_video(:, :, frameIdx), maskArtery);
        ARIvideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
    end
    
    parfor frameIdx = 1:numFrames
        imagesc(ARIvideoRGB(:, :, :, frameIdx));
        title(strcat('Arterial resistivity index value : ', sprintf(" %3.2f", ARI)));
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 18, "points");
        c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
        c.Label.String = 'Arterial resistivity index';
        c.Label.FontSize = 18;
        colormap(cmapARI);
        frame = getframe(f);
        gifWriter.write(frame, frameIdx);
    end
    
    gifWriter.generate();
    gifWriter.delete();
    
end

%% API CALC

mean_vRMS = mean(arterial_signal_smooth);
API = (max_vRMS - min_vRMS) / mean_vRMS;

%% API FIG

% API Graph
graphSignal('6_API', folder, ...
    t, arterial_signal_smooth, '-', cArtery, ...
    t, arterial_signal, ':', cArtery, ...
    Title = sprintf('API = %0.2f', API), Legend = {'Smooth', 'Raw'}, ...
    yLines = [0, min_vRMS, mean_vRMS, max_vRMS], yLineLabels = {'', '', '', ''})

% API Map
[hue_API, sat_API, val_API, cmapAPI] = API2HSVmap(API, M0_ff_image, maskArtery);
APImapRGB = hsv2rgb(hue_API, sat_API, val_API);
APImapRGB = rescale(APImapRGB .* maskArtery) + M0_ff_image .* ~maskArtery;

figure("Visible", "off")
imagesc(APImapRGB);
title(strcat('Arterial Pulsatility Index avg. : ', sprintf(" %3.2f", API)));
axis image
axis off
set(gca, 'LineWidth', 2);
fontsize(gca, 14, "points");
c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
c.Label.String = 'Arterial pulsatility index';
c.Label.FontSize = 14;
colormap(cmapAPI);
exportgraphics(gca, fullfile(ToolBox.PW_path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, 'APImapFig.png')))
exportgraphics(gca, fullfile(ToolBox.PW_path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, 'APImapFig.eps')))

if exportVideos
    % API GIF
    
    ARIvideoRGB = zeros(numX, numY, 3, numFrames);
    gifWriter = GifWriter("API", numFrames);
    f = figure("Visible", "off");
    f.Position = [300, 300, 570, 630];
    
    tic
    parfor frameIdx = 1:numFrames
        [hue_ARI, sat_ARI, val_ARI] = ARI2HSVmap(ARI, M0_ff_video(:, :, frameIdx), maskArtery);
        ARIvideoRGB(:, :, :, frameIdx) = hsv2rgb(hue_ARI, sat_ARI, val_ARI);
    end
    toc
    tic
    parfor frameIdx = 1:numFrames
        imagesc(ARIvideoRGB(:, :, :, frameIdx));
        title(strcat('Arterial resistivity index value : ', sprintf(" %3.2f", ARI)));
        axis image
        axis off
        set(gca, 'LineWidth', 2);
        fontsize(gca, 18, "points");
        c = colorbar('southoutside', 'Ticks', linspace(0, 1, 6));
        c.Label.String = 'Arterial resistivity index';
        c.Label.FontSize = 18;
        colormap(cmapARI);
        frame = getframe(f);
        gifWriter.write(frame, frameIdx);
    end
    toc
    gifWriter.generate();
    gifWriter.delete();
    
end

end
