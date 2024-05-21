function [] = spectrogram(maskArtery, maskBackground, SH_cube, ToolBox)
    %% Variables
    cubeSize = size(SH_cube, 1);
    cubeFreqLength = size(SH_cube, 3);
    cubeFrameLength = size(SH_cube, 4);
    f1 = ToolBox.f1;
    f2 = ToolBox.f2;
    fs = ToolBox.fs;
    batchStride = ToolBox.stride; %convert Hz to s
    % interpolation parameter
    k_int = 4;
    minimumTreshold = 10000; % for the log plotting, we change all smallest value to greater one. In order to avoid log(0).

    %% Resize masks

    ratio = cubeSize / size(maskArtery, 1); % to resize the masks

    maskArtery = logical(imresize(maskArtery, ratio));
    maskBackground = logical(imresize(maskBackground, ratio));

    clear ratio

    %% Video
    video = VideoWriter("stdTemp.avi");
    open(video)

    for i = 1:cubeFrameLength
        standardDev = std(reshape(SH_cube(:, :, :, i), cubeSize, cubeSize, cubeFreqLength), 0, 3);
        figure(1)
        imagesc(standardDev)
        writeVideo(video, getframe(figure(1)))

        if mod(i, 20) == 0
            disp(['Frame : ', num2str(i)])
        end

    end

    close(video)

    video = VideoWriter("stdFreq.avi");
    open(video)

    for i = 1:cubeFreqLength
        standardDev = std(reshape(SH_cube(:, :, i, :), cubeSize, cubeSize, cubeFrameLength), 0, 3);
        figure(1)
        imagesc(standardDev)
        writeVideo(video, getframe(figure(1)))

        if mod(i, 20) == 0
            disp(['Frame : ', num2str(i)])
        end

    end

    close(video)
    clear standardDev

    %% Display
    % Avg Spectrogram plot
    [ArterySpectrum, ~, DeltaSpectrum] = createSpectrum(maskArtery, maskBackground, SH_cube(:, :, :, 1), ToolBox);

    A = zeros((cubeFreqLength - f1) * k_int, cubeFrameLength);
    BG = zeros((cubeFreqLength - f1) * k_int, cubeFrameLength);
    DELTA = zeros(cubeFreqLength * k_int, cubeFrameLength);

    A = zeros(cubeFreqLength * 2 * k_int, cubeFrameLength);
    BG = zeros(cubeFreqLength * 2 * k_int, cubeFrameLength);
    DELTA = zeros(cubeFreqLength * 2 * k_int, cubeFrameLength);

    DELTA_Smooth = zeros(cubeFreqLength, cubeFrameLength);
    A_Smooth = zeros(cubeFreqLength, cubeFrameLength);
    BG_Smooth = zeros(cubeFreqLength, cubeFrameLength);

    frq_shift = linspace(-fs / 2, fs / 2, cubeFreqLength * 2 * k_int)';

    % Figure
    fifig = figure('Name', 'Video Spectrum');
    mamax = 2 * max(ArterySpectrum);
    deltaFirstLimit = DeltaSpectrum(1:round(0.25 * size(DeltaSpectrum, 2)));
    mimin = 0.8 * min(deltaFirstLimit(deltaFirstLimit > 1));

    % Video Creation
    video = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_spectrumDim4')));

    open(video)
    freq = linspace(1, fs / 2 - fs / cubeFreqLength, cubeFreqLength / 2);
    AVG_M0_artery = zeros(1, cubeFrameLength);
    AVG_M0_background = zeros(1, cubeFrameLength);
    AVG_M1M0_artery = zeros(1, cubeFrameLength);
    AVG_M1M0_background = zeros(1, cubeFrameLength);
    AVG_M2M0_artery = zeros(1, cubeFrameLength);
    AVG_M2M0_background = zeros(1, cubeFrameLength);

    for i = 1:cubeFrameLength
        %Spectro
        [ArterySpectrum, BgSpectrum, DeltaSpectrum, ~, Lorenz_Arteries, Lorenz_BKG, M0_artery, M0_background, M1M0_artery, M1M0_background, M2M0_artery, M2M0_background] = createSpectrum(maskArtery, maskBackground, SH_cube(:, :, :, i), ToolBox);
        AVG_M0_artery(i) = sum(M0_artery, "all") / nnz(M0_artery);
        AVG_M0_background(i) = sum(M0_background, "all") / nnz(M0_background);
        AVG_M1M0_artery(i) = sum(M1M0_artery, "all") / nnz(M1M0_artery);
        AVG_M1M0_background(i) = sum(M1M0_background, "all") / nnz(M1M0_background);
        AVG_M2M0_artery(i) = sum(M2M0_artery, "all") / nnz(M2M0_artery);
        AVG_M2M0_background(i) = sum(M2M0_background, "all") / nnz(M2M0_background);

        figure(92)
        loglog(freq, Lorenz_Arteries)
        hold on
        loglog(freq, Lorenz_BKG)
        hold off

        figure(93)
        %smooth
        DELTA_Smooth(:, i) = smoothdata(DeltaSpectrum);
        A_Smooth(:, i) = smoothdata(ArterySpectrum);
        BG_Smooth(:, i) = smoothdata(BgSpectrum);

        %interp00
        tmp_DELTA = interp(DELTA_Smooth(:, i), 2 * k_int);
        tmp_A = interp(A_Smooth(:, i), 2 * k_int);
        tmp_BG = interp(BG_Smooth(:, i), 2 * k_int);
        tmp_DELTA(tmp_DELTA < minimumTreshold) = min(tmp_DELTA(tmp_DELTA > minimumTreshold), [], 'all');

        % disp step
        if mod(i, 20) == 0
            disp(['Frame : ', num2str(i)])
        end

        % PLot
        semilogy(frq_shift, tmp_A, ...
            frq_shift, tmp_BG, ...
            frq_shift, tmp_DELTA);
        ylim([mimin mamax])
        xlim([-fs / 2 fs / 2])
        title('Average spectrum')
        xlabel('frequency (kHz)', 'FontSize', 14);
        ylabel('A.U.', 'FontSize', 14);
        legend('Artery', 'Bg', 'Delta')
        writeVideo(video, getframe(fifig));
        %
        %     DELTA(:,i) = tmp_DELTA(cubeFreqLength*k_int+1:cubeFreqLength*2*k_int,1);
        %     A(:,i)= tmp_A((cubeFreqLength+f1)*k_int+1:cubeFreqLength*2*k_int,1);
        %     BG(:,i)= tmp_BG((cubeFreqLength+f1)*k_int+1:cubeFreqLength*2*k_int,1);

        DELTA(:, i) = tmp_DELTA(:, 1);
        A(:, i) = tmp_A(:, 1);
        BG(:, i) = tmp_BG(:, 1);

    end

    close(video)
    disp('Video creation ok.')
    %clear tmp_BG tmp_A tmp_DELTA DELTA_Smooth A_Smooth BG_Smooth ArterySpectrum BgSpectrum DeltaSpectrum
    %% Figure spectrogram avg

    % DELTA
    yAx = [0 f2];
    xAx = [0 cubeFrameLength * batchStride / (1000 * fs)];

    XAx = linspace(0, xAx(2), size(DELTA, 2));
    YAx = linspace(-yAx(2), yAx(2), size(DELTA, 1));
    figure(38)
    imagesc(XAx, YAx, log(DELTA))
    set(gca, 'YDir', 'normal')
    %ylim([0 f2])
    ylabel('Frequency (kHz)')
    xlabel('Time (s)')
    title('Delta Spectrogram')
    colormap(gray .^ 3)
    pbaspect([cubeFrameLength / (cubeFreqLength * k_int) 1 1]);
    hold on
    p = plot(XAx, mean(log(DELTA), 1));
    p.LineWidth = 2;
    p.Color = "red";
    hold on
    q = plot(XAx, AVG_M1M0_artery - AVG_M1M0_background);
    q.LineWidth = 2;
    q.Color = "magenta";
    hold on
    r = plot(XAx, AVG_M2M0_artery - AVG_M2M0_background);
    r.LineWidth = 2;
    r.Color = "yellow";
    legend('AVG', 'M1/M0', '(M2/M0)^(1/2)')
    hold off

    hold off

    % Artery

    figure(39)
    imagesc(XAx, YAx, log(abs(A + 0.001)))
    set(gca, 'YDir', 'normal')
    %ylim([0 f2])
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    title('artery Spectrogram')
    colormap(gray .^ 3)
    pbaspect([cubeFrameLength / (cubeFreqLength * k_int) 1 1]);
    hold on
    p = plot(XAx, mean(log(A), 1));
    p.LineWidth = 2;
    p.Color = "red";
    hold on
    q = plot(XAx, AVG_M1M0_artery);
    q.LineWidth = 2;
    q.Color = "magenta";
    hold on
    r = plot(XAx, AVG_M2M0_artery);
    r.LineWidth = 2;
    r.Color = "yellow";
    legend('AVG', 'M1/M0', '(M2/M0)^(1/2)')
    hold off

    % Background

    figure(40)
    imagesc(XAx, YAx, log(abs(BG + 0.001)))
    set(gca, 'YDir', 'normal')
    %ylim([0 f2])
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    title('background Spectrogram')
    colormap(gray .^ 3)
    pbaspect([cubeFrameLength / (cubeFreqLength * k_int) 1 1]);
    hold on
    p = plot(XAx, mean(log(BG), 1));
    p.LineWidth = 2;
    p.Color = "red";
    q = plot(XAx, AVG_M1M0_background);
    q.LineWidth = 2;
    q.Color = "magenta";
    hold on
    r = plot(XAx, AVG_M2M0_background);
    r.LineWidth = 2;
    r.Color = "yellow";
    legend('AVG', 'M1/M0', '(M2/M0)^(1/2)')
    hold off

    print('-f38', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_deltaSpectrogram.png')));
    print('-f39', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_arterySpectrogram.png')));
    print('-f40', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_backgroundSpectrogram.png')));
end
