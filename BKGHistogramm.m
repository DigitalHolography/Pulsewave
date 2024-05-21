function [] = BKGHistogramm(data_M2M0, maskBKG, ToolBox)
    %% Init of histogram axis

    [size_data_1, size_data_2, size_data_3] = size(data_M2M0);
    pas = 10;

    freq_artery = round(pas * data_M2M0 .* maskBKG);
    freq_min = min(freq_artery(freq_artery ~= 0), [], 'all');
    freq_max = max(freq_artery, [], 'all');
    yAx = [freq_min freq_max];

    %% Velocity Histogram in arteries
    %FIXME prctile 10% Y = percentil(X,[5 95])

    X = linspace(freq_min, freq_max, freq_max - freq_min + 1);
    n = size(X, 2);
    xAx = [0 n * ToolBox.stride / (1000 * ToolBox.fs)];
    histo_BKG = zeros(size(X, 2), size_data_3);
    %histo_video_artery = zeros(size(X,2),size(dataCubeM2M0,3),3,size(dataCubeM2M0,3));

    f_distrib_freq = figure(159);
    f_distrib_freq.Position(3:4) = [1200 550];
    index_min = find(X == freq_min);
    index_max = find(X == freq_max);
    imagesc(xAx, yAx, histo_BKG(index_min:index_max, :))
    set(gca, 'YDir', 'normal')
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    %colormap("hot")
    f = getframe(gcf);
    [M, N, ~] = size(f.cdata);
    histo_video_artery = zeros(M, N, 3, size_data_3);

    %FIXME avoir une ligne à zéro de trois pixel
    %FIXME getframe pour la couleur
    for t = 1:size_data_3

        for x = 1:size_data_1

            for y = 1:size_data_2

                if (freq_artery(x, y, t) ~= 0)
                    i = find(X == freq_artery(x, y, t));
                    histo_BKG(i, t) = histo_BKG(i, t) + 1;
                end

            end

        end

        %histo_video_artery(:,:,t) = flip(histo_artery,1);
        figure(159)
        imagesc(xAx, yAx, histo_BKG(index_min:index_max, :))
        set(gca, 'YDir', 'normal')
        title("Velocity distribution in arteries")
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        f = getframe(gcf);
        histo_video_artery(:, :, :, t) = imresize(f.cdata, [M N]);

    end

    velocity_dist_arteries = frame2im(getframe(gca));

    w = VideoWriter(fullfile(ToolBox.PW_path_avi, strcat(ToolBox.main_foldername, '_velocity_histogram_artery.avi')));
    tmp = mat2gray(histo_video_artery);
    open(w)

    for j = 1:size(histo_video_artery, 4)
        writeVideo(w, tmp(:, :, :, j));
    end

    close(w);

    print('-f157', '-dpng', fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_velocity_distribution_arteries_fullcycle.png')));

    imwrite(velocity_dist_arteries, fullfile(ToolBox.PW_path_png, strcat(ToolBox.main_foldername, '_Velocity_distribution_in_arteries.png')), 'png');
end
