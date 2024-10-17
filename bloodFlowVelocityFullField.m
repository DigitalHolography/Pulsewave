function bloodFlowVelocityFullField(v_RMS_video, ~, maskArtery, maskVein, M0_disp_video, ToolBox, path)
    %% Velocity funnel Histogram in arteries (exactly the same but with an increasing number of points)
    %FIXME prctile 10% Y = percentil(X,[5 95])

    PW_params = Parameters_json(path);
    veinsAnalysis = PW_params.veins_analysis;

    tic

    ImgM0 = rescale(mean(M0_disp_video, 3));
    [numX, numY, ~] = size(v_RMS_video);

    %% Init of histogram axis

    %%
    radius1 = PW_params.velocityBigRadiusRatio * (numY + numX) / 2;
    radius2 = PW_params.velocitySmallRadiusRatio * (numY + numX) / 2;
    [maskSection] = createMaskSection(ImgM0, maskArtery, radius1, radius2, '_mask_artery_section_velocity_rgb.png', ToolBox, path);
    maskArtery_section = maskArtery & maskSection;

    v_histoArtery = round(v_RMS_video .* maskArtery_section);
    v_minArtery = min(v_histoArtery, [], 'all');
    v_maxArtery = max(v_histoArtery, [], 'all');

    if veinsAnalysis
        v_histoVein = round(v_RMS_video .* maskVein);
        v_minVein = min(v_histoVein, [], 'all');
        v_maxVein = max(v_histoVein, [], 'all');
        v_maxAll = max(v_maxArtery, v_maxVein);
        v_minAll = min(v_minArtery, v_minVein);
    else
        v_maxAll = v_maxArtery;
        v_minAll = v_minArtery;
    end

    v_maxAll_display = round(0.8 * v_maxAll);
    v_minAll_display = round(0.8 * v_minAll);

    yAx = [v_minAll v_maxAll];
    yAx_display = yAx;

    %% Construct velocity map

    radius1 = PW_params.velocityBigRadiusRatio * (numY + numX) / 2;
    radius2 = PW_params.velocitySmallRadiusRatio * (numY + numX) / 2;

    if PW_params.AllCirclesFlag
        nbCircles = PW_params.nbCircles;
        radius0 = radius1;
        radiusmid = (radius1 + radius2) / 2;
        radiusend = radius2;
        deltar = (radiusend - radius0) / nbCircles;
        deltarcentral = (radiusend - radius0) / nbCircles; % two times radius1-radius2 in total

        Color_std = [0.7 0.7 0.7];

        v_RMS_frame = mean(v_RMS_video, 3);

        X = linspace(v_minAll, v_maxAll, v_maxAll - v_minAll + 1);
        histo_artery = zeros(size(X, 2), nbCircles);

        for jj = 1:nbCircles % to parforize (change createMaskSection)
            radiusmid = radius0 + jj * deltar;

            for j = 1:nbCircles

                r1 = radiusmid + j * deltarcentral / 2;
                r2 = radiusmid - j * deltarcentral / 2;
                [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, sprintf('_mask_artery_section_velocity_growing_sections_%d.png', j), ToolBox, path);
                maskArtery_section = maskArtery & maskSection;
                non_zero_points = find(maskArtery_section);
                mean_velocity_in_growing_section(jj, j) = mean(v_RMS_frame(non_zero_points), [1, 2]);
                std_velocity_in_growing_section(jj, j) = std(v_RMS_frame(non_zero_points));
                growing_gap(jj, j) = r2 - r1;
            end

        end

        radiusmid = (radius1 + radius2) / 2;

        for j = 1:nbCircles % to parforize (change createMaskSection)

            r1 = radiusmid + j * deltarcentral / 2;
            r2 = radiusmid - j * deltarcentral / 2;
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, sprintf('_mask_artery_section_velocity_growing_sections_%d.png', j), ToolBox, path);
            maskArtery_section = maskArtery & maskSection;
            v_histoArtery = round(v_RMS_frame .* maskArtery_section);

            for xx = 1:numX

                for yy = 1:numY

                    if (v_histoArtery(xx, yy) ~= 0)
                        i = find(X == v_histoArtery(xx, yy));
                        histo_artery(i, j) = histo_artery(i, j) + 1;
                    end

                end

            end

            histo_artery(:, j) = histo_artery(:, j) / sum(maskArtery_section, [1, 2]);

            number_of_points(j) = sum(maskArtery_section, [1, 2]);

            r1 = radius0 + j * deltar;
            r2 = radius0 + (j - 1) * deltar;
            [maskSection] = createMaskSection(ImgM0, maskArtery, r1, r2, sprintf('_mask_artery_section_velocity_%d.png', j), ToolBox, path);
            maskArtery_section_only = maskArtery & maskSection;

            non_zero_points = find(maskArtery_section_only);
            mean_velocity_in_section_only(j) = mean(v_RMS_frame(non_zero_points), [1, 2]);
            std_velocity_in_section_only(j) = std(v_RMS_frame(non_zero_points));

            rad(j) = r1;

        end

        plot_velocity_funnel_hist = figure(164);

        %
        xAx = number_of_points;

        f_distrib_artery = figure(197);
        f_distrib_artery.Position(3:4) = [500 275];
        index_min = find(X == v_minAll_display);
        index_max = find(X == v_maxAll_display);
        imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
        set(gca, 'YDir', 'normal')
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        colormap("hot")

        ylabel('Velocity (mm.s^{-1})')
        xlabel('Number of pixels')
        title("Velocity distribution in a growing artery section")

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArterieshistogramxNumpoints.png')))

        plot_velocity_funnel = figure(162);

        for i = 1:nbCircles
            plot(growing_gap(i, :), mean_velocity_in_growing_section(i, :), 'LineWidth', 2);
            hold on
        end

        axis tight;
        aa = axis;
        aa(3) = 0;
        aa(4) = 1.3 * aa(4);
        axis(aa);
        hold off

        ylabel('Velocity (mm.s^{-1})')
        xlabel('radius gap (px)')
        title("Velocity in growing arteries sections with the growing section")
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArteriesxradiusgap.png')))

        plot_velocity_in_sections = figure(161);

        curve1 = mean_velocity_in_section_only + 0.5 * std_velocity_in_section_only;
        curve2 = mean_velocity_in_section_only - 0.5 * std_velocity_in_section_only;
        rad2 = [rad, fliplr(rad)];
        inBetween = [curve1, fliplr(curve2)];

        fill(rad2, inBetween, Color_std);
        hold on;
        plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(rad, mean_velocity_in_section_only, '-k', 'LineWidth', 2);
        axis tight;
        hold off

        ylabel('Velocity (mm.s^{-1})')
        xlabel('radius (px)')
        title("Velocity in arteries sections with the radius to center")
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])

        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArteriesxradius.png')))
    end

    %close all

fprintf("- Blood Flow Velocity over the full field timing : %d\n", round(toc))

end
