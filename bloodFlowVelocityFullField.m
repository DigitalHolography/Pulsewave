function bloodFlowVelocityFullField(v_RMS_all, v_one_cycle, maskArtery, maskVein, videoM0, FreqVideoRGB, ToolBox, path)
%% Velocity funnel Histogram in arteries (exactly the same but with an increasing number of points)
    %FIXME prctile 10% Y = percentil(X,[5 95])
    tic
    PW_params = Parameters_json(path);
    Im = rescale(mean(videoM0, 3));
    [Nx, Ny, N_frame] = size(v_RMS_all);
        ImgM0 = rescale(mean(videoM0, 3));


    %% Init of histogram axis

    %%
    [N, M] = size(maskArtery);
    radius1 = PW_params.velocity_bigRadiusRatio * (M + N) / 2;
    radius2 = PW_params.velocity_smallRadiusRatio * (M + N) / 2;    
    [maskSection] = createMaskSection(ImgM0, maskArtery,radius1,radius2,'_mask_artery_section_velocity_rgb.png', ToolBox, path);
    maskArtery_section = maskArtery & maskSection;

    v_histo_artery = round(v_RMS_all .* maskArtery_section);
    v_min_artery = min(v_histo_artery, [], 'all');
    v_max_artery = max(v_histo_artery, [], 'all');

   
    v_max_all = v_max_artery;
    v_min_all = v_min_artery;

    v_max_all_display = round(0.8 * v_max_all);
    v_min_all_display = round(0.8 * v_min_all);

    yAx = [v_min_all v_max_all];
    yAx_display = yAx;
    

    %% Construct velocity map

    ImgM0 = rescale(mean(videoM0, 3));

    N_interp = size(v_one_cycle, 3);
    v_histo_artery = round(v_one_cycle .* maskArtery);
    v_min = min(v_histo_artery, [], 'all');
    v_max = max(v_histo_artery, [], 'all');

    X = linspace(v_min, v_max, v_max - v_min + 1);
    n = size(X, 2);
    histo = zeros(size(X, 2), N_interp);

    radius1 = PW_params.velocity_bigRadiusRatio * (M + N) / 2;
    radius2 = PW_params.velocity_smallRadiusRatio * (M + N) / 2;    
    
    if PW_params.AllCirclesFlag 
        nbCircles = PW_params.nbCircles;
        radius0 = radius1;
        radiusmid = (radius1+radius2)/2;
        radiusend = radius2;
        deltar = (radiusend-radius0)/nbCircles;
        deltarcentral = (radiusend-radius0)/nbCircles; % two times radius1-radius2 in total
    
        Color_std = [0.7 0.7 0.7];
       
        v_RMS_frame = mean(v_RMS_all,3);
    
        X = linspace(v_min_all, v_max_all, v_max_all - v_min_all + 1);
        histo_artery = zeros(size(X, 2), nbCircles);
        for j = 1:nbCircles % to parforize (change createMaskSection)
    
            r1 = radiusmid+j*deltarcentral/2;
            r2 = radiusmid-j*deltarcentral/2;
            [maskSection] = createMaskSection(ImgM0, maskArtery,r1,r2,sprintf('_mask_artery_section_velocity_growing_sections_%d.png',j), ToolBox, path);
            maskArtery_section = maskArtery & maskSection;
            v_histo_artery = round(v_RMS_frame.*maskArtery_section);
            for xx = 1:Nx
    
                for yy = 1:Ny
                    
                    if (v_histo_artery(xx, yy) ~= 0)
                        i = find(X == v_histo_artery(xx, yy));
                        histo_artery(i, j) = histo_artery(i, j) + 1;
                    end
    
                end
    
            end
            histo_artery(:, j) = histo_artery(:, j)/sum(maskArtery_section,[1,2]);
            non_zero_points = find(maskArtery_section);
            mean_velocity_in_growing_section(j) = mean(v_RMS_frame(non_zero_points),[1,2]);
            std_velocity_in_growing_section(j) = std(v_RMS_frame(non_zero_points));
            growing_gap(j) = r1-r2;
            number_of_points(j)=sum(maskArtery_section,[1,2]);
    
            r1 = radius0+j*deltar;
            r2 = radius0+(j-1)*deltar;
            [maskSection] = createMaskSection(ImgM0, maskArtery,r1,r2,sprintf('_mask_artery_section_velocity_%d.png',j), ToolBox, path);
            maskArtery_section_only = maskArtery & maskSection;
            
            non_zero_points = find(maskArtery_section_only);
            mean_velocity_in_section_only(j) = mean(v_RMS_frame(non_zero_points),[1,2]);
            std_velocity_in_section_only(j) = std(v_RMS_frame(non_zero_points));
    
            rad(j)=r1;
    
    
           
        end
        plot_velocity_funnel_hist = figure(164);
    
        
        %
        xAx = number_of_points;
        
    
        f_distrib_artery = figure(197);
        f_distrib_artery.Position(3:4) = [500 275];
        index_min = find(X == v_min_all_display);
        index_max = find(X == v_max_all_display);
        imagesc(xAx, yAx_display, histo_artery(index_min:index_max, :))
        set(gca, 'YDir', 'normal')
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        colormap("hot")
     
        
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Number of pixels')
        title("Velocity distribution in a growing artery section")
        
        exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'bloodFlowVelocity', sprintf('bloodVelocityinArterieshistogramxNumpoints.png')))

        plot_velocity_funnel = figure(162);
    
        
        curve1 = mean_velocity_in_growing_section + 0.5 * std_velocity_in_growing_section;
        curve2 = mean_velocity_in_growing_section - 0.5 * std_velocity_in_growing_section;
        rad2 = [growing_gap, fliplr(growing_gap)];
        inBetween = [curve1, fliplr(curve2)];
        
        fill(rad2, inBetween, Color_std);
        hold on;
        plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
        plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
        plot(rad, mean_velocity_in_section_only, '-k', 'LineWidth', 2);
        axis tight;
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
    toc
end