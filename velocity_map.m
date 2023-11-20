function [] = velocity_map(maskArtery, maskVein, v_RMS, ToolBox)

% TRUE MIN and MAX V_RMS but not realistic

Im = mat2gray(squeeze(mean(v_RMS,3)));

% Ones = ones(size(v_RMS));
% V = mean(v_RMS,3);
% Vmax_Arteries = max(V.*maskArtery,[],'all');
% Vmax_Veins = max(V.*maskVein,[],'all');
% Vmin_Arteries = min(V.*maskArtery+Vmax_Arteries*Ones.*(~maskArtery),[],'all');
% Vmin_Veins = min(V.*maskVein+Vmax_Veins*Ones.*(~maskVein),[],'all');

v_artery = sum(v_RMS.*maskArtery, [1 2])/nnz(maskArtery);
v_vein = sum(v_RMS.*maskVein, [1 2])/nnz(maskVein);
Vmax_Arteries = max(v_artery(:));
Vmax_Veins = max(v_vein(:));
Vmin_Arteries = min(v_artery(:));
Vmin_Veins = min(v_vein(:));

%% Construct velocity map

[hue_artery,sat_artery,val_artery,cmap_artery] = createHSVmap(Im,maskArtery,0,0.18); % 0 / 0.18 for orange-yellow range
[hue_vein,sat_vein,val_vein,cmap_vein] = createHSVmap(Im,maskVein,0.68,0.5); %0.5/0.68 for cyan-dark blue range
val = Im.*(~(maskArtery+maskVein))+val_artery.*maskArtery+val_vein.*maskVein;
flowImageRGB =  hsv2rgb(hue_artery+hue_vein, sat_artery+sat_vein, val);
figure(321)
imshow(flowImageRGB)

%% Construct Velocity video 
flowVideoRGB_one_cycle = zeros(size(v_RMS,1),size(v_RMS,2),3,size(v_RMS,3));

for ii = 1:size(v_RMS,3)
    v = mat2gray(squeeze(v_RMS(:,:,ii)));
    [hue_artery,sat_artery,val_artery,cmap_artery] = createHSVmap(v,maskArtery,0,0.18); % 0 / 0.18 for orange-yellow range
    [hue_vein,sat_vein,val_vein,cmap_vein] = createHSVmap(v,maskVein,0.68,0.5); %0.5/0.68 for cyan-dark blue range
    val = v.*(~(maskArtery+maskVein))+val_artery.*maskArtery+val_vein.*maskVein;
    flowVideoRGB_one_cycle(:,:,:,ii) =   hsv2rgb(hue_artery+hue_vein, sat_artery+sat_vein, val);

    
end

v_histo_artery = round(v_RMS.*maskArtery);
v_min = min(v_histo_artery,[],'all');
v_max = max(v_histo_artery,[],'all');

X = linspace(v_min,v_max,v_max-v_min+1);
n = size(X,2);
histo = zeros(size(X,2),size(v_RMS,3));
for t = 1:size(v_RMS,3)
    for x = 1:size(v_RMS,1)
        for y = 1:size(v_RMS,2)
            if( v_histo_artery(x,y,t) ~= 0) 
            i = find(X == v_histo_artery(x,y,t) );
            histo(i,t) = histo(i,t) + 1;
            end


        end
    end
end 



figure(156)
yAx = [v_min v_max];
xAx = [0 n*ToolBox.stride/(1000*ToolBox.fs)];
imagesc(xAx,yAx,histo)
set(gca,'YDir','normal')
colormap("hot")
ylabel('Velocity (mm.s^{-1})')
xlabel('Time (s)')
title("Velocity distribution in arteries")



v_histo_veins = round(v_RMS.*maskArtery);
v_min = min(v_histo_veins,[],'all');
v_max = max(v_histo_veins,[],'all');

X = linspace(v_min,v_max,v_max-v_min+1);
n = size(X,2);
histo = zeros(size(X,2),size(v_RMS,3));
for t = 1:size(v_RMS,3)
    for x = 1:size(v_RMS,1)
        for y = 1:size(v_RMS,2)
            if( v_histo_veins(x,y,t) ~= 0) 
            i = find(X == v_histo_veins(x,y,t) );
            histo(i,t) = histo(i,t) + 1;
            end


        end
    end
end 



figure(157)
yAx = [v_min v_max];
xAx = [0 n*ToolBox.stride/(1000*ToolBox.fs)];
imagesc(xAx,yAx,histo)
set(gca,'YDir','normal')
colormap("hot")
ylabel('Velocity (mm.s^{-1})')
xlabel('Time (s)')
title("Velocity distribution in arteries")
% h = histogram(v_RMS(:,:,ii).*maskArtery);
% X = h.BinCounts;
% 
% [x,y] = max(X);
% X(y) = 0;
% X = smoothdata(X);
% figure
% plot(X)


% save video
% avi
w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_flowVideo_one_cycle'))) ;
open(w)
for jj = 1:size(flowVideoRGB_one_cycle,4)
    writeVideo(w,squeeze(flowVideoRGB_one_cycle(:,:,:,jj))) ;
end
close(w);
% mp4
w = VideoWriter(fullfile(ToolBox.PW_path_mp4,strcat(ToolBox.main_foldername,'_flowVideo_one_cycle')),'MPEG-4') ;
open(w)
for jj = 1:size(flowVideoRGB_one_cycle,4)
    writeVideo(w,squeeze(flowVideoRGB_one_cycle(:,:,:,jj))) ;
end
close(w);

try
    % Save colorbar
    colorfig = figure(116);
    colorfig.Units = 'normalized';
    colormap(cmap_artery)
    %hCB = colorbar('north');
    hCB = colorbar('north','Ticks',[0,1],'TickLabels',{string(round(Vmin_Arteries,1)),string(round(Vmax_Arteries,1))});
    set(gca,'Visible',false)
    set(gca,'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca,15,"points") ;


    % Save colorbar
    colorfig = figure(117);
    colorfig.Units = 'normalized';
    colormap(cmap_vein)
    %hCB = colorbar('north');
    hCB = colorbar('north','Ticks',[0,1],'TickLabels',{string(round(Vmin_Veins,1)),string(round(Vmax_Veins,1))});
    set(gca,'Visible',false)
    set(gca,'LineWidth', 3);
    hCB.Position = [0.10 0.3 0.81 0.35];
    colorfig.Position(4) = 0.1000;
    fontsize(gca,15,"points") ;


    print('-f116','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_colorbar_velocity_arteries.png'))) ;
    print('-f117','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_colorbar_velocity_veins.png'))) ;

catch
    disp('fail saving colorbars')
end
print('-f156','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_velocity_distribution_arteries.png'))) ;
print('-f157','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_velocity_distribution_veins.png'))) ;
imwrite(flowImageRGB, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_flow_image.png')));

end
