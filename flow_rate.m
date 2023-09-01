function [] = flow_rate(maskArtery, maskVein, maskCRA,maskSectionArtery, v_RMS,dataM0, ToolBox, k,path)

PW_params = Parameters(path);

%maskArtery = imdilate(maskArtery,strel('disk',5));
%FIXME function velocity map

v_RMS_AVG = mean(v_RMS,3);

%%Find the locations of the sections

maskSectionArtery = maskSectionArtery.*maskArtery;

figure(111)
imagesc(maskSectionArtery.*v_RMS_AVG);

maskSectionArtery = bwlabel(maskSectionArtery);

figure(222)
imagesc(maskSectionArtery)

nb_section = max(maskSectionArtery,[],'all');
masksSections = zeros(size(maskArtery,1),size(maskArtery,2),nb_section);

for ii = 1:nb_section
    masksSections(:,:,ii) = (maskSectionArtery == ii);
%     skel = bwskel(logical(masksSections(:,:,ii)));
%     [row,col] = find(skel);
%     x_center = round(mean(row));
%     y_center = round(mean(col));
%     for i = 1:size(masksSections(:,:,ii),1)
%         for j = 1:size(masksSections(:,:,ii),1)
%             if sqrt((i-x_center)^2+(j-y_center)^2)> (0.8*ecart)*size(v_RMS,1)
%                 masksSections(i,j,ii) = 0;
%             end
%         end
%     end
end 

SubImg_locs = zeros(nb_section,2);
SubImg_width = zeros(nb_section,1);
for ii = 1:nb_section
%     skel = bwskel(logical(masksSections(:,:,ii)));
%     [row,col] = find(skel);
    [row,col] = find(masksSections(:,:,ii));
    SubImg_locs(ii,1) = round(mean(row));
    SubImg_locs(ii,2) = round(mean(col));
    SubImg_width(ii) = 0.01*size(maskArtery,1);
    %width(ii) = max(abs(row(1)-row(end)),abs(col(1)-col(end)));
end

%mask_artery = imdilate(maskArtery,strel('disk',5));
mask_artery = maskArtery;
mask_vein = maskVein;
[avg_blood_rate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery,total_blood_volume_rate_artery] = cross_section_analysis_new(SubImg_locs, SubImg_width, mask_artery, v_RMS, PW_params.flowRate_sliceHalfThickness, k,ToolBox,path);
%[avg_blood_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein,total_blood_volume_rate_vein] = cross_section_analysis_new(SubImg_locs, SubImg_width, mask_vein, v_RMS, PW_params.flowRate_sliceHalfThickness, k,ToolBox,path);


maskRGB = ones(size(maskArtery,1), size(maskArtery,2), 3);


%maskRGB = maskRGB .*maskArtery;

% maskRGB(:,:,3) = maskRGB(:,:,3) - maskSectionArtery;
% maskRGB(:,:,2) = maskRGB(:,:,2) - maskSectionArtery;
mean_M0 = mean(dataM0,3);
maskRGB(:,:,3) = mat2gray(mean_M0).*~cross_section_mask_artery;
maskRGB(:,:,2) = mat2gray(mean_M0).*~cross_section_mask_artery;
maskRGB(:,:,1) = mat2gray(mean_M0).*~cross_section_mask_artery+maskRGB(:,:,1).*cross_section_mask_artery;


figure(111111)
imshow(maskRGB)
% max_dataM0 = max(dataM0,[],'all');
% maskRGB_video  = ones(size(maskArtery,1), size(maskArtery,2), 3,size(dataM0,3));
% maskRGB_video(:,:,3,:) = (dataM0/max_dataM0).*~cross_section_mask_artery;
% maskRGB_video(:,:,2,:) = (dataM0/max_dataM0).*~cross_section_mask_artery;
% maskRGB_video(:,:,1,:) = (dataM0/max_dataM0).*~cross_section_mask_artery+maskRGB_video(:,:,1,:).*cross_section_mask_artery;

ratio_etiquette = 1.2;

figure(111222)
imagesc(mean(v_RMS,3).*cross_section_mask_artery)
colormap("gray")
Total_blood_flow_rate = sum(avg_blood_rate_artery);
disp('Total_blood_flow_rate')
disp(Total_blood_flow_rate)

title('Artery Sections' );
axis image
axis off


figure(120)
imshow(maskRGB);
title(['Total blood volume rate : ' num2str(total_blood_volume_rate_artery(1)) ' µL/min (arteries) - ']);
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
Total_blood_volume_rate = getframe(ax,rect);
[M,N,~] = size(Total_blood_volume_rate.cdata);
volume_rate_video_artery = zeros(M,N,3,size(v_RMS,3));

%volume_rate_video_artery = zeros(size(v_RMS,1),size(v_RMS,2),3,size(v_RMS,3));

for tt = 1:size(v_RMS,3)
    imshow(maskRGB);
    x_center = ToolBox.x_barycentre;
    y_center = ToolBox.y_barycentre;

    for ii=1:size(avg_blood_rate_artery,1)
        new_x = x_center + ratio_etiquette*(SubImg_locs(ii,2)-x_center);
        new_y = y_center + ratio_etiquette*(SubImg_locs(ii,1)-y_center);
        text(new_x, new_y, string(round(avg_blood_rate_artery(ii,tt),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
    end

    title(['Total blood volume rate : ' num2str(round(total_blood_volume_rate_artery(tt))) ' µL/min (arteries) - ']);
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    marg = 30;
    pos = ax.Position;
    rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    Total_blood_volume_rate = getframe(ax,rect);
    volume_rate_video_artery(:,:,:,tt) = Total_blood_volume_rate.cdata;


%     figure(115)
%     imshow(maskRGB);
% 
% 
%     for ii=1:size(avg_blood_rate_artery)
%         new_x = x_center + ratio_etiquette*(SubImg_locs(ii,2)-x_center);
%         new_y = y_center + ratio_etiquette*(SubImg_locs(ii,1)-y_center);
%         text(new_x, new_y, string(round(avg_blood_velocity_artery(ii),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
%     end
% 
%     title('Velocity map in arterial vessels (mm.s^-1)');
%     % png
%     % print('-f120','-dpng',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.png'))) ;
%     drawnow
%     ax = gca;
%     ax.Units = 'pixels';
%     marg = 30;
%     pos = ax.Position;
%     rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
%     F_total_blood_flow = getframe(ax,rect);
end

figure(121)
imshow(maskRGB);
x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;

for ii=1:size(avg_blood_rate_artery,1)
    new_x = x_center + ratio_etiquette*(SubImg_locs(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(SubImg_locs(ii,1)-y_center);
    text(new_x, new_y, string(round(mean(avg_blood_rate_artery(ii,:),2))), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end

title(['Total blood volume rate : ' num2str(round(mean(total_blood_volume_rate_artery))) ' µL/min (arteries) - ']);
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_Total_blood_volume_rate = getframe(ax,rect);
volume_rate_video_artery(:,:,:,tt) = F_Total_blood_volume_rate.cdata;


plot_vol_rate_artery = figure(125);
plot_vol_rate_artery.Position(3:4) = [1200 550];
plot(total_blood_volume_rate_artery)
set(gca,'YDir','normal')
ylabel('Total blood volume rate (µL/min (arteries))')
xlabel('Time (s)')
title("Total blood volume rate in arteries")
set(gca,'PlotBoxAspectRatio',[2.5 1 1])



print('-f111222','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Artery_Sections.png'))) ;
print('-f125','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Plot_blood_volume_rate_in_arteries.png'))) ;

%print('-f115','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png'))) ;
%print('-f120','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png'))) ;

%imwrite(F_total_blood_flow.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png')));
imwrite(F_Total_blood_volume_rate.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png')));


w = VideoWriter(fullfile(ToolBox.PW_path_avi,strcat(ToolBox.main_foldername,'_ volume_rate_video_artery.avi')));
tmp = mat2gray( volume_rate_video_artery);
open(w)
for j = 1:size( volume_rate_video_artery,4)
    writeVideo(w,tmp(:,:,:,j)) ;  
end
close(w);
%close all 
end