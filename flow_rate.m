function [] = flow_rate(maskArtery, maskVein, maskCRA,maskSectionArtery, v_RMS, ToolBox, k,path)

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

locs = zeros(nb_section,2);
width = zeros(nb_section,1);
for ii = 1:nb_section
%     skel = bwskel(logical(masksSections(:,:,ii)));
%     [row,col] = find(skel);
    [row,col] = find(masksSections(:,:,ii));
    locs(ii,1) = round(mean(row));
    locs(ii,2) = round(mean(col));
    width(ii) = 0.01*size(maskArtery,1);
    %width(ii) = max(abs(row(1)-row(end)),abs(col(1)-col(end)));
end

%mask_artery = imdilate(maskArtery,strel('disk',5));
mask_artery = maskArtery;
[avg_blood_rate, cross_section_area, avg_blood_velocity, cross_section_mask] = cross_section_analysis_new(locs, width, mask_artery, v_RMS, PW_params.flowRate_sliceHalfThickness, k,ToolBox,path);

%VelocityInArterySections = zeros(size(maskArtery,1),size(maskArtery,2),nb_section);

% for ii = 1:nb_section
%     VelocityInArterySections(:,:,ii) = masksSections(:,:,ii).*v_RMS_AVG;
% end 
% 
% 
% nb_mini_pix = 3*(2*ecart)* size(v_RMS,1);
% 
% avg_blood_rate = zeros(nb_section,1);
% avg_blood_velocity = zeros(nb_section,1);
% size_section =  zeros(nb_section,1);
% 
% for ii = 1:nb_section
%     if nnz(masksSections(:,:,ii)) > nb_mini_pix
%         [avg_blood_rate(ii), avg_blood_velocity(ii), size_section(ii),masksSections(:,:,ii)] = SectionAnalysis(v_RMS_AVG,VelocityInArterySections(:,:,ii),ii);
% 
%        
%     end 
% end

% maskSectionArtery = zeros(size(maskArtery,1),size(maskArtery,2));
% for ii = 1:nb_section
%     maskSectionArtery = maskSectionArtery + masksSections(:,:,ii);
% end
% maskSectionArtery = imdilate(maskSectionArtery,strel('disk',4)).*maskArtery;

maskRGB = ones(size(maskArtery,1), size(maskArtery,2), 3);

maskRGB = maskRGB .*maskArtery;
% maskRGB(:,:,3) = maskRGB(:,:,3) - maskSectionArtery;
% maskRGB(:,:,2) = maskRGB(:,:,2) - maskSectionArtery;

maskRGB(:,:,3) = maskRGB(:,:,3) - cross_section_mask;
maskRGB(:,:,2) = maskRGB(:,:,2) - cross_section_mask;


figure(111111)
imshow(maskRGB)

ratio_etiquette = 1.2;

figure(111222)
imagesc(mean(v_RMS,3).*cross_section_mask)
colormap("gray")
Total_blood_flow_rate = sum(avg_blood_rate);
disp('Total_blood_flow_rate')
disp(Total_blood_flow_rate)

title('Artery Sections' );
axis image
axis off

figure(120)
imshow(maskRGB);
x_center = ToolBox.x_barycentre;
y_center = ToolBox.y_barycentre;

for ii=1:size(avg_blood_rate)
    new_x = x_center + ratio_etiquette*(locs(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(locs(ii,1)-y_center);
    text(new_x, new_y, string(round(avg_blood_rate(ii),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end

title(['Total blood volume rate : ' num2str(Total_blood_flow_rate) ' µL/min (arteries) - ']);
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_Total_blood_volume_rate = getframe(ax,rect);


figure(115)
imshow(maskRGB);


for ii=1:size(avg_blood_rate)
    new_x = x_center + ratio_etiquette*(locs(ii,2)-x_center);
    new_y = y_center + ratio_etiquette*(locs(ii,1)-y_center);
    text(new_x, new_y, string(round(avg_blood_velocity(ii),1)), "FontWeight", "bold","FontSize", 15,   "Color", "white", "BackgroundColor", "black");
end

title('Velocity map in arterial vessels (mm.s^-1)');
% png
% print('-f120','-dpng',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.png'))) ;
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_total_blood_flow = getframe(ax,rect);





print('-f111222','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Artery_Sections.png'))) ;

%print('-f115','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png'))) ;
%print('-f120','-dpng',fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png'))) ;

imwrite(F_total_blood_flow.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_AVG_Velocity_in_arteriy_sections.png')));
imwrite(F_Total_blood_volume_rate.cdata, fullfile(ToolBox.PW_path_png,strcat(ToolBox.main_foldername,'_Total_blood_volume_rate_in_arteries.png')));

%close all 
end