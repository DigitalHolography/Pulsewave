function [flowVideoRGB] = flow_rate(maskArtery, maskVein, maskCRA, v_RMS, one_cycle_dir, filename, k,path)
%SECTION_PLOT Summary of this function goes here
%   Detailed explanation goes here
% k = interpolation 2^k-1 (for size pixel for section calculation)
PW_params = Parameters(path);

one_cycle_dir_png = fullfile(one_cycle_dir, 'png');
one_cycle_dir_eps = fullfile(one_cycle_dir, 'eps');
one_cycle_dir_txt = fullfile(one_cycle_dir, 'txt');
one_cycle_dir_avi = fullfile(one_cycle_dir, 'avi');


%% Define a section (circle)

%FIXME : radius_ratio as an entry param
%[~,radius_ratio] = getPulsewaveParamFromTxt(path,'Radius ratio :');

radius = round(PW_params.radius_ratio* size(v_RMS,1));
%FIXME : anamorphic image
blurred_mask = imgaussfilt(double(mean(v_RMS,3).*double(maskCRA)),round(size(maskCRA,1)/4),'Padding',PW_params.flowRate_gaussFiltPadding);
[~,x_center] = findpeaks(sum(blurred_mask,1));
[~,y_center] = findpeaks(sum(blurred_mask,2));

polygon = nsidedpoly(PW_params.nbSides, 'Center', [x_center, y_center], 'Radius', radius);
points_x = polygon.Vertices(:,1);
points_x(end + 1) = points_x(1);
points_y = polygon.Vertices(:,2);
points_y(end + 1) = points_y(1);
figure(121)
for ii = 1:PW_params.nbSides
    l   = line([points_x(ii), points_x(ii + 1)], [points_y(ii), points_y(ii + 1)]);
    l.Color = 'red';
    l.LineWidth = 2;
end

%Vertices, Edges
[cx, cy, ~] = improfile(maskArtery+maskVein, points_x, points_y);

jj = 0;
% Delete all points which are not in the maskArtery
for ii=1:size(cx, 1)
    ry = round(cy(ii));
    rx = round(cx(ii));
    if (ry > 0 && ry <= size(v_RMS, 1) && rx > 0 && rx <= size(v_RMS, 2))
        jj = jj + 1;
        cy(jj) = ry;
        cx(jj) = rx;
    end
end
if (jj == 0) %If no points, no analysis.
    plot_printed = false;
    return;
end

%% Peak detection for maskVein
figure(154)
imshow(double(maskVein))

[pks_Vein, locs_Vein, width_Vein] = find_cross_section(maskVein, v_RMS, cx, cy, jj);

%% Peak detection for maskArtery
figure(100)
imshow(double(maskArtery))

[pks_Artery, locs_Artery, width_Artery] = find_cross_section(maskArtery, v_RMS, cx, cy, jj);

%% Images HSV Artery-Vein Retina
eta_sat = 0.1;
eta_val = 0.1;
img_v_artery = squeeze(mean(v_RMS,3)) .* maskArtery;


hue = mat2gray(squeeze(mean(v_RMS,3)))*0.18 .* maskArtery; % 0.18 for orange-yellow range in HSV
hue = hue + abs(-mat2gray(squeeze(mean(v_RMS,3)))*0.18 .* maskVein + 0.68 .* maskVein); % x0.18 + 0.5 for cyan-dark blue range in HSV
sat = 1.0 * double(or(maskArtery,maskVein)).* squeeze(mean(v_RMS,3));
val = squeeze(mean(v_RMS,3));
val = mat2gray(val);
tolVal = [0.02, 0.98];
lowhigh = stretchlim(val, tolVal); % adjust contrast a bit
val = mat2gray(imadjust(val, stretchlim(val, tolVal)));
flowMapRGB =  hsv2rgb(hue, sat, val);

%% Video HSV Artery-Vein Retina. Velocity & blood volume rate video
flowVideoRGB = zeros(size(v_RMS,1),size(v_RMS,2),3,size(v_RMS,3));
v_RMS_n = mat2gray(v_RMS);
img_backg = squeeze(mean(v_RMS,3));
img_backg = mat2gray(img_backg);
img_backg = imadjust(img_backg, stretchlim(img_backg, tolVal));
v_artery = sum(v_RMS.*maskArtery, [1 2])/nnz(maskArtery);
v_vein = sum(v_RMS.*maskVein, [1 2])/nnz(maskVein);
Vmax_Arteries = max(v_artery(:));
Vmax_Veins = max(v_vein(:));
Vmin_Arteries = min(v_artery(:));
Vmin_Veins = min(v_vein(:));


adjustedVideo = mat2gray(v_RMS_n);
avgAdjustedVideo = squeeze(mean(adjustedVideo,3));
tolVal = [0.1, 0.99]; 
lowhighVal = stretchlim(avgAdjustedVideo, tolVal); % adjust video contrast a bit 
for ii = 1:size(v_RMS_n,3)
    v = squeeze(v_RMS(:,:,ii));
    hue = mat2gray(v) * 0.18 .* maskArtery; % 0.18 for orange-yellow range in HSV
    hue = hue + abs(-mat2gray(v)*0.18 .* maskVein + 0.68 .* maskVein); % x0.18 + 0.5 for cyan-dark blue range in HSV
    sat = (1.0 - abs(0.5 * mat2gray(v))) .* double(or(maskArtery,maskVein)) ;
    val = mat2gray(v);
    tolVal = [0.02, 0.98];
    lowhigh = stretchlim(val, tolVal); % adjust contrast a bit
    val = mat2gray(imadjust(val, stretchlim(val, tolVal)));
    flowVideoRGB(:,:,:,ii) = hsv2rgb(hue, sat, val);
end
% save video
w = VideoWriter(fullfile(one_cycle_dir_avi,strcat(filename,'_flowVideo'))) ;
open(w)
% flow_video_to_save = mat2gray(flowVideoRGB);% 2nd normalization
for jj = 1:size(flowVideoRGB,4)
    writeVideo(w,squeeze(flowVideoRGB(:,:,:,jj))) ;
end
close(w);

v = squeeze(mean(v_RMS, 3));
hue = mat2gray(v) * 0.18 .* maskArtery; % 0.18 for orange-yellow range in HSV
hue = hue + abs(-mat2gray(v)*0.18 .* maskVein + 0.68 .* maskVein); % x0.18 + 0.5 for cyan-dark blue range in HSV
% sat = 1.0 * double(or(maskArtery,maskVein)) .* mat2gray(v);
sat = (1.0 - abs(0.5 * mat2gray(v))) .* double(or(maskArtery,maskVein)) ;
val = mat2gray(v);
tolVal = [0.02, 0.98];
lowhigh = stretchlim(val, tolVal); % adjust contrast a bit
val = mat2gray(imadjust(val, stretchlim(val, tolVal)));
flowImageRGB = hsv2rgb(hue, sat, val);

figure(321)
imshow(flowImageRGB);

% Save colorbar flow image
color_bounds = [0.1 0.99];
tmp = nonzeros(mat2gray(v).*double(maskArtery));
sat_min_max = stretchlim(tmp, color_bounds);
min_sat = sat_min_max(1);
max_sat = sat_min_max(2);

list = linspace(min_sat, max_sat ,256);
hue = 0.18*list;
sat = 1.0 - 0.5*list;
val = zeros(1,256);
idx = 1;
for ii = list
    if ii < lowhigh(1)
        val(idx) = 0;
    elseif ii > lowhigh(2)
        val(idx) = 1;
    else
        val(idx) = (ii - lowhigh(1)) /(lowhigh(2)-lowhigh(1));
    end
    idx = idx + 1;
end
cmap_arteries = squeeze(hsv2rgb(hue,sat,val));

hue = flip(list)*0.18 + 0.5;
cmap_veins = squeeze(hsv2rgb(hue,sat,val));

colorfig = figure(3210);
%colorbar arteries
colormap(cmap_arteries);
colorfig.Units = 'normalized';
hCB = colorbar('northoutside','Ticks',[0,1],'TickLabels',{string(round(min(Vmin_Arteries,Vmin_Veins),1)), string(round(max(Vmax_Arteries,Vmax_Veins),1))});
set(gca,'Visible',false)
set(gca,'LineWidth', 3);
hCB.Position = [0.1 0.3 0.81 0.35];
colorfig.Position(4) = 0.100;
fontsize(gca,15,"points");
colorTitleHandle = get(hCB,'Title');
titleString = 'Arterial velocity (mm/s)';
set(colorTitleHandle ,'String',titleString);

colorfig2 = figure(3209);
%colorbar veins
colorfig2.Units = 'normalized';
colormap(cmap_veins)
hCB = colorbar('northoutside','Ticks',[0,1],'TickLabels',{string(round(min(Vmin_Arteries,Vmin_Veins),1)), string(round(max(Vmax_Arteries,Vmax_Veins),1))});
set(gca,'Visible',false)
set(gca,'LineWidth', 3);
hCB.Position = [0.1 0.3 0.81 0.35];
colorfig2.Position(4) = 0.100;
fontsize(gca,15,"points");
colorTitleHandle = get(hCB,'Title');
titleString = 'Venous velocity (mm/s)';
set(colorTitleHandle ,'String',titleString);


% Average the blood flow calculation over a rectangle before dividing by the section for blood volume rate
avg_blood_velocity_artery = zeros(length(width_Artery),1);
avg_blood_rate_artery = zeros(length(width_Artery),1);
cross_section_area_artery = zeros(length(width_Artery),1);

avg_blood_velocity_vein = zeros(length(width_Vein),1);
avg_blood_rate_vein = zeros(length(width_Vein),1);
cross_section_area_vein = zeros(length(width_Vein),1);

slice_half_thickness = PW_params.flowRate_sliceHalfThickness; % size of the rectangle area for velocity averaging (in pixel)

%% pour chaque veine jj detectee
[avg_blood_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein] = ...
    cross_section_analysis(locs_Vein, width_Vein, maskVein, cx, cy, v_RMS, slice_half_thickness, k, one_cycle_dir, filename, 'vein',path);
avg_blood_rate_vein_muLmin = avg_blood_rate_vein*60;

%% pour chaque artere ii detectee
[avg_blood_rate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery] = ...
    cross_section_analysis(locs_Artery, width_Artery, maskArtery, cx, cy, v_RMS, slice_half_thickness, k, one_cycle_dir, filename, 'artery',path);
avg_blood_rate_artery_muLmin = avg_blood_rate_artery*60;
%% Display final blood volume rate image
total_blood_rate_artery = sum(avg_blood_rate_artery(:));
total_cross_section_artery = sum(cross_section_area_artery(:));
total_blood_rate_artery_muLmin = total_blood_rate_artery * 60;
disp(['Total cross section of arteries : ' num2str(total_cross_section_artery) ' mm^2']);
disp(['Total blood volume rate in arteries : ' num2str(total_blood_rate_artery) ' mm^3/s']);
disp(['Total blood volume rate in arteries : ' num2str(total_blood_rate_artery_muLmin) ' µL/min']);

total_blood_rate_vein = sum(avg_blood_rate_vein(:));
total_cross_section_vein = sum(cross_section_area_vein(:));
total_blood_rate_vein_muLmin = total_blood_rate_vein * 60;
disp(['Total cross section of veins : ' num2str(total_cross_section_vein) ' mm^2']);
disp(['Total blood volume rate in vein : ' num2str(total_blood_rate_vein) ' mm^3/s']);
disp(['Total blood volume rate in vein : ' num2str(total_blood_rate_vein_muLmin) ' µL/min']);

flowMapArteryRGB = flowMapRGB .* cross_section_mask_artery;
figure(118)
imshow(flowMapArteryRGB)
for ii=1:size(locs_Artery)
    text(cx(locs_Artery(ii)),cy(locs_Artery(ii))+15,string(round(avg_blood_rate_artery(ii),1)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(round(total_blood_rate_artery,1)) ' mm^3/s. '  num2str(round(total_blood_rate_artery_muLmin,1)) ' µL/min']);

flowMapVeinRGB = flowMapRGB .* cross_section_mask_vein;
figure(119)
imshow(flowMapVeinRGB)
for ii=1:size(locs_Vein)
    text(cx(locs_Vein(ii)),cy(locs_Vein(ii))+15,string(round(avg_blood_rate_vein(ii),1)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(round(total_blood_rate_vein,1)) ' (mm^3/s)' num2str(round(total_blood_rate_vein_muLmin,1)) ' (µL/min)']);


maskRGB = ones(size(maskArtery,1), size(maskArtery,2), 3);
maskRGB = maskRGB .* (maskArtery+maskVein);

maskRGB(:,:,1) = maskRGB(:,:,1) - cross_section_mask_vein;
maskRGB(:,:,2) = maskRGB(:,:,2) - cross_section_mask_vein;

maskRGB(:,:,3) = maskRGB(:,:,3) - cross_section_mask_artery;
maskRGB(:,:,2) = maskRGB(:,:,2) - cross_section_mask_artery;

figure(121)
imshow(maskRGB)
for ii = 1:PW_params.nbSides
    l   = line([points_x(ii), points_x(ii + 1)], [points_y(ii), points_y(ii + 1)]);
    l.Color = "#C8CDCD"; % gray line
    l.LineWidth = 2;
end
for ii=1:size(locs_Vein)
    num = string(ii);
    new_x = x_center + 1.2*(cx(locs_Vein(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Vein(ii))-y_center);
    text(new_x, new_y, strcat('V',num), "FontWeight", "bold","FontSize", 20,   "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(locs_Artery)
    num = string(ii);
    new_x = x_center + 1.2*(cx(locs_Artery(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Artery(ii))-y_center);
    text(new_x, new_y, strcat('A',num), "FontWeight", "bold","FontSize", 20,   "Color", "white", "BackgroundColor", "black");
end
% png
% print('-f121','-dpng',fullfile(one_cycle_dir,strcat(filename,'_MaskTopologyAV.png'))) ;
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_MaskTopology = getframe(ax,rect);



figure(120)
imshow(maskRGB);
for ii=1:size(locs_Vein)
    new_x = x_center + 1.2*(cx(locs_Vein(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Vein(ii))-y_center);
    text(new_x, new_y, string(round(avg_blood_rate_vein_muLmin(ii),1)), "FontWeight", "bold", "FontSize", 20,  "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(locs_Artery)
    new_x = x_center + 1.2*(cx(locs_Artery(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Artery(ii))-y_center);
    text(new_x, new_y, string(round(avg_blood_rate_artery_muLmin(ii),1)), "FontWeight", "bold","FontSize", 20,   "Color", "white", "BackgroundColor", "black");
end

title(['Total blood volume rate : ' num2str(round(total_blood_rate_artery_muLmin,1)) ' µL/min (arteries) - ' num2str(round(total_blood_rate_vein_muLmin,1)) ' µL/min (veins)']);
% png
% print('-f120','-dpng',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.png'))) ;
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_Total_blood_volume_rate = getframe(ax,rect);


figure(115)
imshow(maskRGB);
for ii=1:size(locs_Vein)
    new_x = x_center + 1.2*(cx(locs_Vein(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Vein(ii))-y_center);
    text(new_x, new_y, string(round(avg_blood_velocity_vein(ii),1)), "FontWeight", "bold", "FontSize", 20,  "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(locs_Artery)
    new_x = x_center + 1.2*(cx(locs_Artery(ii))-x_center);
    new_y = y_center + 1.2*(cy(locs_Artery(ii))-y_center);
    text(new_x, new_y, string(round(avg_blood_velocity_artery(ii),1)), "FontWeight", "bold","FontSize", 20,   "Color", "white", "BackgroundColor", "black");
end

title('Velocity map in venous & arterial vessels');
% png
% print('-f120','-dpng',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.png'))) ;
drawnow
ax = gca;
ax.Units = 'pixels';
marg = 30;
pos = ax.Position;
rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
F_total_blood_flow = getframe(ax,rect);


% txt file output with measured pulse wave parameters

fileID = fopen(fullfile(one_cycle_dir_txt,strcat(filename,'_pulseWaveParameters.txt')),'a') ;
fprintf(fileID,[...
    'Value of total arterial blood volume rate (µL/min) :\n%d\n' ...
    'Value of total arterial blood volume rate (mm^3/s) :\n%d\n' ...
    'Value of total venous blood volume rate (µL/min) :\n%d\n' ...
    'Value of total venous blood volume rate (mm^3/s) :\n%d\n' ...
    'Total cross section of veins (mm^2) :\n%d\n' ...
    'Total cross section of arteries (mm^2) :\n%d\n'], ...
    total_blood_rate_artery_muLmin, ...
    total_blood_rate_artery, ...
    total_blood_rate_vein_muLmin, ...
    total_blood_rate_vein, ...
    total_cross_section_artery, ...
    total_cross_section_vein);
fclose(fileID) ;

for ii=1:length(avg_blood_rate_artery)
    fileID = fopen(fullfile(one_cycle_dir_txt,strcat(filename,'_pulseWaveOutputParameters.txt')),'a') ;
    fprintf(fileID,[...
        'Artery n°%d : cross_section (mm^2) : \n %d \n ' ...
        'Artery n°%d : vessel diameter (µm) : \n %d \n ' ...
        'Artery n°%d : average velocity (mm/s) : \n %d \n ' ...
        'Artery n°%d : blood rate (µL/min) : \n %d \n '], ...
        ii, ...
        cross_section_area_artery(ii), ...
        ii, ...
        2*sqrt(cross_section_area_artery(ii)/pi)*1000, ... % calculation of the diameter knowing the disc area
        ii, ...
        avg_blood_velocity_artery(ii), ...
        ii, ...
        avg_blood_rate_artery(ii)*60);% mm^3/s -> µL/min 
    fclose(fileID) ;
end

for ii=1:length(avg_blood_rate_vein)
    fileID = fopen(fullfile(one_cycle_dir_txt,strcat(filename,'_pulseWaveParameters.txt')),'a') ;
    fprintf(fileID,[...
        'Vein n°%d : cross_section (mm^2) : \n %d \n ' ...
        'Vein n°%d : vessel diameter (µm) : \n %d \n ' ...
        'Vein n°%d : average velocity (mm/s) : \n %d \n ' ...
        'Vein n°%d : blood rate (µL/min) : \n %d \n '], ...
        ii, ...
        cross_section_area_vein(ii), ...
        ii, ...
        2*sqrt(cross_section_area_vein(ii)/pi)*1000, ... % calculation of the diameter knowing the disc area
        ii, ...
        avg_blood_velocity_vein(ii), ...
        ii, ...
        avg_blood_rate_vein(ii)*60);% mm^3/s -> µL/min 
    fclose(fileID) ;
end

% png
print('-f3209','-dpng',fullfile(one_cycle_dir_png,strcat(filename,'_blood_flow_venous_colorbar.png')));
print('-f3210','-dpng',fullfile(one_cycle_dir_png,strcat(filename,'_blood_flow_arterial_colorbar.png')));
% eps
print('-f3209','-depsc',fullfile(one_cycle_dir_eps,strcat(filename,'_blood_flow_venous_colorbar.eps')));
print('-f3210','-depsc',fullfile(one_cycle_dir_eps,strcat(filename,'_blood_flow_arterial_colorbar.eps')));

imwrite(flowImageRGB, fullfile(one_cycle_dir_png,strcat(filename,'_flow_image.png')));
imwrite(F_total_blood_flow.cdata, fullfile(one_cycle_dir_png,strcat(filename,'_Total_blood_flow.png')));
imwrite(F_MaskTopology.cdata, fullfile(one_cycle_dir_png,strcat(filename,'_MaskTopologyAV.png')));
imwrite(F_Total_blood_volume_rate.cdata, fullfile(one_cycle_dir_png,strcat(filename,'_Total_blood_volume_rate.png')));
list_fig_close = [118,119,121,154,100];
for ii=1:length(list_fig_close)
    close(list_fig_close(ii));
end

end