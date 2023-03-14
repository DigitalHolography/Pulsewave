function [flowVideoRGB] = flow_rate(maskArtery, maskVein, v_RMS, one_cycle_dir, filename)
%SECTION_PLOT Summary of this function goes here
%   Detailed explanation goes here

type_of_selection = "Automatic";
nb_sides = 120;

%% Define a section (circle)

%FIXME : radius_ratio as an entry param
radius_ratio = round(0.15 * size(v_RMS,1));
%FIXME : anamorphic image
blurred_mask = imgaussfilt(double(mean(v_RMS,3).*double(maskArtery+maskVein)),round(size(maskArtery,1)/4),'Padding',0);
[~,x_center] = findpeaks(sum(blurred_mask,1));
[~,y_center] = findpeaks(sum(blurred_mask,2));

polygon = nsidedpoly(nb_sides, 'Center', [x_center, y_center], 'Radius', radius_ratio);
points_x = polygon.Vertices(:,1);
points_x(end + 1) = points_x(1);
points_y = polygon.Vertices(:,2);
points_y(end + 1) = points_y(1);
figure(121)
for ii = 1:nb_sides
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
flowVideoRGB = zeros(size(v_RMS,1),size(v_RMS,2),size(v_RMS,3),3);
v_RMS_n = mat2gray(v_RMS);
for ii = 1:size(v_RMS_n,3)
    v = squeeze(v_RMS_n(:,:,ii));
    img_v_artery = v .* maskArtery;
    hue = v * 0.18 .* maskArtery; % 0.18 for orange-yellow range in HSV
    hue = hue + abs(-v*0.18 .* maskVein + 0.68 .* maskVein); % x0.18 + 0.5 for cyan-dark blue range in HSV
    sat = 1.0 * double(or(maskArtery,maskVein)) .* v;
    val = v;
    tolVal = [0.02, 0.98];
    lowhigh = stretchlim(val, tolVal); % adjust contrast a bit
    val = imadjust(val, stretchlim(val, tolVal));
    tmp = hsv2rgb(hue, sat, val);
    flowVideoRGB(:,:,ii,1) = tmp(:,:,1);
    flowVideoRGB(:,:,ii,2) = tmp(:,:,2);
    flowVideoRGB(:,:,ii,3) = tmp(:,:,3);
end



% Average the blood flow calculation over a rectangle before dividing by the section for blood volume rate
avg_blood_velocity_artery = zeros(length(width_Artery),1);
avg_blood_rate_artery = zeros(length(width_Artery),1);
cross_section_area_artery = zeros(length(width_Artery),1);

avg_blood_velocity_vein = zeros(length(width_Vein),1);
avg_blood_rate_vein = zeros(length(width_Vein),1);
cross_section_area_vein = zeros(length(width_Vein),1);

slice_half_thickness = 10; % size of the rectangle area for velocity averaging (in pixel)

%% pour chaque veine jj detectee
[avg_blood_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein] = ...
    cross_section_analysis(locs_Vein, width_Vein, maskVein, cx, cy, v_RMS, slice_half_thickness);

%% pour chaque artere ii detectee
[avg_blood_rate_artery, cross_section_area_artery, avg_blood_velocity_artery, cross_section_mask_artery] = ...
    cross_section_analysis(locs_Artery, width_Artery, maskArtery, cx, cy, v_RMS, slice_half_thickness);

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
    text(cx(locs_Artery(ii)),cy(locs_Artery(ii))+15,string(round(avg_blood_rate_artery(ii),3)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(total_blood_rate_artery) ' mm^3/s. '  num2str(total_blood_rate_artery_muLmin) ' µL/min']);

flowMapVeinRGB = flowMapRGB .* cross_section_mask_vein;
figure(119)
imshow(flowMapVeinRGB)
for ii=1:size(locs_Vein)
    text(cx(locs_Vein(ii)),cy(locs_Vein(ii))+15,string(round(avg_blood_rate_vein(ii),3)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(total_blood_rate_vein) ' (mm^3/s)' num2str(total_blood_rate_vein_muLmin) ' (µL/min)']);

% maskRGB_artery = zeros(size(maskArtery,1), size(maskArtery,2), 3);
% maskRGB_artery(:,:,:) = maskArtery;
% maskRGB_artery(:,:,2:3) = ~cross_section_mask_artery;
% maskRGB_artery = maskRGB_artery .* maskArtery;
% 
% maskRGB_vein(:,:,3) = cross_section_mask_vein;
% maskRGB_vein = zeros(size(maskVein,1), size(maskVein,2), 3);

maskRGB = ones(size(maskArtery,1), size(maskArtery,2), 3);
maskRGB = maskRGB .* (maskArtery+maskVein);

maskRGB(:,:,1) = maskRGB(:,:,1) - cross_section_mask_vein;
maskRGB(:,:,2) = maskRGB(:,:,2) - cross_section_mask_vein;

maskRGB(:,:,3) = maskRGB(:,:,3) - cross_section_mask_artery;
maskRGB(:,:,2) = maskRGB(:,:,2) - cross_section_mask_artery;

figure(121)
imshow(maskRGB)
for ii = 1:nb_sides
    l   = line([points_x(ii), points_x(ii + 1)], [points_y(ii), points_y(ii + 1)]);
    l.Color = "#C8CDCD"; % gray line
    l.LineWidth = 2;
end
for ii=1:size(locs_Vein)
    num = string(ii);
    text(cx(locs_Vein(ii)),cy(locs_Vein(ii))+15, strcat('V',num), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(locs_Artery)
    num = string(ii);
    text(cx(locs_Artery(ii)),cy(locs_Artery(ii))+15, strcat('A',num), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
% png
print('-f121','-dpng',fullfile(one_cycle_dir,strcat(filename,'_MaskTopologyAV.png'))) ;
% eps
print('-f121','-depsc',fullfile(one_cycle_dir,strcat(filename,'_MaskTopologyAV.eps'))) ;


figure(120)
imshow(maskRGB)
for ii=1:size(locs_Vein)
    text(cx(locs_Vein(ii)),cy(locs_Vein(ii))+15,string(round(avg_blood_rate_vein(ii),3)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
for ii=1:size(locs_Artery)
    text(cx(locs_Artery(ii)),cy(locs_Artery(ii))+15,string(round(avg_blood_rate_artery(ii),3)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "black");
end
title(['Total blood volume rate : ' num2str(total_blood_rate_artery_muLmin) ' µL/min (arteries) - ' num2str(total_blood_rate_vein_muLmin) ' µL/min (veins)']);
% png
print('-f120','-dpng',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.png'))) ;
% eps
print('-f120','-depsc',fullfile(one_cycle_dir,strcat(filename,'_Total_blood_volume_rate.eps'))) ;


% txt file output with measured pulse wave parameters
fileID = fopen(fullfile(one_cycle_dir,strcat(filename,'_pulseWaveParameters.txt')),'a') ;
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

[avg_blood_rate_vein, cross_section_area_vein, avg_blood_velocity_vein, cross_section_mask_vein] = ...
    cross_section_analysis(locs_Vein, width_Vein, maskVein, cx, cy, v_RMS, slice_half_thickness);

for ii=1:length(avg_blood_rate_artery)
    fileID = fopen(fullfile(one_cycle_dir,strcat(filename,'_pulseWaveParameters.txt')),'a') ;
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
    fileID = fopen(fullfile(one_cycle_dir,strcat(filename,'_pulseWaveParameters.txt')),'a') ;
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

end


