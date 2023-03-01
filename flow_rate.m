function [plot_printed] = flow_rate(maskArtery, maskArteryInPlane, v_RMS)
%SECTION_PLOT Summary of this function goes here
%   Detailed explanation goes here

type_of_selection = "Automatic";
nb_sides = 120;

tmp = interp2(squeeze(v_RMS(:,:,1)),1);
v_RMS_interp = zeros(size(tmp,1),size(tmp,2),size(v_RMS,3));
for jj=1:size(v_RMS,3)
    v_RMS_interp(:,:,jj) = interp2(v_RMS(:,:,jj),1);
end

%FIXME rename instead of copy
v_RMS = v_RMS_interp;
clear('v_RMS_interp');

maskArtery = logical(mat2gray(imread( ...
    fullfile("C:\Users\Rakushka\Downloads",...
    strcat("210218_GUN0323_OD_ONH_0_DopplerRMS",'_maskArteryHD_2.png')),'png')...
    ));

figure(100)
imshow(double(maskArtery))
points_selected = 0;
points_x = [];
points_y = [];

blurred_mask = imgaussfilt(double(mean(v_RMS,3).*maskArtery),round(size(maskArtery,1)/4));

[~,x_center] = findpeaks(sum(blurred_mask,1));
[~,y_center] = findpeaks(sum(blurred_mask,2));

polygon = nsidedpoly(nb_sides, 'Center', [x_center, y_center], 'Radius', 150);
points_x = polygon.Vertices(:,1);
points_x(end + 1) = points_x(1);
points_y = polygon.Vertices(:,2);
points_y(end + 1) = points_y(1);
for ii = 1:nb_sides
    l = line([points_x(ii), points_x(ii + 1)], [points_y(ii), points_y(ii + 1)]);
    l.Color = 'red';
    l.LineWidth = 3;
end

%Vertices, Edges

[cx, cy, ~] = improfile(maskArtery, points_x, points_y);

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


v_RMS_masked = v_RMS .* maskArtery;

new_plot = zeros(jj, size(v_RMS_masked, 3));
for ii=1:jj
    new_plot(ii, :) = squeeze(v_RMS_masked(round(cy(ii)), round(cx(ii)), :));
end

param_peak = std(new_plot(:));

plot_values = squeeze(mean(new_plot, 2));
[pks,locs,width,~] = findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);

filtered_pks = pks;
filtered_locs = locs;
filtered_width = width;
for ii=1:size(width)
    if width(ii)<3 % 3 pixels threshold to delete spurious peaks
        filtered_width(ii) = [];
        filtered_locs(ii) = [];
        filtered_pks(ii) = [];
    end
end
pks = filtered_pks;
locs = filtered_locs;
width = filtered_width;


% Average the blood flow calculation over a circle before dividing by the section
avg_blood_flow = zeros(length(width),1);
avg_blood_rate = zeros(length(width),1);
for ii=1:size(width)
    %     abs_x = point_x(locs(ii));
    %     abs_y = point_y(locs(ii));
    center = [cx(locs(ii)) cy(locs(ii))];
    radii = 10; % nb pixel
    [xx,yy] = meshgrid(1:size(blurred_mask,1),1:size(blurred_mask,2));
    circular_mask = false(size(blurred_mask,1),size(blurred_mask,2));
    circular_mask = circular_mask | hypot(xx - center(1), yy - center(2)) <= radii;
    mask_slices = circular_mask.*maskArtery;
    %FIXME : show green slices of selected area
    %FIXME mean(v_RMS,3)
    avg_blood_flow(ii) = sum(nnz(mean(v_RMS,3).*mask_slices))/numel(nnz(mean(v_RMS,3).*mask_slices));
    avg_blood_rate(ii) = avg_blood_flow(ii)*pi*(width(ii)*0.0102/2)^2; % mm^3/s % 0.0102mm = size 1 pixel
end

figure(100)
for ii=1:size(locs)
    text(cx(locs(ii)),cy(locs(ii)),string(round(avg_blood_rate(ii),3)), "FontWeight", "bold", "Color", "white", "BackgroundColor", "blue");
end
title("Blood rate (mm^3/s)")

figure(101)
plot(plot_values);
findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
% text(locs,pks,num2str((1:numel(pks))'))
text(locs,pks,string(round(avg_blood_rate,3)))
title("Peaks of luminosity")
pbaspect([1.618 1 1]);

figure(102)
h = pcolor(new_plot);
set(h, 'EdgeColor', 'none');
yticks(locs);
yticklabels(1:size(locs));
title("Luminosity evolution");
pbaspect([1.618 1 1]);
axis tight;
plot_printed = true;

for ii=1:size(locs)
    figure(ii)
    plot(new_plot(locs(ii),:));
end

%%


subImgHW = round(min([size(maskArtery,1) size(maskArtery,2)]*0.08))

%coords img du peak
% cx(locs(1))
% cy(locs(1))
img = squeeze(mean(v_RMS,3)) .* maskArtery;

for ii = 1:size(locs)
    %FIXME bords d IMG
    xRange = round((-subImgHW/2:1:subImgHW/2) + cx(locs(ii)));
    yRange = round((-subImgHW/2:1:subImgHW/2) + cy(locs(ii)));
    subImg = img(yRange,xRange);

    figure(1000*ii)
    imagesc(subImg);
    %make disk mask
    %FIXME img anamorphique
    subImg = cropCircle(subImg);
    figure(1001*ii)
    imagesc(subImg);


    theta = linspace(0,180,181);
    projx = zeros(size(subImg,1),length(theta));
    projy = zeros(size(subImg,2),length(theta));
    for tt = 1:length(theta)
        tmpImg = imrotate(subImg,theta(tt),'bilinear','crop');
        projx(:,tt) = squeeze(sum(tmpImg,1));
        projy(:,tt) = squeeze(sum(tmpImg,2));
    end

    max_projx = max(max(projx)); 
    [x_max,~] = find(projx == max_projx) %x_max angle de rotation pour une coupe normale au vaisseau
    subImg = imrotate(subImg,x_max,'bilinear','crop');
    section_cut = projx(x_max,:);
    [~,~,width(ii),~] = findpeaks(section_cut,size(projx,1),"MinPeakProminence",std(section_cut));
    width(ii)
end

end


