function [plot_printed] = section_plot(hologram_loaded, picture, hologram, log_label,type_of_selection,nb_sides)
%SECTION_PLOT Summary of this function goes here
%   Detailed explanation goes here
if (hologram_loaded == false)
    plot_printed = false;
    return;
end




<<<<<<< HEAD
figure(100)
=======
figure(1)
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0
imshow(picture)
log_label.Value = "Please choose a point to analyse";
points_selected = 0;
points_x = [];
points_y = [];


if strcmp(type_of_selection,'Manual') 
    log_label.Value = "Waiting for input (right click to stop)";
    [x_click, y_click, button] = ginput(1);
    while (button ~= 3)
        points_x(end + 1) = round(x_click);
        points_y(end + 1) = round(y_click);
        if (button == 32)
            points_x(end) = points_x(1);
            points_y(end) = points_y(1);
        end
        if (points_selected ~= 0)
            l = line([points_x(end - 1), points_x(end)], [points_y(end - 1), points_y(end)]);
            l.Color = 'red';
            l.LineWidth = 3;
        end
        points_selected = points_selected + 1;
        [x_click, y_click, button] = ginput(1);
    end
else
    [x_click_center, y_click_center] = ginput(1);
    [x_click_radius, y_click_radius] = ginput(1);
    radius = [round(x_click_center), round(y_click_center);round(x_click_radius),round(y_click_radius)];
    polygon = nsidedpoly(nb_sides, 'Center', [x_click_center, y_click_center], 'Radius', sqrt((radius(1,1)-radius(2,1))^2+(radius(1,2)-radius(2,2))^2));
    points_x = polygon.Vertices(:,1);
    points_x(end + 1) = points_x(1);
    points_y = polygon.Vertices(:,2);
    points_y(end + 1) = points_y(1);
    for i = 1:nb_sides
        l = line([points_x(i), points_x(i + 1)], [points_y(i), points_y(i + 1)]);
        l.Color = 'red';
        l.LineWidth = 3;
    end
end
log_label.Value = "Calculations...";
[cx, cy, ~] = improfile(picture, points_x, points_y);
reference = mat2gray(hologram);
j = 0;
% Delete all points which are not in the picture
for i=1:size(cx, 1)
    ry = round(cy(i));
    rx = round(cx(i));
    if (ry > 0 && ry <= size(reference, 1) && rx > 0 && rx <= size(reference, 2))
        j = j + 1;
        cy(j) = ry;
        cx(j) = rx;
    end
end
if (j == 0) %If no points, no analysis.
    log_label.Value = "No points to analyse !";
    plot_printed = false;
    return;
end
new_picture = zeros(j, size(hologram, 3));
for i=1:j
    new_picture(i, :) = squeeze(reference(round(cy(i)), round(cx(i)), :));
end

param_peak = std(new_picture(:));

log_label.Value = "Printing plot...";
plot_values = squeeze(mean(new_picture, 2));
[pks,locs] = findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
for i=1:size(locs)
    text(cx(locs(i)),cy(locs(i)),num2str(i), "FontWeight", "bold", "Color", "white", "BackgroundColor", "blue");
end

<<<<<<< HEAD
figure(101)
=======
figure(2)
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0
plot(plot_values);
findpeaks(plot_values,size(plot_values, 2),'MinPeakProminence',param_peak);
text(locs,pks,num2str((1:numel(pks))'))
title("Peaks of luminosity")
pbaspect([1.618 1 1]);

<<<<<<< HEAD
figure(102)
=======
figure(3)
>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0
h = pcolor(new_picture);
set(h, 'EdgeColor', 'none');
yticks(locs);
yticklabels(1:size(locs));
title("Luminosity evolution");
pbaspect([1.618 1 1]);
axis tight;
log_label.Value = strcat("Plot printed !");
plot_printed = true;
<<<<<<< HEAD

for i=1:size(locs)
    figure(i)
    plot(new_picture(locs(i),:));
end

end


=======
end

>>>>>>> 03c83fcfd740c99fc55cc6bae0ee1a3c04bd84e0
