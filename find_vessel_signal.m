function signal = find_vessel_signal(mask, v_RMS, circle_in_img_x, circle_in_img_y, jj)

% find_cross_section detects each peak in a circular patern 
%   pks(ii) : list of peak values of v_RMS
%   locs(ii) : location of peaks in v_RMS
%   width(ii) : cross_section width of segmented vessels
%   circle_in_img_x(ii) : x coords of the ii point in the circular section
%   circle_in_img_y(ii) : y coords of the ii point in the circular section
%   jj : list of points in the image %%(FIXME)

v_RMS_masked = v_RMS .* mask;

circle_path_signal = zeros(jj, size(v_RMS_masked, 3));
circle_path_mask = zeros(jj, 1);
for ii=1:jj
    circle_path_signal(ii, :) = squeeze(v_RMS_masked(round(circle_in_img_y(ii)), round(circle_in_img_x(ii)), :));
    circle_path_mask(ii) = squeeze(mask(round(circle_in_img_y(ii)), round(circle_in_img_x(ii))));
end

% STD along time for any circular point in the image
param_peak = std(circle_path_signal(:));

%FIXME : plot_values another name
plot_values = squeeze(mean(circle_path_signal, 2));
[pks,locs,width,~] = findpeaks(plot_values,1:size(plot_values, 2),'MinPeakProminence',param_peak);


signal = zeros(length(pks), 1);

for i = 1 : length(pks)
    half_width = floor(width/2);
    signal(i) = squeeze(mean(plot_values(locs(i)-half_width:locs(i)+half_width)));
end

% artery_counter = 0;
% for ii=2:jj
%     if (circle_path_mask(ii) - circle_path_mask(ii - 1)) == -1
%         artery_counter = artery_counter + 1;
%     elseif (circle_path_mask(ii) - circle_path_mask(ii - 1)) == 1
% 
%     else
%         1;
%     end
% end

%FIXME : plot_values another name
% plot_values = squeeze(mean(circle_path_signal, 2));
% plot_values = movmean(plot_values, 4);
% [pks , locs , width , ~] = findpeaks(plot_values,1:size(plot_values, 2),'MinPeakProminence',param_peak);

%FIXME : filtrage pics secondaires
% filtered_pks = pks;
% filtered_locs = locs;
% filtered_width = width;
% for ii=1:size(width)
%     if width(ii)<3 % 3 pixels threshold to delete spurious peaks
%         filtered_width(ii) = [];
%         filtered_locs(ii) = [];
%         filtered_pks(ii) = [];
%     end
% end
% pks = filtered_pks;
% locs = filtered_locs;
% width = filtered_width;

% figure(102)
% h = pcolor(circle_path_signal);
% set(h, 'EdgeColor', 'none');
% yticks(locs);
% yticklabels(1:size(locs));
% title("Luminosity evolution");
% pbaspect([1.618 1 1]);
% axis tight;

end