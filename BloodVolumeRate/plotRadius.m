function [mean_BvrT, mean_std_BvrT] = plotRadius(vr_avg_r, vr_std_r, fullTime, rad, index_start, index_end, name)

ToolBox = getGlobalToolBox;

numCircles = size(vr_avg_r, 1);
Color_std = [0.7 0.7 0.7];

BvrR = sum(vr_avg_r, 2); % sum of sections for each circle
std_BvrR = sqrt(sum(vr_std_r .^ 2, 2)); % sqrt of the sum of variances

figure("Visible","off");
mean_BvrR = squeeze(mean(BvrR(:, :, index_start:index_end), 3))';
mean_std_BvrR = squeeze(rms(std_BvrR(:, :, index_start:index_end), 3))'; % quadratic mean
curve1 = mean_BvrR + 0.5 * mean_std_BvrR;
curve2 = mean_BvrR - 0.5 * mean_std_BvrR;
rad2 = [rad, fliplr(rad)];
inBetween = [curve1, fliplr(curve2)];

fill(rad2, inBetween, Color_std);
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
plot(rad, mean_BvrR, '-k', 'LineWidth', 2);
yline(mean(mean_BvrR), '--k', 'LineWidth', 2);
legend({'', '', '', '', sprintf('mean = %0.2f µL/min', mean(mean_BvrR)), '', ''});

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0 , 1.07 * axP(4)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('radius in pixels')
title("Time average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'LineWidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_mean_%s_radius.png", ToolBox.main_foldername, name)))

figure("Visible","off");

hold on;

for circleIdx = 1:numCircles
    plot(fullTime, squeeze(BvrR(circleIdx, :, :)), 'LineWidth', 2);
end

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3) , 1.07 * axP(4)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial variations of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_variance_%s_time.png", ToolBox.main_foldername, name)))

figure("Visible", "off");

mean_BvrT = squeeze(mean(BvrR, 1))'; % mean of all circle's sum of all sections
mean_BvrT_value = mean(mean_BvrT(index_start:index_end)); % average in time of this signal
max_BvrT_value = max(mean_BvrT(index_start:index_end));
mean_std_BvrT = squeeze(rms(std_BvrR, 1))';
mean_std_BvrT_value = rms(mean_BvrT(index_start:index_end));

hold off
curve1 = mean_BvrT + 0.5 * mean_std_BvrT;
curve2 = mean_BvrT - 0.5 * mean_std_BvrT;
ft2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

hold on;
fill(ft2, inBetween, Color_std);
yline(0, 'k-', 'LineWidth', 2)
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, mean_BvrT, '-k', 'LineWidth', 2);
yline(mean_BvrT_value, '--k', 'LineWidth', 2)

plot(fullTime(index_start), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(index_end), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(index_start:index_end), repmat(1.07 * max_BvrT_value, index_end - index_start + 1), '-k', 'LineWidth', 2);

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3) , 1.07 * axP(4)])
box on

hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title(sprintf("Total blood volume rate (Avg. %0.2f µL/min)", mean_BvrT_value))
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_allrad_%s_time.png", ToolBox.main_foldername, name)))

fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername, '_', 'PW_main_outputs', '.txt')), 'a');
fprintf(fileID, 'Mean Blood Volume Rate %s : %f (µL/min) \r\n', name,mean_BvrT_value);
fprintf(fileID, 'Std Blood Volume Rate %s : %f (µL/min) \r\n', name,mean_std_BvrT_value);
fclose(fileID);

end