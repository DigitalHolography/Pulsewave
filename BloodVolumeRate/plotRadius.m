function [mean_BvrT, mean_std_BvrT] = plotRadius(vr_avg_r, vr_std_r, rad, index_start, index_end, name)

ToolBox = getGlobalToolBox;

numCircles = size(vr_avg_r, 2);
BvrR = sum(vr_avg_r, 2);
std_BvrR = sqrt(sum(vr_std_r .^ 2, 2)); % sqrt of the sum of variances

figure("Visible","off");
mean_BvrR = squeeze(mean(BvrR(:, :, index_start:index_end), 3))';
mean_std_BvrR = squeeze(rms(std_BvrR(:, :, index_start:index_end), 3))'; % quadratic mean
curve1 = mean_BvrR + 0.5 * mean_std_BvrR;
curve2 = mean_BvrR - 0.5 * mean_std_BvrR;
rad2 = [rad, fliplr(rad)];
inBetween = [curve1, fliplr(curve2)]';

fill(rad2, inBetween, Color_std);
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2);
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2);
plot(rad, mean_BvrR, '-k', 'LineWidth', 2);
yline(mean(mean_BvrR), '--k', 'LineWidth', 2);
legend({'', '', '', '', sprintf('mean = %0.2f µL/min', mean(mean_BvrR)), '', ''});

axis tight;
aa = axis;
aa(3) = -5;
aa(4) = 95;
axis(aa);
hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('radius in pixels')
title("Time average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'LineWidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_meanVolumeRate_%s_radius.png", ToolBox.main_foldername, name)))

figure("Visible","off");

hold on;

for circleIdx = 1:numCircles
    plot(fullTime, squeeze(BvrR(circleIdx, :, :)), 'LineWidth', 2);
end

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial variations of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_varianceVolumeRate_%s_time.png", ToolBox.main_foldername, name)))

figure("Visible","off");

mean_BvrT = squeeze(mean(BvrR, 1))';
mean_BvrT_value = mean(mean_BvrT(index_start:index_end));
max_BvrT_value = max(mean_BvrT(index_start:index_end));
mean_std_BvrT = squeeze(rms(std_BvrR, 1))';

hold off
curve1 = mean_BvrT + 0.5 * mean_std_BvrT;
curve2 = mean_BvrT - 0.5 * mean_std_BvrT;
ft2 = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)]';

fill(ft2, inBetween, Color_std);
hold on;
yline(0, 'k-', 'LineWidth', 2)
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, mean_BvrT, '-k', 'LineWidth', 2);
yline(mean_BvrT_value, '--k', 'LineWidth', 2)

plot(fullTime(index_start), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_end), 1.07 * max_BvrT_value, 'k|', 'MarkerSize', 10);
plot(fullTime(index_start:index_end), repmat(1.07 * max_BvrT_value, index_end - index_start + 1), '-k');
legend({'', '', '', '', '', sprintf('mean = %0.2f µL/min', mean_BvrT_value), '', ''});

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4) * 1.2])
box on

hold off

ylabel('Blood Volume Rate (µL/min)')
xlabel('time (s)')
title("Radial average of Blood Volume Rate")
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
set(gca, 'Linewidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_volumeRate_allrad_%s_time.png", ToolBox.main_foldername, name)))

end