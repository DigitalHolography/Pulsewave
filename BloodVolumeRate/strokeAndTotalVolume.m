function strokeAndTotalVolume(mean_BvrT, mean_std_BvrT, systolesIndexes, fullTime, numInterp)

ToolBox = getGlobalToolBox;

figure("Visible","off");

[interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(mean_BvrT, systolesIndexes, numInterp, mean_std_BvrT);
dt = (fullTime(2) - fullTime(1));
pulseTime = dt * (1:numInterp) * avgLength / numInterp;

[~, amin] = min(interp_BvrT);
[~, amax] = max(interp_BvrT);
cshiftn = numInterp - amin;

hold off

% Retinal Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime, fliplr(pulseTime)];
inBetween = [curve1, fliplr(curve2)]';
cRose = [254, 191, 210] / 255;
fill(ft2, inBetween, cRose, 'EdgeColor', 'none');
xline(pulseTime(end), 'k--', 'LineWidth', 2)

% Remaining Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
curve1 = curve1(1:min(amax + cshiftn,numInterp));
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime(1:min(amax + cshiftn,numInterp)), fliplr(pulseTime(1:min(amax + cshiftn,numInterp)))];
inBetween = [curve1, fliplr(curve2)]';
cCrimson = [222, 49, 99] / 255;
xline(pulseTime(min(amax + cshiftn,numInterp)), 'k--', 'LineWidth', 2)
fill(ft2, inBetween, cCrimson, 'EdgeColor', 'none');

% Grey STD and Signal
interp_BvrT2 = repmat(interp_BvrT, 1, 3);
interp_std_BvrT2 = repmat(interp_std_BvrT, 1, 3);
pulseTime2 = dt * (-numInterp + 1:numInterp * 2) * avgLength / numInterp;

hold on
curve1 = circshift(interp_BvrT2, cshiftn) + 0.5 * circshift(interp_std_BvrT2, cshiftn);
curve2 = circshift(interp_BvrT2, cshiftn) - 0.5 * circshift(interp_std_BvrT2, cshiftn);
ft2 = [pulseTime2, fliplr(pulseTime2)];
inBetween = [curve1, fliplr(curve2)]';
cSTD = [0.7 0.7 0.7];
fill(ft2, inBetween, cSTD, 'EdgeColor', 'none');
plot(pulseTime2, circshift(interp_BvrT2, cshiftn), '-k', 'LineWidth', 2);

yline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4)])
xlim([pulseTime(1)-1/2 * pulseTime(end), 3/2 * pulseTime(end)])
box on

ylabel('Blood Volume Rate (µL/min)')
xlabel('Time (s)')
ccinterpBvrT = circshift(interp_BvrT, cshiftn);
dt2 = pulseTime2(2) - pulseTime2(1);
stroke_volume_value = sum(ccinterpBvrT(1:min(amax + cshiftn,numInterp))) * dt2 / 60 * 1000; % in nL
total_volume_value = sum(ccinterpBvrT) * dt2 / 60 * 1000;
title(sprintf("Retinal Stroke Volume : %02.0f nL and Total Volume : %02.0f nL", stroke_volume_value, total_volume_value));
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
box on
set(gca, 'LineWidth', 2)

exportgraphics(gca, fullfile(ToolBox.PW_path_png, 'volumeRate', sprintf("%s_%s", ToolBox.main_foldername, 'strokeAndTotalVolume.png')))
end