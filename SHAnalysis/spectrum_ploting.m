function spectrum_ploting(SH, mask_artery, mask_background, fs, f1, f2)
%SPECTRUM_PLOTING Summary of this function goes here
%   Detailed explanation goes here

SH_artery = SH .* mask_artery;
SH_background = SH .* mask_background;

% %Calcul AVG background spectrum
spectrumAVG_background = squeeze(sum(SH_background, [1 2])) / nnz(SH_background(:, :, 1));
% spectrumAVG_background_calcul_omega = sum(SH_background,[1 2])/nnz(SH_background(:,:,1));
% momentM2_background = moment2(spectrumAVG_background_calcul_omega, f1, f2, fs, size(SH_background,3), 0);
% momentM0_background = moment0(spectrumAVG_background_calcul_omega, f1, f2, fs, size(SH_background,3), 0);
% omegaRMS_background=squeeze(sqrt(momentM2_background./momentM0_background));
% omegaRMS_index_background=omegaRMS_background*size(SH_artery,3)/fs;
% I_omega_background=log(spectrumAVG_background(round(omegaRMS_index_background)));
%
% %Calcul of Omega AVG artery
spectrumAVG_artery = squeeze(sum(SH_artery, [1 2])) / nnz(SH_artery(:, :, 1));
% spectrumAVG_artery_calcul_omega = sum(SH_artery,[1 2])/nnz(SH_artery(:,:,1));
% momentM2 = moment2(spectrumAVG_artery_calcul_omega, f1, f2, fs, size(SH_artery,3), 0);
% momentM0 = moment0(spectrumAVG_artery_calcul_omega, f1, f2, fs, size(SH_artery,3), 0);
% omegaRMS=squeeze(sqrt(momentM2./momentM0));
% omegaRMS_index=omegaRMS*size(SH_artery,3)/fs;
% I_omega=log(spectrumAVG_artery(round(omegaRMS_index)));

momentM2 = moment2(SH, f1, f2, fs, size(SH_artery, 3), 0);
momentM0 = moment0(SH, f1, f2, fs, size(SH_artery, 3), 0);
momentM2M0 = sqrt(momentM2 ./ mean(momentM0, [1 2]));
momentM2M0 = ff_correction(momentM2M0, round(0.15 * size(SH, 1)));

%Calcul AVG background spectrum
momentM1 = moment1(SH, f1, f2, fs, size(SH_background, 3), 0);
[~, v] = min(imgaussfilt(momentM1, 10), [], 'all');
[y, x] = ind2sub([size(SH, 1) size(SH, 2)], v);
[X, Y] = meshgrid(1:size(SH, 1), 1:size(SH, 2));
circle_mask1 = sqrt((X - x) .^ 2 + (Y - y) .^ 2) <= size(SH, 2) / 2;
circle_mask1 = ~(circle_mask1)';
mask_background = mask_background .* circle_mask1;
% SH_background = SH_background .* circle_mask1;
% momentM2_background = moment2(SH_background, f1, f2, fs, size(SH_background,3), 0);
% momentM0_background = moment0(SH_background, f1, f2, fs, size(SH_background,3), 0);
% momentM0_background(momentM0_background==0) = 1;
% momentM2M0_background = sqrt(momentM2_background ./ momentM0_background);
% momentM2M0_background = flat_field_correction(momentM2M0_background, 153);
momentM2M0_background = momentM2M0 .* mask_background;
omegaRMS_background = sum(momentM2M0_background, [1 2]) / nnz(SH_background(:, :, 1));
omegaRMS_index_background = omegaRMS_background * size(SH_artery, 3) / fs;
I_omega_background = log(spectrumAVG_background(round(omegaRMS_index_background)));

%Calcul of Omega AVG artery
% momentM2_artery = moment2(SH_artery, f1, f2, fs, size(SH_artery,3), 0);
% momentM0_artery = moment0(SH_artery, f1, f2, fs, size(SH_artery,3), 0);
% momentM0_artery(momentM0_artery==0) = 1;
% momentM2M0_artery = sqrt(momentM2_artery ./ momentM0_artery);
% momentM2M0_artery = flat_field_correction(momentM2M0_artery, 153);
momentM2M0_artery = momentM2M0 .* mask_artery;
omegaRMS = sum(momentM2M0_artery, [1 2]) / nnz(SH_artery(:, :, 1));
omegaRMS_index = omegaRMS * size(SH_artery, 3) / fs;
I_omega = log(spectrumAVG_artery(round(omegaRMS_index)));

disp(omegaRMS);
disp(omegaRMS_background);

axis_x = linspace(-fs / 2, fs / 2, size(SH_artery, 3));

% omegaAVG_left = [0.4 0.6];
% omegaAVG_right = [0.4 0.4];
figure(33533)

p_artery = plot(axis_x, fftshift(log(spectrumAVG_artery)), 'red', 'LineWidth', 1, 'DisplayName', 'Arteries');
ylim([.99 * log(min(spectrumAVG_artery)) .7 * log(max(spectrumAVG_artery))])
xlim([-fs / 2 fs / 2])
rectangle('Position', [-f1 .95 * log(min(spectrumAVG_artery)) 2 * f1 (.7 * log(max(spectrumAVG_artery)) - .95 * log(min(spectrumAVG_artery)))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [-fs / 2 .95 * log(min(spectrumAVG_artery)) (fs / 2 - f2) (.7 * log(max(spectrumAVG_artery)) - .95 * log(min(spectrumAVG_artery)))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [f2 .95 * log(min(spectrumAVG_artery)) (fs / 2 - f2) (.7 * log(max(spectrumAVG_artery)) - .95 * log(min(spectrumAVG_artery)))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
om_RMS_line = line([-omegaRMS omegaRMS], [I_omega I_omega]);
om_RMS_line.Color = 'red';
om_RMS_line.LineStyle = '-';
om_RMS_line.Marker = '|';
om_RMS_line.MarkerSize = 12;
om_RMS_line.LineWidth = 1;
om_RMS_line.Tag = 'f RMS';
text(0, I_omega, 'f_{RMS}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

xline(f1, '--')
xline(f2, '--')
xline(-f1, '--')
xline(-f2, '--')
xticks([-f2 -f1 0 f1 f2])
xticklabels({num2str(round(-f2, 1)), num2str(round(-f1, 1)), '0', num2str(round(f1, 1)), num2str(round(f2, 1))})
% annotation('textarrow',omegaAVG_left,omegaAVG_right,'String','omega','FontSize',13,'Linewidth',2)
title('Average spectrum')

% plot(fullTime,fullArterialPulse,'-k', fullTime,fullBackgroundSignal,':k', fullTime, fullVenousSignal, '-.k', 'LineWidth',2) ;
% title('arterial pulse waveform and background signal'); % averaged outside of segmented vessels
fontsize(gca, 12, "points");
xlabel('frequency (kHz)', 'FontSize', 14);
ylabel('log S', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 1);
uistack(p_artery, 'top');
uistack(gca, 'top');

hold on
p_background = plot(axis_x, fftshift(log(spectrumAVG_background)), 'black--', 'LineWidth', 1, 'DisplayName', 'Background');
om_RMS_line = line([-omegaRMS_background omegaRMS_background], [I_omega_background I_omega_background]);
om_RMS_line.Color = 'black';
om_RMS_line.LineStyle = '-';
om_RMS_line.Marker = '|';
om_RMS_line.MarkerSize = 12;
om_RMS_line.LineWidth = 1;
om_RMS_line.Tag = 'f RMS';
text(0, I_omega_background, 'f_{RMS background}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
legend('', '', '', '', '', 'Arteries', 'Background');

% print('-f33533','-depsc', 'C:\Users\Vladikavkaz\Pictures\article_Pulse\spectrum.eps') ;
% print('-f33533','-dpng', 'C:\Users\Vladikavkaz\Pictures\article_Pulse\spectrum.png') ;

% x = [0.74 0.79];    % adjust length and location of arrow
% y = [0.3 0.3];      % adjust hieght and width of arrow
% annotation('textarrow',x,y,'String',' Growth ','FontSize',13,'Linewidth',2)

end
