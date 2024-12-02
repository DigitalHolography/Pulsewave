function [spec_plot, delt_spec_plot] = showSpectrum(ToolBox, maskArtery, maskBackground, maskSection, SH)
% Shows the spectrum of the preview batch if any

fs = ToolBox.fs;
f1 = ToolBox.f1;
f2 = ToolBox.f2;

SH = abs(SH) .^ 2;

spectrum = squeeze(mean(SH .* maskSection, [1, 2])); % The full spectrum of the power doppler image

j_win = size(SH, 3);
fullfreq = fftshift([0:j_win / 2 - 1, -j_win / 2:-1] * fs * 1000 / j_win);

% fitting to a lorentizian
exclude = abs(fullfreq) < f1 * 1000;

arterySpectrum = squeeze(sum(SH .* maskArtery .* maskSection, [1, 2]) / nnz(maskArtery .* maskSection));
backgroundSpectrum = squeeze(sum(SH .* maskBackground .* maskSection, [1, 2]) / nnz(maskBackground .* maskSection));
deltaSpectrum = arterySpectrum - backgroundSpectrum;

spec_plot = figure(57);

x = fullfreq;
y = double(10 * log(fftshift(spectrum / sum(spectrum(fftshift(~exclude))))));

fitEqn = @(a, b, x) ...
    10 * log(1 ./ ((1 + abs(x / a) .^ b)) / sum(1 ./ ((1 + abs(fullfreq(~exclude) / a) .^ b))));
l = fittype(fitEqn);
f = fit(x', double(y), l, Exclude = exclude);

plot(x / 1000, double(y), 'k-', 'LineWidth', 2); hold on;
plot(x / 1000, feval(f, x), 'k--', 'LineWidth', 2); hold on;
xline(-f1, 'k--', 'LineWidth', 2)
xline(-f2, 'k--', 'LineWidth', 2)
xline(f1, 'k--', 'LineWidth', 2)
xline(f2, 'k--', 'LineWidth', 2)
hold off
legend('avg spectrum', ['fit model 1/(1+(x/a)^b)', 'a = ', num2str(f.a), ' b = ', num2str(f.b)]);
title('Spectrum');
fontsize(gca, 14, "points");
xlabel("Frequency (kHz)", 'FontSize', 14);
ylabel("S(f) (dB)", 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);

delt_spec_plot = figure(59);

ya = double(10 * log(fftshift(arterySpectrum / sum(arterySpectrum(fftshift(~exclude))))));
yb = double(10 * log(fftshift(backgroundSpectrum / sum(backgroundSpectrum(fftshift(~exclude))))));
yd = double(10 * log(fftshift(deltaSpectrum / sum(deltaSpectrum(fftshift(~exclude))))));

lorentzEqn = @(a, b, x) ...
    10 * log(1 ./ ((1 + abs(x / a) .^ b)) / sum(1 ./ ((1 + abs(fullfreq(~exclude) / a) .^ b))));
l = fittype(lorentzEqn);
fl = fit(x', double(yd), l, Exclude = exclude); % fit a model to the delta spectrum.

plot(x / 1000, double(ya), 'k-', 'LineWidth', 2); hold on;
plot(x / 1000, double(yb), 'k--', 'LineWidth', 2); hold on;
plot(x / 1000, double(yd), 'k:', 'LineWidth', 2); hold on;
plot(x / 1000, feval(fl, x), 'k.', 'LineWidth', 2); hold on;
xline(-f1, 'k--', 'LineWidth', 2)
xline(-f2, 'k--', 'LineWidth', 2)
xline(f1, 'k--', 'LineWidth', 2)
xline(f2, 'k--', 'LineWidth', 2)
hold off
legend('artery spectrum', 'background spectrum', 'delta spectrum', ['lorentz model 1/(1+(x/a)^2)', 'a = ', num2str(fl.a), ' b = ', num2str(fl.b)]);
title('Spectrum');
fontsize(gca, 14, "points");
xlabel("Frequency (kHz)", 'FontSize', 14);
ylabel("S(f) (dB)", 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2);

% figure(58)
% spectrum_angle = squeeze(sum(SH_angle .* circle, [1, 2]) / nnz(circle));

% plot(fullfreq / 1000, 180 / pi * fftshift(spectrum_angle), 'k-', 'LineWidth', 2)
% hold on
% xline(time_transform.f1, 'k--', 'LineWidth', 2)
% xline(time_transform.f2, 'k--', 'LineWidth', 2)
% xline(-time_transform.f1, 'k--', 'LineWidth', 2)
% xline(-time_transform.f2, 'k--', 'LineWidth', 2)
% hold off
% title('Spectrum');
% fontsize(gca, 14, "points");
% xlabel("Frequency (kHz)", 'FontSize', 14);
% ylabel("arg(S(f)) (°)", 'FontSize', 14);
% pbaspect([1.618 1 1]);
% set(gca, 'LineWidth', 2);
% axis tight;

end
