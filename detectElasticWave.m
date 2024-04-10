function detectElasticWave(video, mask_artery, maskCRA)
% choose a number of dominant frequencies to analyse
n_freq = 2;
pixel_size = 6/1535;
video_arteries = video .* mask_artery;

%% find dominant frequencies
spectrum_arteries = fft(video_arteries, [], 3);
fft_arteries_pos = spectrum_arteries(:,:,1:floor(size(video,3)/2));
spectrum_arteries_pos = squeeze(mean(abs(fft_arteries_pos), [1 2]));

[pks,locs] = findpeaks(spectrum_arteries_pos);
[~, I] = sort(pks, 'descend');
locs = locs(I);
% dominant_freq = locs(1:n_freq);

%% find center CRA
blurred_mask = imgaussfilt(double(mean(video,3).*double(maskCRA)),round(size(maskCRA,1)/4),'Padding',0);
[~,x_center] = findpeaks(sum(blurred_mask,1));
[~,y_center] = findpeaks(sum(blurred_mask,2));

% video = circshift(video, -(ceil(size(video, 1)/2) - x_center), 1);
% video = circshift(video, -(ceil(size(video, 2)/2) - y_center), 2);

%r = 100; % radius
%x = linspace(-r, r, size(video, 1)); % x values
x = 1:size(video, 1);
y = 1:size(video, 2);
%y = linspace(-r, r, size(video, 2)); % y values
[X, Y] = meshgrid(x, y);
D = sqrt((X - floor(length(x)/2)).^2 + (Y - floor(length(x)/2)).^2);

meanZ = zeros((floor(length(x)/2) - 1), n_freq);
ring_size = zeros((floor(length(x)/2) - 1), n_freq);


%% filter signal at given frequencies
%FIXME for now we take a given thickness of the band
% make sure dominant freq-5 is bigger than 0
for jj = 1 : n_freq
    dominant_freq = locs(jj);
    [a,b] = butter(10,[dominant_freq-3 dominant_freq+3]/floor(size(video,3)/2));
    [h,~] = freqz(b,a,floor(size(video,3)/2));
    H = repmat(h,1,size(video, 1), size(video, 2));
    H = permute(H, [2 3 1]);
    fft_arteries_pos_filt = fft_arteries_pos .* H;

    signal_artery = ifft(fft_arteries_pos_filt, [], 3);
%     signal_artery = ifft(fft_arteries_pos, [], 3);

    %% find correlations in filtered signal
    SFFT_signal_artery = (fft2(signal_artery));
    sspectrum_artery = squeeze(mean(abs(SFFT_signal_artery), 3));


    Z = (sspectrum_artery);
    %tmp = zeros(floor(length(x)/2)-1, size(Z,3));
    
    for i = 1 : floor(length(x)/2)-1
%         idx = find(D >= x(i), 1);
%         if ~isempty(idx)
            ring = (D >= x(i)) .* (D < x(i+1));
            if sum(ring, "all")>0
                %meanZ(i, jj) = squeeze(sum(Z.*ring, "all")/nnz(ring));
                ring_size(i,jj) = nnz(ring);
                meanZ(i, jj) = squeeze(sum(Z.*ring, "all"));
                %tmp(i,pp) = squeeze(sum(Z(:,:,pp).*ring, "all"));
            end
%         end
    end
end

[x, y] = convertk2lambda(meanZ(:,1), pixel_size);


figure(1)
plot(x, y);
1;




end