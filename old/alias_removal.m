function [img_svded] = alias_removal(img)
%FIXME : hanning window
%FIXME : use zeros in the FFT domain near 0 frequency + add butterworth filter (butter on matlab)

% img is a gray image for now (2D-matrix)

img = gpuArray(img);
img = single(img);
img_org = img;
img = img - mean(img,"all");

%% FFT section
fft_img = fftshift(fft2(img));
figure(1)
imagesc(log10(abs(fft_img)));
fft_img_filtered = zeros(size(fft_img));
fft_img_filtered(100:end-100,100:end-100) = fft_img(100:end-100,100:end-100);
%ring test
border1 = 150;
border2 = 200;
mask = ones(size(fft_img));
mask(border1:end-border1,border1:end-border1) = zeros(size(fft_img,1)-2*border1+1);
mask(border2:end-border2,border2:end-border2) = ones(size(fft_img,1)-2*border2+1);
fft_img_filtered_ring = fft_img.*mask;
figure(2)
imagesc(mask)


figure(3)
imagesc(log10(abs(fft_img_filtered_ring)));

img_filtered = abs(ifft2(fftshift(fft_img_filtered_ring)));
figure(4)
imagesc(mat2gray(img_filtered));
axis off
axis image
colormap gray
colorbar

figure(5) % original img
imagesc(mat2gray(img_org));
axis off
axis image
colormap gray
colorbar

% %% FIND LOCATION OF ALIASES
% c = xcorr2(img);
% figure(1)
% imagesc(c);
% figure(2)
% c_lineshape = mean(c(1022:1024,:),1);
% plot(1:size(c,2),c_lineshape);
% %FIXME : fix the parameters in the findpeaks
% %FIXME : erease the peaks near the center of the lineshape
% [peaks,locs] = findpeaks(c_lineshape);
% disp(locs);
% disp(peaks);
% locs = [locs(4),locs(8)];
% pks = [peaks(4),peaks(8)];
% 
% %% SVD FILTERING
% reshape(img_org,[size(img_org,1)*size(img_org,2),1]);
% A = zeros(size(img_org,1)*size(img_org,2),3);
% for ii=0
%     tmp_img_0 = circshift(img_org,[0,ii]);
%     tmp_img_1 = circshift(img_org,[0,locs(1)+ii]);
%     tmp_img_2 = circshift(img_org,[0,locs(2)+ii]);
%     A(:,1) = reshape(tmp_img_0,[size(tmp_img_0,1)*size(tmp_img_0,2),1]);
%     A(:,2) = reshape(tmp_img_1,[size(tmp_img_1,1)*size(tmp_img_1,2),1]);
%     A(:,3) = reshape(tmp_img_2,[size(tmp_img_2,1)*size(tmp_img_2,2),1]);
% end
% [U,S,V] = svd(A,"econ");
% size(S);
% Threshold1 = 2;
% Threshold2 = 3;
% img_svded = U*S(:,Threshold1:Threshold2)*V(:,Threshold1:Threshold2)';
% img_svded = reshape(img_svded,[size(img_org,1),size(img_org,2),size(S,1)]);
% 
% figure(1)
% for ii=1:3
%     imagesc(img_svded(:,:,mod(ii,3)+1));
%     axis off
%     axis image
%     colormap gray
%     colorbar 
%     pause(0.3)
% end
% 
% cor_factor = sqrt(mean(pks)/max(peaks));
% B = img_org-cor_factor*img_svded(:,:,2)-cor_factor*img_svded(:,:,3);
% 
% figure(3) % corrected img
% imagesc(mat2gray(B));
% axis off
% axis image
% colormap gray
% colorbar
% 
% figure(4) % original img
% imagesc(mat2gray(img_org));
% axis off
% axis image
% colormap gray
% colorbar

end
