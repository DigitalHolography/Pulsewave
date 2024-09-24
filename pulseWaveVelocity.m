function PWV = pulseWaveVelocity(U,mask,ToolBox,path)
% Computes the pulse wave velocity based on a cross correlation computation

% U(x,y,t) usually M0
% center the [x,y] barycenter (the center of the CRA)
PW_params = Parameters_json(path);

implay(rescale(U).*mask);

[N, M] = size(mask);
N_frame = size(U, 3);
[x, y] = meshgrid(1:M, 1:N);

timePeriod = ToolBox.stride / ToolBox.fs / 1000;

x_bary = ToolBox.x_barycentre;
y_bary = ToolBox.y_barycentre;

% radii approach
% m = floor((N+M)/2/10);
% 
% list_radius = linspace(0,1,m+1);
% U_r = zeros([m,N_frame]);
% 
% parfor i=1:m
%     c1 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i)* (N+M)/2;
%     c2 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i+1)* (N+M)/2;
%     U_r(i,:) = (sum(U.*(xor(c1,c2)&mask),[1,2])/nnz(xor(c1,c2)&mask));
% end
% 
% U_r = U_r(~isnan(U_r));
% U_r = reshape(U_r,[],N_frame);

dxx = 5;
skel = bwskel(mask);
grid = ones([M,N])<0;
grid(1:dxx:end,:)=true;
grid(:,1:dxx:end)=true;
interpoints = grid& skel; % get points interpolating with grid

numpoints = sum(interpoints,'all'); % get the number of points

[interpoints_y,interpoints_x] = ind2sub(size(interpoints),find(interpoints)); % y first
k = dsearchn([interpoints_x,interpoints_y],[x_bary,y_bary]); % get the nearest point to the center

absx = zeros(numpoints); % x and y position register
absy = zeros(numpoints);
abs_dist = zeros(numpoints); % vessel curvilign absis

absx(1) = interpoints_x(k); % nearest point to the center
absy(1) = interpoints_y(k);
 (k) = [];
interpoints_y(k) = [];
abs_dist(1) = 0;

for kb = 2:numpoints
    k = dsearchn([interpoints_x,interpoints_y],[absx(kb-1),absy(kb-1)]);
    absx(kb) = interpoints_x(k);
    absy(kb) = interpoints_y(k);
    interpoints_x(k) = []; % deleting point from list
    interpoints_y(k) = [];
end
for kb = 2:numpoints
    abs_dist(kb) = abs_dist(kb-1) + sqrt((absx(kb)-absx(kb-1))^2+(absy(kb)-absy(kb-1))^2)*PW_params.cropSection_pixelSize / 2 ^ PW_params.k;
end

figure(73)
plot(abs_dist);

L = single(zeros(size(mask)));
U_x = single(zeros([numpoints,N_frame]));
for i=1:numpoints
    sk_mask = ones(size(mask))<0;
    sk_mask(absy(i),absx(i)) = true;
    sectio = imdilate(sk_mask, strel('disk', floor(N * PW_params.masks_radius/5)))&mask;
    L(sectio) = i;
    U_x(i,:)=squeeze(mean(U.*sectio,[1,2]));
end



figure(74);
imagesc(L)
title('Selected sections along the artery')

figure(75);
imagesc(U_x)


ft_U_x = fft(U_x,[],2);
ph = angle(ft_U_x);
% xc = xcorr(U_x')'; % calculates all the time cross correlations between all the sections
% midpoint = round(numpoints/2);
% rr = ones(2*numpoints-1)<0;
% rr(numpoints:end)=true;
% for i=1:2*numpoints-1
%     xc_a = zeros(2*N_frame-1);
%     for i=1:numpoints
%         if rr(i)
%             xc_a = xc_a + xcorr(U_x(i,:),U_x(numpoints-i,:));
%         end
%         xc_a = xc_a/sum(rr(1:numpoints));
%     end
%     rr = circshift(rr,-i);
% 
%     xc_averaged(i,:) = mean(xc(rows,:),1); % averages all cross correlations between sections with an DX=(midpoint-i) distance between them
figure(99)
imagesc(U_x);
Ux = rescale(U_x,0,1,'InputMin',min(U_x,[],2),'InputMax',max(U_x,[],2));
figure(100)
imagesc(Ux);
% U_x = hilbert(U_x);
% xc = reshape(xcorr(U_x')',numpoints,numpoints,[]); % calculates all the time cross correlations between all the sections
% 
% rr = ones(3*numpoints-1);
% rr(2*numpoints-1) = 1;
% for i=1:2*numpoints-1
%     r = rr(numpoints:2*numpoints-1);
%     c = rr(1:numpoints);
%     toep = toeplitz(c,r);
%     rr = circshift(rr,-1);
%     xc_averaged(i,:) = mean(xc.*toep ,[1,2]);
% end
% figure(76);
% imagesc((xc_averaged));
% end
Ux = Ux - mean(Ux,2);
hU_x = hilbert(Ux')';
figure(101)
plot(Ux(50,:));hold on;
plot(real(hUx(50,:)));
plot(imag(hUx(50,:)));
title('rescaled and centered U(50,t) and its hilbert transform')
hxc = reshape(xcorr(hU_x')',numpoints,numpoints,[]); % calculates all the time cross correlations between all the sections

rr = ones(3*numpoints-1);
rr(2*numpoints-1) = 1;
for i=1:2*numpoints-1
    r = rr(numpoints:2*numpoints-1);
    c = rr(1:numpoints);
    toep = toeplitz(c,r);
    rr = circshift(rr,-1);
    hxc_averaged(i,:) = mean(real(hxc).*toep ,[1,2]);
end
figure(76);
imagesc((hxc_averaged));
end

