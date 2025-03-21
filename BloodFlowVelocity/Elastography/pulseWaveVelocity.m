function PWV = pulseWaveVelocity(U, mask)
% Computes the pulse wave velocity based on a cross correlation computation
% U is the field over which we compute the velocity and mask is the mask of
% the selected retinal artery

% U(x,y,t) usually M0
% center the [x,y] barycenter (the center of the CRA)
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

implay(rescale(U) .* mask);

[numX, numY] = size(mask);
N_frame = size(U, 3);

x_bary = ToolBox.x_barycenter;
y_bary = ToolBox.y_barycenter;

% radii approach
% m = floor((numX+numY)/2/10);
%
% list_radius = linspace(0,1,m+1);
% U_r = zeros([m,N_frame]);
%
% parfor i=1:m
%     c1 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i)* (numX+numY)/2;
%     c2 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i+1)* (numX+numY)/2;
%     U_r(i,:) = (sum(U.*(xor(c1,c2)&mask),[1,2])/nnz(xor(c1,c2)&mask));
% end
%
% U_r = U_r(~isnan(U_r));
% U_r = reshape(U_r,[],N_frame);

%% create a grid of points to select points along the skeletton of the artery mask
%
dxx = 5;
skel = bwskel(mask);
grid = ones([numY, numX]) < 0;
grid(1:dxx:end, :) = true;
grid(:, 1:dxx:end) = true;
interpoints = grid & skel; % get points interpolating with grid

numpoints = sum(interpoints, 'all'); % get the number of points

%% register the positions of points stating by the one closest to the CRA then going from closest to closest

[interpoints_y, interpoints_x] = ind2sub(size(interpoints), find(interpoints)); % y first
k = dsearchn([interpoints_x, interpoints_y], [x_bary, y_bary]); % get the nearest point to the center

absx = zeros([1, numpoints]); % x and y position register
absy = zeros([1, numpoints]);
abs_dist = zeros([1, numpoints]); % vessel curvilign absis

absx(1) = interpoints_x(k); % nearest point to the center
absy(1) = interpoints_y(k);
interpoints_x(k) = [];
interpoints_y(k) = [];
abs_dist(1) = 0;

for kb = 2:numpoints
    k = dsearchn([interpoints_x, interpoints_y], [absx(kb - 1), absy(kb - 1)]);
    absx(kb) = interpoints_x(k);
    absy(kb) = interpoints_y(k);
    interpoints_x(k) = []; % deleting point from list
    interpoints_y(k) = [];
end

for kb = 2:numpoints
    abs_dist(kb) = abs_dist(kb - 1) + sqrt((absx(kb) - absx(kb - 1)) ^ 2 + (absy(kb) - absy(kb - 1)) ^ 2) * params.px_size;
end

figure(73)
plot(abs_dist);

L = single(zeros(size(mask)));
U_x = single(zeros([numpoints, N_frame]));

for i = 1:numpoints
    sk_mask = ones(size(mask)) < 0;
    sk_mask(absy(i), absx(i)) = true;
    sectio = imdilate(sk_mask, strel('disk', floor(numX * params.masks_radius / 5))) & mask;
    L(sectio) = i;
    U_x(i, :) = squeeze(mean(U .* sectio, [1, 2]));
end

figure(74);
imagesc(L)
title('selected sections along the artery')

figure(75);
imagesc(U_x)

Ux = filloutliers(U_x', "linear", "movmedian", 20, ThresholdFactor = 1)';
figure(76);
imagesc(Ux)

Ux = rescale(Ux, 0, 1, 'InputMin', min(Ux, [], 2), 'InputMax', max(Ux, [], 2));
figure(77); imagesc(Ux);

% Ux = rescale(U_x,0,1,'InputMin',min(Ux,[],2),'InputMax',max(Ux,[],2));
% figure(100);imagesc(Ux);

ft_Ux = fft(Ux, [], 2);
ph = angle(ft_Ux);
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

% Ux = rescale(Ux,0,1,'InputMin',min(U_x,[],1),'InputMax',max(U_x,[],1));
% figure(100)
% imagesc(Ux);
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
%     xc_averaged(i,:) = mean(xc.*toep ,[1,2]);1
% end
% figure(76);
% imagesc((xc_averaged));
% end
Ux = Ux - mean(Ux, 2);
hUx = hilbert(Ux')';
figure(101)
plot(Ux(50, :)); hold on;
plot(real(hUx(50, :)));
plot(imag(hUx(50, :)));
title('rescaled and centered U(50,t) and its hilbert transform')

hxc = reshape(real(xcorr(exp(1j * angle(hUx')))'), [], numpoints, numpoints); % calculates all the time cross correlations between all the sections

imagesc(mean(squeeze(hxc(:, :, :)), 3))

% hxc = permute(hxc,[2,3,1]);
% rr = ones([1 3*numpoints-1]);
% rr(2*numpoints-1) = 1;
% for i=1:2*numpoints+1
%     r = rr(numpoints:2*numpoints-1);
%     c = rr(1:numpoints);
%     toep = toeplitz(c,r);
%     rr = circshift(rr,-1);
%     hxc_averaged(i,:) = mean((hxc).*toep ,[1,2]);
% end
% figure(76);
% imagesc((hxc_averaged));
end
