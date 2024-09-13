function VOP = pulse_wave_velocity(U,mask,center)
% Computes the pulse wave velocity based on a cross correlation computation

% U(x,y,t) usually M0
% center the [x,y] barycenter (the center of the CRA)
[N, M] = size(mask);
N_frame = size(U, 3);
[x, y] = meshgrid(1:M, 1:N);

x_bary = center(1);
y_bary = center(2);

m = floor((N+M)/2/10);

list_radius = linspace(0,1,m+1);
U_r = zeros([m,N_frame]);

parfor i=1:m
    c1 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i)* (N+M)/2;
    c2 = sqrt((x-x_bary).^2+(y-y_bary).^2)<list_radius(i+1)* (N+M)/2;
    U_r(i,:) = (sum(U.*(xor(c1,c2)&mask),[1,2])/nnz(xor(c1,c2)&mask));
end

U_r = U_r(~isnan(U_r));
U_r = reshape(U_r,[],N_frame);

xcorr2(U_r)

end