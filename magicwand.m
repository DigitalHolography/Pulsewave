function mask = magicwand(im, tolerance, conn)
% ENTRY :
% im = maskArtery (binarized)
% conn = 8 by default (fonction bwlabel)
A = logical(im);
B = double(bwlabel(A,conn));
figure(666)
imagesc(B)
colorbar
maxiB = max(B(:));
B = B/maxiB;
B = B>tolerance;
mask1 = xor(A,B);

AA = flip(A,2);
BB = double(bwlabel(AA,conn));
figure(777)
imagesc(BB)
colorbar
maxiBB = max(BB(:));
BB = BB/maxiBB;
BB = BB>tolerance;
mask2 = xor(AA,BB);

mask = mask1.*flip(mask2,2);
figure(888)
imagesc(mask)
colorbar

end